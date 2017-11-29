 #include "Penelope.h"
#include "SourceHead.h"
#include "egs_input.h"
#include "egs_base_geometry.h"
#include "egs_functions.h"

#define ti tp[it]
#define SOFT false;

// Begin: essential for command line call
bool b_abort = false;
bool b_thread_active = false;

void executeJob(const char* configFileName, MPS& configMacro);

extern "C" __declspec(dllexport) void startSimulation(const char* configFileName, MPS& configMacro, bool bWait) //launch a thread to do the simulation
{
	b_thread_active = true;
	std::thread thd(executeJob, configFileName, std::ref(configMacro));
	if (bWait) thd.join();
	else thd.detach();
}

extern "C" __declspec(dllexport) void stopSimulation()
{
	if (b_thread_active) b_abort = true; //only if the thread is on
	b_thread_active = false;
}
// End: essential for command line call

// This function will allow larger turning angle than the segment approximation
void moveEP(ThreadParameters& p, double ds, double Rb, double t)
{
	double arg = ds / Rb; // the rotation angle

	double f1, f2;
	if (arg < 0.2) // use power series to approximate f1 and f2 when arg is small
	{
		double arg2 = arg*arg;
		f1 = -0.5*arg2 + 1.0 / 24.0*arg2*arg2;  // for 0.2, relative error is 2.2e-6
		f2 = arg - 1.0 / 6.0*arg*arg2;          // for 0.2, relative error is 1.3e-5, absolute error is 2.6e-6
	}
	else
	{
		f1 = cos(arg) - 1;
		f2 = sin(arg);
	}

	p.x.x += Rb*(-f1*t*p.v.y + f2*p.v.x);
	p.x.y += Rb*(f1*t*p.v.x + f2*p.v.y);
	p.x.z += ds*p.v.z;

	double vxt = p.v.x;
	p.v.x += f1*vxt - f2*t*p.v.y;
	p.v.y += f1*p.v.y + f2*t*vxt;
	//p.v.z += 0
}

void executeJob(const char* configFileName, MPS& configMacro) //execute one job according to the config file
{
	RunTimeCounter totalTime; // Counter the total running time

	//parse the whole config file
	ConfigFile cf(configFileName);

	//find out where is the config file located, and do macro replacement
	string cdir(configFileName);
	size_t pos = cdir.find_last_of("\\/");
	cdir.erase(++pos, string::npos);
	cf.macroReplace(string("$cdir$"), cdir);
	cf.macroReplace(configMacro);

	//initialize the log 
	string logDir, logAppend, logDescription;
	cf.getValue("log file directory", logDir);
	cf.getValue("log file append", logAppend);
	cf.getValue("log description", logDescription);
	if (0 == logAppend.compare("no")) logAppend = "w";
	else logAppend = "a";
	Log.setLogName(0, logAppend, logDir);//start log recording for this job
	Log("The job config file name = %s", configFileName);
	if (logDescription.compare("NA") != 0) Log("Short description: %s", logDescription.c_str());
	Log("Start log time = %s\n\n", Log.timeNow());

	//get the thread number
	int NThread = 1; // it is the logic CPU core number of this device
	if (!cf.getValue("NThread", NThread))  NThread = get_thread_num(); //the number will take over the default one

	//initialize the source
	ConfigFile *subCF = cf.getBlock("SOURCEHEAD");
	SourceHead_Init(subCF, NThread);
	double Emax = SourceHead_Emax();

	int AirID = -1;
	//initialize the phantom using EGS++ geometry library
	//egsWarningVisible();
	subCF = cf.getBlock("PHANTOM");
	string phtFile;
	subCF->getValue("geometry file", phtFile);
	EGS_Input* input = new EGS_Input;
	input->setContentFromFile(phtFile.c_str());
	EGS_Input *geom_input = input->takeInputItem("geometry definition");
	if (!geom_input) {
		exitApp("No geometry definition in the input file\n");
	}
	EGS_BaseGeometry* g = EGS_BaseGeometry::createGeometry(geom_input);
	if (!g) {
		exitApp("Failed to construct the simulation geometry\n");
	}
	delete geom_input;
	EGS_BaseGeometry::describeGeometries();
	egsInformation("\nThe simulation geometry is of type %s and has the "
		"name '%s'\n\n", g->getType().c_str(), g->getName().c_str());
	vector<int> ids; //extract the material ID (the name format should be ID_denotation)
	for (int j = 0; j < g->nMedia(); j++)
	{
		const char *medname = g->getMediumName(j);
		try {
			int penid = std::atoi(medname);
			if (penid == 104) AirID = j; // the material id of air in PENELOPE database is 104
			ids.push_back(penid);
		}
		catch (std::exception&) {
			egsFatal("Invalid material name format, i.e. ID_denotation\n");
		}
	}
	double Bz = 0;
	subCF->getValue("magnetic field", Bz); // z component of the magnetic field
	double rf = 100 / (lightSpeed*fabs(Bz)); //radius factor= 100/C/B

	//initialize PENELOPE
	subCF = cf.getBlock("PENELOPE");
	PENELOPE PE;
	PE.init(subCF, NThread, Emax, ids);
	Material* mat = PE.getMat();
	ThreadParameters* tp = PE.getTP();

	//initialize the array to score the dose
	double* AirDE = new double[NThread];
	for (int i = 0; i < NThread; ++i) AirDE[i] = 0;
	double dose = 0, uncertainty = 0;

	double fSIMU = 1;
	int NBatch = 50;
	cf.getValue("NSIMU", fSIMU);
	cf.getValue("NBatch", NBatch);
	double hBatch = fSIMU / NBatch;
	double fNSIMU = hBatch / NThread;
	/***************************************************************************************/
	
	const double A_BIT_FAR = 1e-3; // unit cm, add this distance whenever the particle crosses a boundary
	const double MIN_CURVE_DISTANCE = 10e-2; // unit cm, the minimum distance to the surface within one material. Current value = 100 micro meter
	const double ARG_MAX = 0.05; // the max angle that one line segment approximates
	
	Log("\nCalculating dose with %d thread(s), please wait patiently...\n\n", NThread);
	RunTimeCounter rc; // counting for the total simulation time
	RunTimeCounter percentSpeed;
	long loopCounter = 0; // for debug purpose only
	for (int ib = 0; ib < NBatch; ++ib)
	{
#ifdef USE_OPENMP // all thread works together to finish one batch
#pragma omp parallel num_threads(NThread)
#endif
		{
			// For efficiency concern, I will use goto statements. Sorry for messing up the logic
			int it = omp_get_thread_num(); //thread index starting from 0
			if (0 == it) percentSpeed.start();
			double hist = 0;
			for (; hist <= fNSIMU; ++hist)
			{
				for (int nSample = 1; nSample > 0;)
				{
					nSample = SourceHead_Sample(it, ti.E, ti.weight, ti.x, ti.v, ti.type);
					if (nSample < 0) break;
					/*************************** Begin: simulate this new particle ************************/

					//propagate the particle to the material by straight line
					double t = 1e30;
					int imed = -1; // Current material index, used in mat[imed] functions calls
					ti.RegionID = -1; // Current region index
					// advance the particle for distance t and return region and material index
					ti.RegionID = g->howfar(ti.RegionID, (const EGS_Vector &)ti.x, (const EGS_Vector &)ti.v, t, &imed);
					if (ti.RegionID < 0) continue; // The particle won't intersect with the target object
					t += A_BIT_FAR; // We move a little bit further to reduce possible rounding problems
					ti.x += ti.v*t; // Update the position

					++loopCounter; // for debug

					bool bNeedRefill = false;
					ti.Mode = SOFT; // Begin with the SOFT knock
					double DS = 0; // Jump distance
					int iNreg, iNmed; // The new region and material ID
					double Rb, Bzt; // used in MoveEP and line segment approximation

					while (true) // Loop to simulate one history
					{
						DS = mat[imed].jump(it); // jump some distance
						// g->hownear() seems less expensive than g->howfar() ???
						double perp = g->hownear(ti.RegionID, (const EGS_Vector &)ti.x);

						if (ti.type != photon)
						{
							Rb = rf*sqrt(ti.E*(ti.E + TEs));
							Bzt = (ti.type == positron) ? -1 : 1;
							if (Bz < 0) Bzt = -Bzt;
						}

						if (DS < perp) // It will remain in the same region
						{
							if (ti.type == photon) ti.x += ti.v*DS;
							else moveEP(ti, DS, Rb, Bzt);
						}
						else // Less frequent but more expensive call
						{
							if (ti.type == photon) // It's photon, so we will call g->howfar() to proceed
							{
HowfarMove:						t = DS; // Loop entery that need call g->howfar() directly
								iNreg = g->howfar(ti.RegionID, (const EGS_Vector &)ti.x, (const EGS_Vector &)ti.v, t, &iNmed);
								if (iNreg < 0) bNeedRefill = true; //  The current particle escaped the geometry object
								else
								{
									if (iNreg == ti.RegionID) ti.x += ti.v*DS; // Not crossing a boundary, move DS
									else // Crossed a region, either goto CaseNewMat or CaseSameMat
									{
										ti.RegionID = iNreg; // Update region id
										t += A_BIT_FAR; // we move a little bit further to reduce possible rounding problems
										ti.x += ti.v*t; // Update position: move distance t
										if (iNmed == imed) // Although crossing a boundary, the material is the same
										{
											DS -= t; // The distance to move must be deduced by t
											goto HowfarMove;
										}
										else // Crossing a boundary and entering a different material
										{
											imed = iNmed; // Update media id
#ifdef DEBUG
											if (imed < 0)
											{
												printf("imed<0 error!\n");
												exitApp("error1");
											}
#endif
											ti.Mode = SOFT;
											DS = mat[imed].jump(it);
											goto HowfarMove;
										}
									}
								}
							}
							else // It's electron/positron, we will call g->hownear() iterately to advance in curve
							{								
								bool bMoveStyle = true; // true, curve movement; false, line segment approximation
								moveEP(ti, perp, Rb, Bzt);
								DS -= perp;
								perp = g->hownear(ti.RegionID, (const EGS_Vector &)ti.x); // update the perpendicular distance

								while (true)
								{	
									// will do direct curve calculation
									if (bMoveStyle) 
									{
										moveEP(ti, perp, Rb, Bzt);
										DS -= perp;
										perp = g->hownear(ti.RegionID, (const EGS_Vector &)ti.x); // update the perpendicular distance

										if (perp > DS) // can finish the moving without crossing any boundary
										{
											moveEP(ti, DS, Rb, Bzt);
											break; 
										}
										if (perp < MIN_CURVE_DISTANCE) bMoveStyle = false; // will use line segment to approach the boundary
									}
									// use the line segment approximation
									else
									{										
										double s = Rb*ARG_MAX;
										t = s; // test if it will intersect with the boundary by calling g->howfar()
										iNreg = g->howfar(ti.RegionID, (const EGS_Vector &)ti.x, (const EGS_Vector &)ti.v, t, &iNmed);
										if (t >= DS) // exhaust DS in one step, then we can stop
										{
											ti.x += ti.v*DS;
											break;
										}
										// else need more step(s)
										if (iNreg < 0) //  The current particle escaped the geometry object
										{
											bNeedRefill = true;
											break;
										}
										else
										{
											if (iNreg == ti.RegionID) ti.x += ti.v*s; // Not crossing a boundary
											else // Crossed a region
											{
												ti.RegionID = iNreg; // Update region id
												if (iNmed == imed) // Although crossing a boundary, the materials are the same
												{
													ti.x += ti.v*s; // use up the s distance. Since s is usually very small,
																	// we think s won't cause to cross another boundary
												}
												else // Crossing a boundary and entering a different material
												{
													imed = iNmed; // Update media id
													ti.Mode = SOFT;
													t += A_BIT_FAR; // increase the moving distance a little bit to reduce rounding problems
													ti.x += ti.v*t; // Update position: move distance t
													DS = mat[imed].jump(it);
													continue; // this will keep on using line segment approximation
												}
											}
										}

										DS -= s; // DS is reduced by one segment s
										perp = g->hownear(ti.RegionID, (const EGS_Vector &)ti.x);
										if (perp > MIN_CURVE_DISTANCE) bMoveStyle = true;
										else // continue with the segment approximation
										{
											// change the direction slightly
											double f1 = cos(ARG_MAX) - 1; // This should be optimized by compiler automatically
											double f2 = sin(ARG_MAX);
											double vxt = ti.v.x;
											ti.v.x += f1*vxt - f2*Bzt*ti.v.y;
											ti.v.y += f1*ti.v.y + f2*Bzt*vxt;
										}
									}	
								}

							}
						}


						if (!bNeedRefill) // Knock to loss energy and change direction
						{
							double DE = mat[imed].knock(it);

							if (imed == AirID) AirDE[it] += DE; // Record the energy deposited into air cavity
							if (ti.E == 0) bNeedRefill = true;
						}

						if (bNeedRefill)
						{
							if (ti.secStack.empty()) break; // Break simulating one history
							ti.secStack.pop(ti.E, ti.weight, ti.RegionID, ti.x, ti.v, ti.type);
							imed = g->medium(ti.RegionID); // Get the material id for the stored particle

#ifdef DEBUG
							if (imed < 0)
							{
								printf("imed<0 error!\n");
								exitApp("error2");
							}
#endif
							bNeedRefill = false;
							ti.Mode = SOFT;
						}
					}
					/***************************  End: simulate this new particle  ************************/
				}
			}
			if (0 == it) Log("Thread 0: processed ------------------------ %3d%%,   speed = %d h/s\n",
				int((ib + 1) * 100 / (double)NBatch), int(hist / percentSpeed.stop(true))); //report progress
		}
		double DeltaE = 0;
		for (int i = 0; i < NThread; ++i)
		{
			DeltaE += AirDE[i];
			AirDE[i] = 0;
		}
		dose += DeltaE;
		uncertainty += DeltaE*DeltaE;
	}

	delete[] AirDE;

	Log("Thread 0: processed ------------------------ 100%%, speed = %d h/s\n", int(fSIMU / NThread / rc.stop(true)));
	Log("\nHistory number on this node is %.0f due to process division\n", fSIMU);
	Log("\nThe simulation costs %f seconds, ie %f minutes\n", rc.stop(true), rc.stop(true) / 60);

	/*****************************************************************************************************/
	string outname;
	cf.getValue("output file name", outname);
	Log("\nThe total run time from entering main() is %f s\n\n\n\n", totalTime.stop());

	Log("\nThe energy in air cavity is %g Gy*cm^3", dose*1.602e-16*SourceHead_NormFactor() / fSIMU);
	Log("The uncertainty in air cavity is %f", sqrt(uncertainty / (dose*dose) - 1.0 / NBatch));
	/*******************clean up work for this run***********************/
	SourceHead_Delete();
}
