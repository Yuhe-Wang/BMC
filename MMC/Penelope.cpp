
#define CallStyle //_stdcall //if use Compaq Visual Fortran compiled pinit.dll
#include "Penelope.h"

//#define STEP_DEBUG_STATUS
#define STEP_NUM 1700000
//#define ACCU_DEBUG
#define INCIDENT_NUM 1000000
/************************Start: common data, unchanged during simulation **************************/
double ELtab[NEGrid]; //log(E) grid table
double Etab[NEGrid];
double ELow, ELowL, EUp; //energy up and low bound
double DLE; //delta log(E)
double LDE; //linear delta E

double FLCOH;
double EPH[8000], XPH[8000][10];
int IPHF[99], IPHL[99], NPHS[99];
double ET[15000], FR[15000];
int IAL[15000], IS1[15000], IS2[15000], IFIRST[99][9], ILAST[99][9];
double XESI[6000][9], XPSI[6000][9];
int IESIF[99], NSESI[99], IPSIF[99], NSPSI[99];
double EB[99][30];
double BET[6];

ThreadParameters* tp = NULL;
SFloat *density = NULL;
SFloat maxDensity = 1;
/************************End: common data, unchanged during simulation **************************/

void calcG12(double b, double &g1, double &g2) //used in pair production
{
	double b2 = b*b;
	double barctan = 0;
	if (b > 1e-10) barctan = b*atan(1.0 / b);
	else barctan = b* PI / 2.0;
	double lnb = 2 * log(1 + b2);
	double A = b2*(4 - 4 * barctan - 3 * log(1.0 + 1.0 / b2));
	g1 = 7.0 / 3.0 - lnb - 6 * barctan - A;
	g2 = 11.0 / 6.0 - lnb - 3 * barctan + 0.5*A;
}

/************************Start: implemtents of class Material **************************/
void Material::load(int mati)
{


    M = mati;
    char name[20];
    sprintf(name, "Materials/mat%d.dat", M);

    FILE *fp = NULL;

    fp = fopen(name, "rb");                                                                                                                                                            
	if (fp == NULL) exitApp("Cannot open Material file!");
	fread(&M, sizeof(int), 1, fp);
	fread(&Eabs[electron], sizeof(double), 1, fp);
	fread(&Eabs[photon], sizeof(double), 1, fp);
	fread(&Eabs[positron], sizeof(double), 1, fp);
	fread(&C1, sizeof(double), 1, fp);
	fread(&C2, sizeof(double), 1, fp);
	fread(&Wcc, sizeof(double), 1, fp);
	fread(&Wcr, sizeof(double), 1, fp);
	fread(&Ecutr, sizeof(double), 1, fp);
	fread(&Zt, sizeof(double), 1, fp);
	fread(&At, sizeof(double), 1, fp);
	fread(&massDensity, sizeof(double), 1, fp);
	fread(&numDensity, sizeof(double), 1, fp);
	fread(&elemNum, sizeof(int), 1, fp);
	fread(weight, sizeof(double), 30, fp);
	fread(IZ, sizeof(int), 30, fp);
	fread(SGRA, sizeof(double), 200, fp);
	fread(SGCO, sizeof(double), 200, fp);
	fread(SGPH, sizeof(double), 200, fp);
	fread(SGPP, sizeof(double), 200, fp);
	fread(X2COH, sizeof(double), 241, fp);
	fread(PDCOH, sizeof(double), 241, fp);
	fread(FCO, sizeof(double), 64, fp);
	fread(UICO, sizeof(double), 64, fp);
	fread(FJ0, sizeof(double), 64, fp);
	fread(KZCO, sizeof(int), 64, fp);
	fread(KSCO, sizeof(int), 64, fp);
	fread(&NOSCCO, sizeof(int), 1, fp);
	fread(&ZEQPP, sizeof(double), 1, fp);
	fread(&Fpp, sizeof(double), 1, fp);
	fread(&Bpp, sizeof(double), 1, fp);
	fread(SEHEL, sizeof(double), 200, fp);
	fread(SEHIN, sizeof(double), 200, fp);
	fread(SEISI, sizeof(double), 200, fp);
	fread(SEHBR, sizeof(double), 200, fp);
	fread(SETOT, sizeof(double), 200, fp);
	fread(CSTPE, sizeof(double), 200, fp);
	fread(RSTPE, sizeof(double), 200, fp);
	fread(DEL, sizeof(double), 200, fp);
	fread(W1E, sizeof(double), 200, fp);
	fread(W2E, sizeof(double), 200, fp);
	fread(T1E, sizeof(double), 200, fp);
	fread(T2E, sizeof(double), 200, fp);
	fread(RNDCE, sizeof(double), 200, fp);
	fread(AE, sizeof(double), 200, fp);
	fread(BE, sizeof(double), 200, fp);

	fread(SPHEL, sizeof(double), 200, fp);
	fread(SPHIN, sizeof(double), 200, fp);
	fread(SPISI, sizeof(double), 200, fp);
	fread(SPHBR, sizeof(double), 200, fp);
	fread(SPAN, sizeof(double), 200, fp);
	fread(SPTOT, sizeof(double), 200, fp);
	fread(CSTPP, sizeof(double), 200, fp);
	fread(RSTPP, sizeof(double), 200, fp);
	fread(W1P, sizeof(double), 200, fp);
	fread(W2P, sizeof(double), 200, fp);
	fread(T1P, sizeof(double), 200, fp);
	fread(T2P, sizeof(double), 200, fp);
	fread(RNDCP, sizeof(double), 200, fp);
	fread(AP, sizeof(double), 200, fp);
	fread(BP, sizeof(double), 200, fp);

	fread(&EELMAX, sizeof(double), 1, fp);
	fread(&PELMAX, sizeof(double), 1, fp);

	fread(XSE, sizeof(double), 128 * 200, fp);
	fread(PSE, sizeof(double), 128 * 200, fp);
	fread(ASE, sizeof(double), 128 * 200, fp);
	fread(BSE, sizeof(double), 128 * 200, fp);
	fread(ITLE, sizeof(int), 128 * 200, fp);
	fread(ITUE, sizeof(int), 128 * 200, fp);

	fread(XSP, sizeof(double), 128 * 200, fp);
	fread(PSP, sizeof(double), 128 * 200, fp);
	fread(ASP, sizeof(double), 128 * 200, fp);
	fread(BSP, sizeof(double), 128 * 200, fp);
	fread(ITLP, sizeof(int), 128 * 200, fp);
	fread(ITUP, sizeof(int), 128 * 200, fp);

	fread(&meanExcitE, sizeof(double), 1, fp);
	fread(&OP2, sizeof(double), 1, fp);

	fread(F, sizeof(double), 64, fp);
	fread(Ui, sizeof(double), 64, fp);
	fread(WRI, sizeof(double), 64, fp);
	fread(KZ, sizeof(int), 64, fp);
	fread(KS, sizeof(int), 64, fp);

	fread(&NOSC, sizeof(int), 1, fp);

	for (int i = 0; i < 64; ++i)
	{
		for (int j = 0; j < 200; ++j) fread(&EINAC[j][i], sizeof(double), 1, fp);
	}
	for (int i = 0; i < 64; ++i)
	{
		for (int j = 0; j < 200; ++j) fread(&PINAC[j][i], sizeof(double), 1, fp);
	}

	fread(PBcut, sizeof(double), 200, fp);
	fread(WBcut, sizeof(double), 200, fp);

	for (int i = 0; i < 32; ++i)
	{
		for (int j = 0; j < 200; ++j) fread(&PDFB[j][i], sizeof(double), 1, fp);
	}
	for (int i = 0; i < 32; ++i)
	{
		for (int j = 0; j < 200; ++j) fread(&PACB[j][i], sizeof(double), 1, fp);
	}

	for (int i = 0; i < 6; ++i)
	{
		for (int j = 0; j < 21; ++j)
		{
			for (int k = 0; k < 4; ++k) fread(&BP1[i][j][k], sizeof(double), 1, fp);
		}
	}
	for (int i = 0; i < 6; ++i)
	{
		for (int j = 0; j < 21; ++j)
		{
			for (int k = 0; k < 4; ++k) fread(&BP2[i][j][k], sizeof(double), 1, fp);
		}
	}
	

	fclose(fp);
}

void Material::start()
{
	int it = omp_get_thread_num();//thread index
	Particle &p = tp[it].p;
	KnockMode &Mode = tp[it].Mode;
	if (p.E<1.00000001*ELow || p.E>0.99999999*EUp) exitApp("E out of range in function start()");
	Mode = SOFT;
}

double Material::jump(double DSMax, int it)
{
	//int it = omp_get_thread_num(); //thread index
	PENRNG &rng = tp[it].rng;
	Particle &p = tp[it].p;
	int &IE = tp[it].IE; //left index of interval
	double &XIE = tp[it].XIE; //portion in that interval
	//inverse mean free path related parameters
	double &p1 = tp[it].p1;
	double &p2 = tp[it].p2;
	double &p3 = tp[it].p3;
	double &p4 = tp[it].p4;
	double &p5 = tp[it].p5;
	double &p6 = tp[it].p6;
	double &st = tp[it].st;

	double &W1 = tp[it].W1;
	double &W2 = tp[it].W2;
	double &T1 = tp[it].T1;
	double &T2 = tp[it].T2;
	double &DS1 = tp[it].DS1;
	double &DST = tp[it].DST;
	KnockMode &Mode = tp[it].Mode;
	bool &KSoftI = tp[it].KSoftI;
	bool &KSoftE = tp[it].KSoftE;
	bool &KDelta = tp[it].KDelta;

	//DSMax, and return value both have unit cm
	//locate the energy interval first
	double XE = (log(p.E) - ELowL) / DLE;
	IE = (int)XE;
	XIE = XE - IE;
	
	//need to track the Positions of each particles.
	if (p.type == photon)
	{
		p1 = exp((1 - XIE)*SGRA[IE] + XIE*SGRA[IE + 1]); //N*sigma for Rayleigh scattering
		p2 = exp((1 - XIE)*SGCO[IE] + XIE*SGCO[IE + 1]); //N*sigma for Compton scattering
		p3 = SGPH[IE];
		if (p.E < 1.023e6) p4 = 0;
		else p4 = exp((1 - XIE)*SGPP[IE] + XIE*SGPP[IE + 1]);
		st = p1 + p2 + p3 + p4; //sum of N*sigma for 4 mechanism, inverse of lamda
		return -log(rng()) / st; //sample the jump distance
	}
	else if (p.type == electron)
	{
		p2 = exp((1 - XIE)*SEHEL[IE] + XIE*SEHEL[IE + 1]);
		p3 = exp((1 - XIE)*SEHIN[IE] + XIE*SEHIN[IE + 1]);
		p4 = exp((1 - XIE)*SEHBR[IE] + XIE*SEHBR[IE + 1]);
		p5 = exp((1 - XIE)*SEISI[IE] + XIE*SEISI[IE + 1]);
		if (Mode == HARD) return DS1;//hard event
		//the rest is soft event
		if (W1E[IE + 1] > -78.3)
		{
			W1 = exp((1 - XIE)*W1E[IE] + XIE*W1E[IE + 1]);
			W2 = exp((1 - XIE)*W2E[IE] + XIE*W2E[IE + 1]);
		}
		else
		{
			W1 = 0;
			W2 = 0;
		}
		if (T1E[IE + 1] > -78.3)
		{
			T1 = exp((1 - XIE)*T1E[IE] + XIE*T1E[IE + 1]);
			T2 = exp((1 - XIE)*T2E[IE] + XIE*T2E[IE + 1]);
		}
		else
		{
			T1 = 0;
			T2 = 0;
		}
		st = p2 + p3 + p4 + p5;
		double DSMaxp = DSMax;
		if (W1 > 1e-20)
		{
			KSoftI = true;
			//get the real maximum of DS
			double DSmc = 4 / st;
			if (DSmc < DSMaxp) DSMaxp = DSmc;
			else if (DSMaxp < 1e-8) DSMaxp = DSmc;
			//The value of DSMAXP is randomized to eliminate dose artifacts at the end of the first step
			DSMaxp = (0.5 + rng()*0.5)*DSMaxp;
			//Upper bound for the interaction probability along the step (including artificial energy straggling).
			double EDE0 = W1*DSMaxp;
			double VDE0 = W2*DSMaxp;
			double KW1E = (W1E[IE + 1] - W1E[IE]) / (Etab[IE + 1] - Etab[IE]);
			double KW2E = (W2E[IE + 1] - W2E[IE]) / (Etab[IE + 1] - Etab[IE]);
			double EDEM = EDE0*max(1 - 0.5*KW1E*EDE0, 0.75);
			double VDEM = VDE0*max(1 - (KW1E + 0.5*KW2E)*EDE0, 0.75);
			double W21 = VDEM / EDEM;
			double ELower;
			if (EDEM > 9.0*W21) ELower = max(p.E - (EDEM + 3.0*sqrt(VDEM)), Eabs[electron]);
			else if (EDEM > 3.0*W21) ELower = max(p.E - (EDEM + sqrt(3 * VDEM)), Eabs[electron]);
			else ELower = max(p.E - 1.5*(EDEM + W21), Eabs[electron]);
			double XE1 = (log(ELower) - ELowL) / DLE;
			int IE1 = (int)XE1;
			double XIE1 = XE1 - IE1;
			st = max(st, exp(SETOT[IE1] + (SETOT[IE1 + 1] - SETOT[IE1])*XIE1));
		}
		else KSoftI = false;

		if (T1 > 1e-20) KSoftE = true;
		else KSoftE = false;

		DST = -log(rng()) / st;
		if (DST < DSMaxp) KDelta = false;
		else
		{
			DST = DSMaxp;
			KDelta = true;
		}
		if (KSoftE == false && KSoftI == false) //hard interaction
		{
			Mode = HARD;
			DS1 = 0;
			return DST;
		}
		else //soft interaction
		{
			double DS = DST*rng();
			DS1 = DST - DS;
			return DS;
		}
	}
	else //p.type ==positron
	{
		p2 = exp((1 - XIE)*SPHEL[IE] + XIE*SPHEL[IE + 1]);
		p3 = exp((1 - XIE)*SPHIN[IE] + XIE*SPHIN[IE + 1]);
		p4 = exp((1 - XIE)*SPHBR[IE] + XIE*SPHBR[IE + 1]);
		p5 = exp((1 - XIE)*SPISI[IE] + XIE*SPISI[IE + 1]);
		p6 = exp((1 - XIE)*SPAN[IE] + XIE*SPAN[IE + 1]);
		if (Mode == HARD) return DS1;//hard event
		//the rest is soft event
		if (W1P[IE + 1] > -78.3)
		{
			W1 = exp((1 - XIE)*W1P[IE] + XIE*W1P[IE + 1]);
			W2 = exp((1 - XIE)*W2P[IE] + XIE*W2P[IE + 1]);
		}
		else
		{
			W1 = 0;
			W2 = 0;
		}
		if (T1P[IE + 1] > -78.3)
		{
			T1 = exp((1 - XIE)*T1P[IE] + XIE*T1P[IE + 1]);
			T2 = exp((1 - XIE)*T2P[IE] + XIE*T2P[IE + 1]);
		}
		else
		{
			T1 = 0;
			T2 = 0;
		}
		st = p2 + p3 + p4 + p5 + p6;
		double DSMaxp = DSMax;
		if (W1 > 1e-20)
		{
			KSoftI = true;
			//get the real maximum of DS
			double DSmc = 4 / st;
			if (DSmc < DSMaxp) DSMaxp = DSmc;
			else if (DSMaxp < 1e-8) DSMaxp = DSmc;
			//The value of DSMAXP is randomized to eliminate dose artifacts at the end of the first step
			DSMaxp = (0.5 + rng()*0.5)*DSMaxp;
			//Upper bound for the interaction probability along the step (including artificial energy straggling).
			double EDE0 = W1*DSMaxp;
			double VDE0 = W2*DSMaxp;
			double KW1P = (W1P[IE + 1] - W1P[IE]) / (Etab[IE + 1] - Etab[IE]);
			double KW2P = (W2P[IE + 1] - W2P[IE]) / (Etab[IE + 1] - Etab[IE]);
			double EDEM = EDE0*max(1 - 0.5*KW1P*EDE0, 0.75);
			double VDEM = VDE0*max(1 - (KW1P + 0.5*KW2P)*EDE0, 0.75);
			double W21 = VDEM / EDEM;
			double ELower;
			if (EDEM>9.0*W21) ELower = max(p.E - (EDEM + 3.0*sqrt(VDEM)), Eabs[positron]);
			else if (EDEM > 3.0*W21) ELower = max(p.E - (EDEM + sqrt(3 * VDEM)), Eabs[positron]);
			else ELower = max(p.E - 1.5*(EDEM + W21), Eabs[positron]);
			double XE1 = (log(ELower) - ELowL) / DLE;
			int IE1 = (int)XE1;
			double XIE1 = XE1 - IE1;
			st = max(st, exp(SPTOT[IE1] + (SPTOT[IE1 + 1] - SPTOT[IE1])*XIE1));
		}
		else KSoftI = false;

		if (T1 > 1e-20) KSoftE = true;
		else KSoftE = false;

		DST = -log(rng()) / st;
		if (DST < DSMaxp) KDelta = false;
		else
		{
			DST = DSMaxp;
			KDelta = true;
		}
		if (KSoftE == false && KSoftI == false) //hard interaction
		{
			Mode = HARD;
			DS1 = 0;
			return DST;
		}
		else //soft interaction
		{
			double DS = DST*rng();
			DS1 = DST - DS;
			return DS;
		}
	}
}

double Material::knock(int it)
{
	//int it = omp_get_thread_num(); //thread index
	PENRNG &rng = tp[it].rng;
	Particle &p = tp[it].p;
	int &IE = tp[it].IE; //left index of interval
	double &XIE = tp[it].XIE; //portion in that interval
	double &p1 = tp[it].p1;
	double &p2 = tp[it].p2;
	double &p3 = tp[it].p3;
	double &p4 = tp[it].p4;
	double &p5 = tp[it].p5;
	double &p6 = tp[it].p6;
	double &st = tp[it].st;
	double &DST = tp[it].DST;
	KnockMode &Mode = tp[it].Mode;
	bool &KSoftI = tp[it].KSoftI;
	bool &KSoftE = tp[it].KSoftE;
	bool &KDelta = tp[it].KDelta;
	double &W1 = tp[it].W1;
	double &W2 = tp[it].W2;
	double &T1 = tp[it].T1;
	double &T2 = tp[it].T2;

	double DE = 0; //energy deposited by the particle in the material

	if (p.type == photon)
	{
		SFloat den = density[p.iabsv];
		if (den == 0) den = 1.29e-3f;//set the density to be same as air
		double sm = st *(maxDensity / den);
		double sst = sm*rng();
		if (sst > st) return 0; //delta interaction
		else if (sst < p1) return GRA(); //GRA scattering
		else if (sst < p1 + p2)  return GCO(it); //Compton scattering
		else if (sst < p1 + p2 + p3) return GPH(); //Photoelectric absorption 
		else return GPP(); //Electron-positron pair production
	}
	else if (p.type == electron)
	{
		if (Mode == SOFT) //soft event will happen
		{
			Mode = HARD; //reset to hard event test
			if (KSoftI == false) DE = 0;
			else
			{
				double EDE0 = W1*DST;
				double VDE0 = W2*DST;
				double K1 = (W1E[IE + 1] - W1E[IE]) / (Etab[IE + 1] - Etab[IE]);
				double K2 = (W2E[IE + 1] - W2E[IE]) / (Etab[IE + 1] - Etab[IE]);
				double EDE = EDE0*max(1 - 0.5*K1*EDE0, 0.75); //average of DE
				double VDE = VDE0*max(1 - (K1 + 0.5*K2)*EDE0, 0.75); //variance of DE
				//Generation of random values DE distributed from 0 to infinity
				//with mean EDE and variance VDE
				double sigma = sqrt(VDE);
				if (sigma < 0.333333333*EDE) //sampling from Gaussian distribution (Box-Muller method) to accelerate the calculation
				{
					double r = 0;
					do
					{
						r = sqrt(-2.0*log(rng()))*sin(2 * PI*rng())* TRUNC;
					} while (fabs(r) > 3);
					DE = EDE + r*sigma;
				}
				else
				{
					double r = rng();
					double EDE2 = EDE*EDE;
					double VDE3 = 3 * VDE;
					if (EDE2 < VDE3)
					{
						double pnull = (VDE3 - EDE2) / (VDE3 + 3.0*EDE2);
						if (r < pnull) DE = 0;
						else DE = 1.5*(EDE + VDE / EDE)*(r - pnull) / (1.0 - pnull);
					}
					else DE = EDE + (2.0*r - 1.0)*sqrt(VDE3);
				}
			}

			p.E -= DE;
			if (p.E < Eabs[electron])
			{
				DE = p.E + DE; //lost all energy
				p.E = 0;
				return DE;
			}
			if (KSoftE == false) return DE;

			// Bielajew's randomly alternate hinge
			if (rng() > 0.5&&DE > 1.0e-3)
			{
				//calculate the energy index
				double XE = (log(p.E) - ELowL) / DLE;
				IE = (int)XE;
				XIE = XE - IE;
				if (T1E[IE + 1] > -78.3)
				{
					T1 = exp(T1E[IE] + XIE*(T1E[IE + 1] - T1E[IE]));
					T2 = exp(T2E[IE] + XIE*(T2E[IE + 1] - T2E[IE]));
				}
				else T1 = T2 = 0;
				if (T1 < 1e-35) return DE;
			}
			// 1st and 2nd moments of the angular distribution
			double emu1 = 0.5*(1.0 - exp(-DST*T1)); //average mu 1
			double emu2 = emu1 - (1.0 - exp(-DST*T2)) / 6.0; //average mu 2
			//double numerator = 2 * emu1 - 3 * emu2;
			double denominator = 1 - 2 * emu1;
			double b = (2 * emu1 - 3 * emu2) / denominator;
			double a = denominator + b;
			double r = rng();
			double cosTheta;
			if (r < a) cosTheta = 1 - 2 * b*(r / a);
			else cosTheta = 1 - 2 * (b + (1 - b)*(r - a) / (1 - a));
			p.changeByCos(cosTheta, 2*PI*rng());;
			return DE;
		}
		else // a hard event will happen
		{
			Mode = SOFT; //why set back to SOFT ???
			if (KDelta) return 0.0;//the maximum allowed step length is exceeded
			//p2 hard elastic, p3 hard inelastic, p4 hard bremsstrahlung, p5 inner-shell impact
			//p2,p3,p4, p5 have been calculated in function jump()
			double st_now = p2 + p3 + p4 + p5;
			double sst = max(st, st_now)*rng();
			if (sst < p2) //elastic scattering
			{
				double pc = (1 - XIE)*RNDCE[IE] + XIE*RNDCE[IE + 1];
				double cosTheta = 0;
				if (p.E > EELMAX) cosTheta = EEL(exp((1 - XIE)*AE[IE] + XIE*AE[IE + 1]), (1 - XIE)*BE[IE] + XIE*BE[IE + 1], pc);
				else cosTheta = EELd(pc, it);
				p.changeByCos(cosTheta, 2*PI*rng());;
				DE = 0; // no energy lost
				return DE;
			}
			else if (sst < p2 + p3)  return EIN();
			else if (sst < p2 + p3 + p4)  return EBR();
			else if (sst < st_now) return ESI();
			else return 0.0;// it may happen if st > st_now
		}

	}
	else //p.type ==positron
	{
		if (Mode == SOFT) //soft event will happen
		{
			Mode = HARD; //reset to hard event test
			if (KSoftI == false) DE = 0;
			else
			{
				double EDE0 = W1*DST;
				double VDE0 = W2*DST;
				double K1 = (W1P[IE + 1] - W1P[IE]) / (Etab[IE + 1] - Etab[IE]);
				double K2 = (W2P[IE + 1] - W2P[IE]) / (Etab[IE + 1] - Etab[IE]);
				double EDE = EDE0*max(1 - 0.5*K1*EDE0, 0.75); //average of DE
				double VDE = VDE0*max(1 - (K1 + 0.5*K2)*EDE0, 0.75); //variance of DE
				//Generation of random values DE distributed from 0 to infinity
				//with mean EDE and variance VDE
				double sigma = sqrt(VDE);
				if (sigma < 0.333333333*EDE) //sampling from Gaussian distribution (Box-Muller method) to accelerate the calculation
				{
					double r = 0;
					do
					{
						r = sqrt(-2.0*log(rng()))*sin(2 * PI*rng())* TRUNC;
					} while (fabs(r) > 3);
					DE = EDE + r*sigma;
				}
				else
				{
					double r = rng();
					double EDE2 = EDE*EDE;
					double VDE3 = 3 * VDE;
					if (EDE2 < VDE3)
					{
						double pnull = (VDE3 - EDE2) / (VDE3 + 3.0*EDE2);
						if (r < pnull) DE = 0;
						else DE = 1.5*(EDE + VDE / EDE)*(r - pnull) / (1.0 - pnull);
					}
					else DE = EDE + (2.0*r - 1.0)*sqrt(VDE3);
				}
			}

			p.E -= DE;
			if (p.E < Eabs[positron])
			{
				DE = p.E + DE + TEs; //lost all energy
				PANR();
				p.E = 0;
				return DE;
			}
			if (KSoftE == false) return DE;

			// Bielajew's randomly alternate hinge
			if (rng() > 0.5 && DE > 1.0e-3)
			{
				double XE = (log(p.E) - ELowL) / DLE;
				IE = (int)XE;
				XIE = XE - IE;
				if (T1E[IE + 1] > -78.3)
				{
					T1 = exp((1 - XIE)*T1P[IE] + XIE*T1P[IE + 1]);
					T2 = exp((1 - XIE)*T2P[IE] + XIE*T2P[IE + 1]);
				}
				else T1 = T2 = 0;
				if (T1 < 1e-35) return DE;
			}
			// 1st and 2nd moments of the angular distribution
			double emu1 = 0.5*(1.0 - exp(-DST*T1)); //average mu 1
			double emu2 = emu1 - (1.0 - exp(-DST*T2)) / 6.0; //average mu 2
			//double numerator = 2 * emu1 - 3 * emu2;
			double denominator = 1 - 2 * emu1;
			double b = (2 * emu1 - 3 * emu2) / denominator;
			double a = denominator + b;
			double r = rng();
			double cosTheta;
			if (r < a) cosTheta = 1 - 2 * b*(r / a);
			else cosTheta = 1 - 2 * (b + (1 - b)*(r - a) / (1 - a));
			p.changeByCos(cosTheta, 2*PI*rng());;
			return DE;
		}
		else // a hard event will happen
		{
			Mode = SOFT; //why set back to SOFT ???
			if (KDelta) return 0.0;//the maximum allowed step length is exceeded
			//p2 hard elastic, p3 hard inelastic, p4 hard bremsstrahlung, p5 inner-shell impact, p6 annihilation
			//p2,p3,p4, p5 have been calculated in function jump()
			double st_now = p2 + p3 + p4 + p5 + p6;
			double sst = max(st, st_now)*rng();

			if (sst < p2) //elastic scattering
			{
				double pc = (1 - XIE)*RNDCP[IE] + XIE*RNDCP[IE + 1];
				double cosTheta = 0;
				if (p.E > PELMAX) cosTheta = EEL(exp((1 - XIE)*AP[IE] + XIE*AP[IE + 1]), (1 - XIE)*BP[IE] + XIE*BP[IE + 1], pc);
				else cosTheta = PELd(pc);
				p.changeByCos(cosTheta, 2*PI*rng());;
				DE = 0; // no energy lost
				return DE;
			}
			else if (sst < p2 + p3)  return PIN();
			else if (sst < p2 + p3 + p4)  return EBR();
			else if (sst < p2 + p3 + p4 + p5) return PSI();
			else if (sst < st_now) return PAN();
			else return 0.0;// it may happen if st > st_now
		}
	}
}

double Material::GRA() //Rayleigh scattering, depends on X2COH[241], PDCOH[241], FLCOH, see eq(2.18) in PENELOPE 2005
{
	int it = omp_get_thread_num(); //thread index
	Particle &p = tp[it].p;
	PENRNG &rng = tp[it].rng;

	//no energy loss, only Direction change
	double cosTheta = 0;
	double x2max = 2.0*log(41.2148*p.E / Es); //define log(xmax^2), xmax=41.2148*E/Es
	int jm = 0, j = 0, jt = 0;
	if (x2max < X2COH[1]) jm = 0;
	else if (x2max > X2COH[239]) jm = 239;
	else jm = (int)((x2max - X2COH[0]) / FLCOH); //j must belong to [0,239]
	//linear interpolation of PImax, actually they're all log values
	double PImax = PDCOH[jm] + (PDCOH[jm + 1] - PDCOH[jm]) / (X2COH[jm + 1] - X2COH[jm]) * (x2max - X2COH[jm]);
	do
	{
		double ru = PImax + log(rng());
		//binary search for the corresponding x^2 of ru
		j = 0;
		int jr = jm + 1;
		do
		{
			jt = (j + jr) / 2;
			if (ru > PDCOH[jt]) j = jt;
			else jr = jt;
		} while (jr - j > 1);

		double dPI = PDCOH[j + 1] - PDCOH[j];
		double x2log = 0;
		if (dPI > 1e-12) x2log = X2COH[j] + (X2COH[j + 1] - X2COH[j]) / dPI *(ru - PDCOH[j]) - x2max;
		else x2log = X2COH[j] - x2max;
		cosTheta = 1.0 - 2.0*exp(x2log);
	} while (rng() > 0.5*(1.0 + cosTheta*cosTheta)); //rejection method

	p.changeByCos(cosTheta, 2*PI*rng());;//only need to change the Direction of the scattered photon

	return 0; // means there's no energy lost
}

double Material::GCO(int it)
{
	//int it = omp_get_thread_num(); //thread index
	Particle &p = tp[it].p;
	PENRNG &rng = tp[it].rng;
	SecondaryStack &secStack = tp[it].secStack;

	const double D2 = 1.4142135623731, D1 = 1.0 / D2, D12 = 0.5;

	double Ek = p.E / Es;
	double Ek2 = Ek + Ek + 1.0; //2k+1
	double Ek3 = Ek*Ek; //k^2
	double Ek1 = Ek3 - Ek2 - 1.0; //k^2-2k-2
	double tmin = 1 / Ek2; // 1/(2k+1)
	double tmin2 = tmin*tmin; // 1/(2k+1)^2
	double a1 = log(Ek2);
	double a2 = a1 + 2.0*Ek*(1.0 + Ek)*tmin2; //a1+a2 in the manual

	double Ep = 0; //energy of the scattered photon
	double cosTheta = 0; //Direction of the scattered photon
	int iSh = -1; //index of the shell ionized, starting from 0

	double tau = 0;
	double Tcos = 0;
	double S = 0; // S(E,theta)

	if (p.E > 5e6) //No Doppler broadening for E greater than 5 MeV
	{
		do
		{
			//sample tau
			do
			{
				if (rng()*a2 < a1) //i=1
					tau = pow(tmin, rng());
				else
					tau = sqrt(1.0 + rng()*(tmin2 - 1.0));
				Tcos = (1.0 + tau*(Ek1 + tau*(Ek2 + tau*Ek3))) / (Ek3*tau*(1.0 + tau*tau));
			} while (rng()>Tcos);
			Ep = tau*p.E;
			cosTheta = 1.0 - (1.0 - tau) / (Ek*tau);

			//sample which shell excited
			Tcos = Zt*rng();
			S = 0;
			iSh = -1; //for search judgment
			for (int i = 0; i < NOSCCO; ++i) //for each shell, starting from 0;
			{
				S += FCO[i];
				if (S > Tcos) { iSh = i; break; }
			}
			if (iSh == -1) iSh = NOSCCO - 1; //avoid some critical problem;
		} while (Ep > p.E - UICO[iSh]);
	}
	else //should take Doppler broadening into account
	{
		//********sample Direction of the scattered photon
		//calculate S(E,theta=pi), i.e. S0
		double S0 = 0;
		double ni[64];
		double pac[64]; //accumulation funcion of fi*theta(E-Ui)*ni, used to sample which shell being excited
		for (int i = 0; i <NOSCCO; ++i) //for each shell
		{
			if (p.E > UICO[i]) //like theta(E-Ui)
			{
				double aux = p.E*(p.E - UICO[i])*2.0; //2E*(E-Ui), auxiliary variable
				double PzJi0 = FJ0[i] * (aux - Es*UICO[i]) / (Es*sqrt(aux + aux + UICO[i] * UICO[i]));
				if (PzJi0 > 0) ni[i] = 1 - 0.5*exp(D12 - (D1 + D2*PzJi0)*(D1 + D2*PzJi0));
				else ni[i] = 0.5 * exp(D12 - (D1 - D2*PzJi0)*(D1 - D2*PzJi0));
				S0 += FCO[i] * ni[i];
			}
		}


		double cdt1 = 0;
		do
		{
			// Sampling tau
			if (rng()*a2 < a1) tau = pow(tmin, rng());
			else tau = sqrt(1.0 + rng()*(tmin2 - 1.0));
			cdt1 = (1.0 - tau) / (Ek*tau); //1-cos(theta)
			//calculate S(E, theta)
			S = 0;
			for (int i = 0; i <NOSCCO; ++i) //for each shell
			{
				if (p.E > UICO[i]) //like theta(E-Ui)
				{
					double aux = p.E*(p.E - UICO[i])*cdt1; //2E*(E-Ui), auxiliary variable
					double PzJi0 = FJ0[i] * (aux - Es*UICO[i]) / (Es*sqrt(aux + aux + UICO[i] * UICO[i]));
					if (PzJi0 > 0) ni[i] = 1 - 0.5*exp(D12 - (D1 + D2*PzJi0)*(D1 + D2*PzJi0));
					else ni[i] = 0.5 * exp(D12 - (D1 - D2*PzJi0)*(D1 - D2*PzJi0));
					S += FCO[i] * ni[i];
					pac[i] = S;
				}
				else pac[i] = S - 1e-6; //avoid sampling the E<Ui shell
			}
			Tcos = S*(1.0 + tau*(Ek1 + tau*(Ek2 + tau*Ek3))) / (Ek3*tau*(1.0 + tau*tau));
		} while (rng()*S0 > Tcos);
		cosTheta = 1.0 - cdt1; //got the Direction of the scattered photon

		//***************sampling pz;
		double pz = 0;
		double Fmax = 0;
		double Fpz = 0;
		double RS = 0;

		do
		{
			do
			{
				//sample which shell is excited
				iSh = -1; //the index of excited shell
				RS = S*rng(); //random sampled part of S
				for (int i = 0; i < NOSCCO; ++i)
				{
					if (RS < pac[i]) { iSh = i; break; }
				}
				if (iSh == -1) iSh = NOSCCO - 1; //avoid some critical problem;

				//sample pz using formula (2.56)
				double A = rng()*ni[iSh];
				if (A < 0.5) pz = (D1 - sqrt(D12 - log(A + A))) / (D2*FJ0[iSh]);
				else pz = (sqrt(D12 - log(2.0 - A - A)) - D1) / (D2*FJ0[iSh]);
			} while (pz < -1.0);

			//calculate Fmax
			double xqc = 1.0 + tau*(tau - 2.0*cosTheta);
			double AF = sqrt(xqc)*(1.0 + tau*(tau - cosTheta) / xqc);
			if (AF > 0.0) Fmax = 1.0 + AF*0.2;
			else Fmax = 1.0 - AF*0.2;
			Fpz = 1.0 + AF*max(min(pz, 0.2), -0.2);
		} while (Fmax*rng() > Fpz); //rejection selecting pz

		//******************calculate Energy of the scattered photon Ep
		double tz = pz*pz;
		double b1 = 1.0 - tz*tau*tau;
		double b2 = 1.0 - tz*tau*cosTheta;
		int signpz = (pz > 0) ? 1 : -1;
		Ep = p.E*(tau / b1)*(b2 + signpz*sqrt(fabs(b2*b2 - b1*(1.0 - tz))));
	}
	//by here, we have got Ep, cosTheta, and iSh(start from 0)

	double DE = p.E - Ep;//Energy loss of incident photon

	//*************calculate the energy of secondary electron
	double EE; //energy of secondary electron
	if (KSCO[iSh] < 10)
	{
		if (UICO[iSh] > Ecutr) EE = DE - UICO[iSh];
		else EE = DE;
	}
	else EE = DE;

	//***********by here we got Ep, cosTheta, DE, EE, cosThetaE, iSh
	int IZA = KZCO[iSh];
	int ISA = KSCO[iSh];
	double phi = 2 * PI*rng();
	if (IZA > 0 && ISA < 10) relax(IZA, ISA);


	if (EE>Eabs[electron]) // need to simulate the secondary electron
	{
		//*************calculate the Direction of the secondary electron
		double cosThetaE = 0;
		double Q2 = p.E*p.E + Ep*(Ep - 2 * p.E*cosTheta);
		if (Q2 > 1e-12) cosThetaE = (p.E - Ep*cosTheta) / sqrt(Q2);
		else cosThetaE = 1.0;

		Particle pnew = p;
		pnew.type = electron;
		pnew.changeByCos(cosThetaE, phi + PI);
		pnew.E = EE;
		//store the secondary electron
		secStack.push(pnew);
	}

	if (Ep < Eabs[photon]) // can stop the simulation
	{
		DE = p.E;
		p.E = 0.0;
	}
	else //continue simulation
	{
		p.changeByCos(cosTheta, phi);
		p.E = Ep;
	}

	return DE;
}

double Material::GPH()
{
	int it = omp_get_thread_num(); //thread index
	Particle &p = tp[it].p;
	int &IE = tp[it].IE;
	SecondaryStack &secStack = tp[it].secStack;
	PENRNG &rng = tp[it].rng;
	//to calculate
	int IZZ = 0;
	int ISH = 0;
	double EE = 0;//energy of photoelectric electron
	double EL = log(p.E);


	double pac[35]; //accumulation probability for selecting which element is simulated
	int iEel[35]; //index of corresponding energy for each element

	double ptot = 0.0;

	int iL, iR, iM;
	double dEg = 0;
	double pL = 0;
	for (int i = 0; i < elemNum; ++i)
	{
		IZZ = IZ[i];
		//binary search to locate the energy
		iL = IPHF[IZZ - 1] - 1;
		iR = IPHL[IZZ - 1] - 1;
		do
		{
			iM = (iL + iR) / 2;
			if (EL >EPH[iM]) iL = iM;
			else iR = iM;
		} while (iR - iL > 1);
		iEel[i] = iL; //store the location for future use
		dEg = EPH[iL + 1] - EPH[iL];

		if (dEg > 1e-15)	 pL = XPH[iL][0] + (XPH[iL + 1][0] - XPH[iL][0]) / dEg *(EL - EPH[iL]);
		else pL = XPH[iL][0];
		ptot += weight[i] * exp(pL);
		pac[i] = ptot;
	}

	//****  Sample the active element
	int iex = -1; //index of excited element
	double randp = rng()*SGPH[IE] / numDensity;
	for (int i = 0; i < elemNum; ++i)
	{
		if (randp < pac[i])
		{
			iex = i;
			IZZ = IZ[iex];
			break;
		}
	}
	if (iex == -1) //delta interaction
	{
		IZZ = ISH = 0;
		EE = 0;

		return 0.0; // nothing happen
	}
	//****  Sample the active shell of that element
	int iexsh = -1;
	if (NPHS[IZZ - 1] > 1) //need to select which shell is activated
	{
		int iE = iEel[iex];
		dEg = EPH[iE + 1] - EPH[iE];
		//element::gXSect &g=mo[iex].gXS[iE];
		double pt = 0;
		//prepare accumulation table
		if (dEg > 1e-15) //linear interpolation
		{
			ptot = exp(XPH[iE][0] + (XPH[iE + 1][0] - XPH[iE][0]) / dEg*(EL - EPH[iE]));
			for (int j = 1; j <= NPHS[IZZ - 1]; ++j)
			{
				pL = XPH[iE][j] + (XPH[iE + 1][j] - XPH[iE][j]) / dEg*(EL - EPH[iE]);
				pt += exp(pL);
				pac[j - 1] = pt;
			}
		}
		else //don't need linear interpolation
		{
			ptot = exp(XPH[iE][0]);
			for (int j = 0; j <NPHS[IZZ - 1]; ++j)
			{
				pt += exp(XPH[iE][j + 1]);
				pac[j] = pt;
			}
		}

		//do the selection
		randp = rng()*ptot;
		for (int i = 0; i < NPHS[IZZ - 1]; ++i)
		{
			if (randp < pac[i])
			{
				iexsh = i;
				ISH = i + 1;
				break;
			}
		}
		if (iexsh == -1) ISH = 10;
	}
	else ISH = 1;

	//by here, we've got IZZ, ISH, calculate EE next
	if (ISH < 10)
	{
		double EBB = EB[IZZ - 1][ISH - 1];
		if (EBB > Ecutr) EE = p.E - EBB;
		else
		{
			EE = p.E;
			ISH = 10;
		}
	}
	else EE = p.E;

	//by here, we've got IZZ, ISH, EE, calculate scattering angle of the electron
	if (EE > Eabs[electron])
	{
		double cosTheta = 0; //scattered Direction
		if (EE > 1e9) cosTheta = 1.0;
		else
		{
			double gamma = 1 + EE / Es;
			double gamma2 = gamma*gamma;
			double beta = sqrt(1 - 1.0 / gamma2);
			double A = 1 / beta - 1.0;
			double Aux = 0.5*beta*gamma*(gamma - 1)*(gamma - 2); //auxiliary variable
			double g0 = 2.0*(1.0 / A + Aux);
			double nu = 0, xi = 0, gnu = 0;
			do
			{
				xi = rng();
				nu = 2.0 * A / ((A + 2.0)*(A + 2.0) - 4.0 * xi)*(2.0 * xi + (A + 2.0)*sqrt(xi));
				gnu = (2.0 - nu)*(1 / (A + nu) + Aux);
			} while (rng()*g0 > gnu);
			cosTheta = 1 - nu;
		}
		Particle pe = p; //information of the electron
		pe.type = electron;
		pe.E = EE;
		pe.changeByCos(cosTheta, 2*PI*rng());;
		secStack.push(pe);//store the photo-electron in stack
	}

	if (ISH < 9) //the inner shell is excited, and need relax simulation
	{
		relax(IZZ, ISH);
	}

	double DE = p.E;
	p.E = 0; //means to end the simulation of this photon
	return DE;
}

double Material::GPP()
{
	int it = omp_get_thread_num();//thread index
	Particle &p = tp[it].p;
	SecondaryStack &secStack = tp[it].secStack;
	PENRNG &rng = tp[it].rng;

	double rki = Es / p.E;
	double eps = 0; //energy loss fraction from photon to electron
	if (p.E < 1.1e6) //very close to the lower limit
	{
		eps = rki + (1.0 - 2.0*rki)*rng();
	}
	else
	{
		double a = ZEQPP*alpha;
		double a2 = a*a;
		double t = sqrt(2 * rki);
		double t2 = t*t;
		double F0 = (-1.774 - 1.210e1*a + 1.118e1*a2)*t
			+ (8.523 + 7.326e1*a - 4.441e1*a2)*t2
			- (1.352e1 + 1.211e2*a - 9.641e1*a2)*t*t2
			+ (8.946 + 6.205e1*a - 6.341e1*a2)*t2*t2;
		double g0 = Fpp + F0;
		double bm = 4 * rki / Bpp;
		double g1, g2;
		calcG12(bm, g1, g2);
		double phi1m = g1 + g0;
		double phi2m = g2 + g0;
		double rx = 0.5 - rki;
		double P1 = 2.0 / 3.0*rx*rx*phi1m;
		P1 = P1 / (P1 + phi2m);

		double b = 0;
		//sample eps
		while (true)
		{
			if (rng() <= P1) //case 1
			{
				double rnd = 2 * rng() - 1;
				if (rnd < 0)  eps = 0.5 - rx*pow(fabs(rnd), 1.0 / 3.0);
				else eps = 0.5 + rx*pow(rnd, 1.0 / 3.0);
				b = rki / (Bpp*eps*(1 - eps));
				calcG12(b, g1, g2);
				double phi1 = max(g1 + g0, 0.0);
				if (rng()*phi1m <= phi1) break;
			}
			else //case 2
			{
				eps = rki + 2.0*rx*rng();
				b = rki / (Bpp*eps*(1 - eps));
				calcG12(b, g1, g2);
				double phi2 = max(g2 + g0, 0.0);
				if (rng()*phi2m <= phi2) break;
			}
		}
	}

	//here we've got eps
	double cosThetaE = 0, cosThetaP; //Direction of the electron/positron
	double beta = 0;
	//produced electron
	double Ee = eps*p.E - Es; //kinetic energy of electron
	if (Ee > Eabs[electron]) //need record new electron
	{
		cosThetaE = 2.0 * rng() - 1.0;
		beta = sqrt(Ee*(Ee + TEs)) / (Ee + Es);
		cosThetaE = (cosThetaE + beta) / (cosThetaE*beta + 1);
	}
	else cosThetaE = 1;

	//produced positron
	double Ep = (1 - eps)*p.E - Es; //kinetic energy of positron
	if (Ep > Eabs[positron]) //need record new positron
	{
		cosThetaP = 2.0 * rng() - 1.0;
		beta = sqrt(Ep*(Ep + TEs)) / (Ep + Es);
		cosThetaP = (cosThetaP + beta) / (cosThetaP*beta + 1);
	}
	else cosThetaP = 1;

	double DE = p.E;
	Particle pnew;
	if (Ee > Eabs[electron]) //need record new electron
	{
		pnew = p; //copy some information
		pnew.type = electron;
		pnew.changeByCos(cosThetaE, 2 * PI*rng());
		pnew.E = Ee;
		secStack.push(pnew);
	}

	if (Ep > Eabs[positron]) //need record new positron
	{
		pnew = p;
		pnew.type = positron;
		pnew.changeByCos(cosThetaP, 2 * PI*rng());
		pnew.E = Ep;
		secStack.push(pnew);
#ifdef _DEBUG
		p.u = pnew.u;
		p.v = pnew.v;
		p.w = pnew.w;
#endif

		DE = DE - TEs;// the positron carries a latent energy of 2*Es
	}
	else  PANR(); // the kinetic energy is too low, simulate the static annihilation 

	p.E = 0; //stop tracking the photon
	return DE;
}

double Material::EEL(double A, double B, double rndc)
{
	int it = omp_get_thread_num();//thread index
	PENRNG &rng = tp[it].rng;
	//used for both electron and positron
	//return the cosTheta
	double A1 = A + 1;
	double B1;
	double mu;
	if (B>0) //case I
	{
		double muav = A*A1*log(A1 / A) - A;
		B1 = 1 - B;
		double rnd0 = B1*A1*muav / (A + muav);
		double rnd = rndc + rng()*(1 - rndc);
		if (rnd < rnd0) mu = rnd*A / (B1*A1 - rnd);
		else if (rnd > rnd0 + B)
		{
			double rndmB = rnd - B;
			mu = rndmB*A / (B1*A1 - rndmB);
		}
		else mu = muav;
	}
	else //case II
	{
		double BB = -B;
		B1 = 1.0 - BB;
		double muc = rndc*A / (B1*A1 - rndc);
		double PW = B1*A*(1.0 - muc) / (A + muc);
		if (rng()*(BB + PW) < BB) mu = 0.5*(1.0 + sqrt(rng()));
		else
		{
			double rndrc = rng()*(1 - muc);
			mu = (A*rndrc + A1*muc) / (A1 - rndrc);
		}
	}

	return 1 - 2 * mu;
}

double Material::EELd(double rndc, int it)
{
	PENRNG &rng = tp[it].rng;
	//int it = omp_get_thread_num(); //thread index
	int &IE = tp[it].IE; //left index of interval
	double &XIE = tp[it].XIE; //portion in that interval

	double mu;
	int JE;
	if (rng() < XIE) JE = IE + 1;
	else JE = IE;
	double ru = rndc + rng()*(1 - rndc); //range from rndc to 1, hard event
	int itn = (int)(ru*(NP - 1));
	int i = ITLE[itn][JE] - 1;
	int j = ITUE[itn][JE] - 1;
	if (j - i >= 2)
	{
		int k;
		do
		{
			k = (i + j) / 2;
			if (ru > PSE[k][JE]) i = k;
			else j = k;
		} while (j - i > 1);
	}
	double rr = ru - PSE[i][JE];
	if (rr > 1e-16)
	{
		double xx = XSE[i][JE];
		double aa = ASE[i][JE];
		double bb = BSE[i][JE];
		double d = PSE[i + 1][JE] - PSE[i][JE];
		double cd = (1 + aa + bb)*d;
		mu = xx + (cd*rr / (d*d + (aa*d + bb*rr)*rr))*(XSE[i + 1][JE] - xx);
	}
	else mu = XSE[i][JE];

	return 1 - 2 * mu;
}

double Material::PELd(double rndc)
{
	int it = omp_get_thread_num(); //thread index
	int &IE = tp[it].IE; //left index of interval
	double &XIE = tp[it].XIE; //portion in that interval
	PENRNG &rng = tp[it].rng;

	double mu;
	int JE;
	if (rng() < XIE) JE = IE + 1;
	else JE = IE;
	double ru = rndc + rng()*(1 - rndc); //range from rndc to 1, hard event
	int itn = (int)(ru*(NP - 1));
	int i = ITLP[itn][JE] - 1;
	int j = ITUP[itn][JE] - 1;
	if (j - i >= 2)
	{
		int k;
		do
		{
			k = (i + j) / 2;
			if (ru > PSP[k][JE]) i = k;
			else j = k;
		} while (j - i > 1);
	}
	double rr = ru - PSP[i][JE];
	if (rr > 1e-16)
	{
		double xx = XSP[i][JE];
		double aa = ASP[i][JE];
		double bb = BSP[i][JE];
		double d = PSP[i + 1][JE] - PSP[i][JE];
		double cd = (1 + aa + bb)*d;
		mu = xx + (cd*rr / (d*d + (aa*d + bb*rr)*rr))*(XSP[i + 1][JE] - xx);
	}
	else mu = XSP[i][JE];

	return 1 - 2 * mu;
}

double Material::EIN()
{
	int itt = omp_get_thread_num(); //thread index
	Particle &p = tp[itt].p;
	SecondaryStack &secStack = tp[itt].secStack;
	int &IE = tp[itt].IE; //left index of interval
	double &XIE = tp[itt].XIE; //portion in that interval
	PENRNG &rng = tp[itt].rng;

	//Random sampling of hard inelastic collisions of electrons.
	//Sternheimer-Liljequist GOS model
	//DE energy deposited, EP energy of the scattered electron, cosTheta angle of the scattered electron
	//ES energy of the secondary electron, cosThetaS angle of the secondary electron
	double DE, EP, cosTheta, ES, cosThetaS;//outputs
	double E = p.E; //current energy of electron
	double delta = (1 - XIE)*DEL[IE] + XIE*DEL[IE + 1];
	if (p.E < Wcc) //delta interaction
	{
		DE = 0;
		return DE;
	}

	int JE;
	if (rng() < XIE) JE = IE + 1;
	else JE = IE;

	//selecting the active oscillator, i0 is the selected index
	double r = rng();
	int i0 = 0;
	int j0 = NOSC;
	int it = 0;
	do
	{
		it = (i0 + j0) / 2;
		if (r > EINAC[JE][it]) i0 = it;
		else j0 = it;
	} while (j0 - i0 > 1);

	//constants
	double RB = E + TEs;
	double gamma = 1 + E / Es;
	double gamma2 = gamma*gamma;
	double beta2 = 1.0 - 1.0 / gamma2;
	double amol = (1.0 - 1.0 / gamma)*(1.0 - 1.0 / gamma);

	//partial cross sections of the active oscillator
	double Wi = WRI[i0];
	double cps, cp, qm, xhdl, xhdt;
	//distant excitation
	if (Wi > Wcc && Wi < E)
	{
		cps = E*RB;
		cp = sqrt(cps);
		double xhdt0 = max(log(gamma2) - beta2 - delta, 0.0);
		if (Wi>1.0e-6 * E)
		{
			double cpp = sqrt((E - Wi)*(E - Wi + TEs));
			qm = sqrt((cp - cpp)*(cp - cpp) + Es*Es) - Es;
		}
		else
		{
			qm = Wi*Wi / (beta2 * TEs);
			qm = qm*(1.0 - qm / TEs);
		}
		if (qm < Wi)
		{
			xhdl = log(Wi*(qm + TEs) / (qm*(Wi + TEs))) / Wi;
			xhdt = xhdt0 / Wi;
		}
		else
		{
			qm = Wi;
			xhdl = 0.0;
			xhdt = 0.0;
		}
	}
	else
	{
		qm = Wi;
		cps = 0.0;
		cp = 0.0;
		xhdl = 0.0;
		xhdt = 0.0;
	}

	//Close collisions
	double WMaxC = 0.5*E;
	double WCL = max(Wcc, Wi);
	double RCL = WCL / E;
	double xhc; //hard cross section of close interaction
	if (WCL < WMaxC)
	{
		double RL1 = 1.0 - RCL;
		xhc = (amol*(0.5 - RCL) + 1.0 / RCL - 1.0 / RL1 + (1.0 - amol)*log(RCL / RL1)) / E;
	}
	else xhc = 0;

	double xhtot = xhc + xhdl + xhdt;//total cross section
	if (xhtot < 1e-35) return 0.0; //no interaction will happen

	double rxh = rng()*xhtot;
	if (rxh < xhc)// close collision is chosen
	{
		double A = 5.0*amol;
		double ARCL = A*0.5*RCL;
		double FB, RK, RK2, RKF, PHI;
		do
		{
			FB = (1.0 + ARCL)*rng();
			if (FB < 1.0) RK = RCL / (1.0 - FB*(1.0 - (RCL + RCL)));
			else RK = RCL + (FB - 1.0)*(0.5 - RCL) / ARCL;
			RK2 = RK*RK;
			RKF = RK / (1.0 - RK);
			PHI = 1.0 + RKF*RKF - RKF + amol*(RK2 + RKF);
		} while (rng()*(1.0 + A*RK2)>PHI);

		// Energy and scattering angle(primary electron).
		DE = RK*E;
		EP = E - DE;
		cosTheta = sqrt(EP*RB / (E*(RB - DE)));
		// Energy and emission angle of the delta ray.
		if (KS[i0] < 10) //?????
		{
			if (Ui[i0] > Ecutr) ES = DE - Ui[i0];
			else ES = DE;
		}
		else ES = DE;
		cosThetaS = sqrt(DE*RB / (E*(DE + TEs)));
	}
	else if (rxh < xhc + xhdl) //Hard distant longitudinal collision
	{
		DE = Wi;
		EP = E - DE;
		double QS = qm / (1.0 + qm / TEs);
		double Q = QS / (pow((QS / DE)*(1.0 + DE / TEs), rng()) - (QS / TEs));
		double QTREV = Q*(Q + TEs);
		double CPPS = EP*(EP + TEs);
		cosTheta = (CPPS + cps - QTREV) / (2.0*cp*sqrt(CPPS));
		if (cosTheta>1.0) cosTheta = 1.0;
		if (KS[i0] < 10) ES = DE - Ui[i0];
		else ES = DE;
		cosThetaS = 0.5*(DE*(E + RB - DE) + QTREV) / sqrt(cps*QTREV);
		if (cosThetaS>1.0) cosThetaS = 1.0;
	}
	else //Hard distant transverse collision
	{
		DE = Wi;
		EP = E - DE;
		cosTheta = 1.0;
		if (KS[i0] < 10) //?????
		{
			if (Ui[i0] > Ecutr) ES = DE - Ui[i0];
			else ES = DE;
		}
		else ES = DE;
		cosThetaS = 0.5;
	}


	double phi = 2 * PI*rng();
	if (ES > Eabs[electron]) //record secondary electron
	{
		Particle psec = p;
		psec.changeByCos(cosThetaS, phi + PI);
		psec.E = ES;
		psec.type = electron;
		secStack.push(psec);
	}
	if (EP > Eabs[electron])
	{
		p.E = EP;
		p.changeByCos(cosTheta, phi);
		return DE;
	}
	DE = E;
	p.E = 0.0; //lost all energy
	return DE;
}

double Material::PIN()
{
	int itt = omp_get_thread_num(); //thread index
	Particle &p = tp[itt].p;
	SecondaryStack &secStack = tp[itt].secStack;
	int &IE = tp[itt].IE; //left index of interval
	double &XIE = tp[itt].XIE; //portion in that interval
	PENRNG &rng = tp[itt].rng;

	//Random sampling of hard inelastic collisions of electrons.
	//Sternheimer-Liljequist GOS model
	//DE energy deposited, EP energy of the scattered electron, cosTheta angle of the scattered electron
	//ES energy of the secondary electron, cosThetaS angle of the secondary electron
	double DE, EP, cosTheta, ES, cosThetaS;//outputs
	double E = p.E; //current energy of electron
	double delta = (1 - XIE)*DEL[IE] + XIE*DEL[IE + 1];
	if (p.E < Wcc) //delta interaction
	{
		DE = 0;
		return DE;
	}

	int JE;
	if (rng() < XIE) JE = IE + 1;
	else JE = IE;

	//selecting the active oscillator, i0 is the selected index
	double r = rng();
	int i0 = 0;
	int j0 = NOSC;
	int it = 0;
	do
	{
		it = (i0 + j0) / 2;
		if (r > PINAC[JE][it]) i0 = it;
		else j0 = it;
	} while (j0 - i0 > 1);

	//constants
	double RB = E + TEs;
	double gamma = 1 + E / Es;
	double gamma2 = gamma*gamma;
	double beta2 = 1.0 - 1.0 / gamma2;
	double G12 = (gamma + 1.0)*(gamma + 1.0);
	double amol = (E / (E + Es))*(E / (E + Es));
	double BHA1 = amol*(2.0*G12 - 1.0) / (gamma2 - 1.0);
	double BHA2 = amol*(3.0 + 1.0 / G12);
	double BHA3 = amol*2.0*gamma*(gamma - 1.0) / G12;
	double BHA4 = amol*(gamma - 1.0)*(gamma - 1.0) / G12;

	//partial cross sections of the active oscillator
	double Wi = WRI[i0];
	double cps, cp, qm, xhdl, xhdt;
	//distant excitation
	if (Wi > Wcc && Wi < E)
	{
		cps = E*RB;
		cp = sqrt(cps);
		double xhdt0 = max(log(gamma2) - beta2 - delta, 0.0);
		if (Wi>1.0e-6 * E)
		{
			double cpp = sqrt((E - Wi)*(E - Wi + TEs));
			qm = sqrt((cp - cpp)*(cp - cpp) + Es*Es) - Es;
		}
		else
		{
			qm = Wi*Wi / (beta2 * TEs);
			qm = qm*(1.0 - qm / TEs);
		}
		if (qm < Wi)
		{
			xhdl = log(Wi*(qm + TEs) / (qm*(Wi + TEs))) / Wi;
			xhdt = xhdt0 / Wi;
		}
		else
		{
			qm = Wi;
			xhdl = 0.0;
			xhdt = 0.0;
		}
	}
	else
	{
		qm = Wi;
		cps = 0.0;
		cp = 0.0;
		xhdl = 0.0;
		xhdt = 0.0;
	}

	//Close collisions
	double WMaxC = E;
	double WCL = max(Wcc, Wi);
	double RCL = WCL / E;
	double xhc; //hard cross section of close interaction
	if (WCL < WMaxC)
	{
		double RL1 = 1.0 - RCL;
		xhc = ((1.0 / RCL - 1.0) + BHA1*log(RCL) + BHA2*RL1
			+ (BHA3 / 2.0)*(RCL*RCL - 1.0) + (BHA4 / 3.0)*(1.0 - RCL*RCL*RCL)) / E;
	}
	else xhc = 0.0;

	double xhtot = xhc + xhdl + xhdt;//total cross section
	if (xhtot < 1e-35) return 0.0; //no interaction will happen

	double rxh = rng()*xhtot;
	if (rxh < xhc)// close collision is chosen
	{
		double RK, PHI;
		do
		{
			RK = RCL / (1.0 - rng()*(1.0 - RCL));
			PHI = 1.0 - RK*(BHA1 - RK*(BHA2 - RK*(BHA3 - BHA4*RK)));
		} while (rng() > PHI);

		// Energy and scattering angle(primary electron).
		DE = RK*E;
		EP = E - DE;
		cosTheta = sqrt(EP*RB / (E*(RB - DE)));
		// Energy and emission angle of the delta ray.
		if (KS[i0] < 10) //?????
		{
			if (Ui[i0] > Ecutr) ES = DE - Ui[i0];
			else ES = DE;
		}
		else ES = DE;
		cosThetaS = sqrt(DE*RB / (E*(DE + TEs)));
	}
	else
	{
		DE = Wi;
		EP = E - DE;
		if (rxh < xhc + xhdl) //Hard distant longitudinal collision
		{
			double QS = qm / (1.0 + qm / TEs);
			double Q = QS / (pow((QS / DE)*(1.0 + DE / TEs), rng()) - (QS / TEs));
			double QTREV = Q*(Q + TEs);
			double CPPS = EP*(EP + TEs);
			cosTheta = (CPPS + cps - QTREV) / (2.0*cp*sqrt(CPPS));
			if (cosTheta > 1.0) cosTheta = 1.0;
			if (KS[i0] < 10) ES = DE - Ui[i0];
			else ES = DE;
			cosThetaS = 0.5*(DE*(E + RB - DE) + QTREV) / sqrt(cps*QTREV);
			if (cosThetaS>1.0) cosThetaS = 1.0;
		}
		else //Hard distant transverse collision
		{
			cosTheta = 1.0;
			if (KS[i0] < 10) //?????
			{
				if (Ui[i0] > Ecutr) ES = DE - Ui[i0];
				else ES = DE;
			}
			else ES = DE;
			cosThetaS = 0.5;
		}
	}
	double phi = 2 * PI*rng();
	if (ES > Eabs[electron]) //record secondary electron
	{
		Particle psec = p;
		psec.changeByCos(cosThetaS, phi + PI);
		psec.E = ES;
		psec.type = electron;
		secStack.push(psec);
	}
	if (EP > Eabs[positron])
	{
		p.E = EP;
		p.changeByCos(cosTheta, phi);
		return DE;
	}
	DE = E + TEs;
	p.E = 0.0; //lost all energy
	PANR();
	return DE;
}

double Material::EBR()
{
	int it = omp_get_thread_num();//thread index
	Particle &p = tp[it].p;
	SecondaryStack &secStack = tp[it].secStack;
	int &IE = tp[it].IE; //left index of interval
	double &XIE = tp[it].XIE; //portion in that interval
	PENRNG &rng = tp[it].rng;

	if (p.E < Wcr) return 0.0;// no interaction happen
	int i, j, k, m; //index of energy
	if (rng() < XIE) i = IE + 1;
	else i = IE;

	double W = 0;//(reduced) photon energy
	double pt, W1, W2, bj, aj, pmax;
	do
	{
		pt = PBcut[i] + rng()*(PACB[i][NBW - 1] - PBcut[i]); //uniform random number from PACBcut[i] to PACB[i][NWB-1]
		//locate where the kappa is
		j = 0;
		k = NBW - 1;
		do
		{
			m = (j + k) / 2;
			if (pt > PACB[i][m]) j = m;
			else k = m;
		} while (k - j > 1);
		//j = binarySearch(&PACB[i][0], NBW, pt);
		//
		W1 = WB[j];
		W2 = WB[j + 1];
		bj = (PDFB[i][j + 1] - PDFB[i][j]) / (W2 - W1);
		aj = PDFB[i][j] - bj*W1;
		if (W1 < WBcut[i]) W1 = WBcut[i];
		if (W2 < W1)
		{
			Log("PENELOPE EBR Warning: Conflicting end-point values in EBR()");
			W = W1;
			break;
		}
		pmax = max(aj + bj*W1, aj + bj*W2);
		do //rejection sampling
		{
			W = W1*pow(W2 / W1, rng());
		} while (rng()*pmax > aj + bj*W);
		W *= p.E; //this is the radiated photon energy
	} while (W < Wcr); //radiated energy should be larger than Wcr, may not end???


	//here we got the radiated energy W. Next, do the angular sampling
	if (W > Eabs[photon]) //sample the direction of the photon
	{
		double beta = sqrt(p.E*(p.E + TEs)) / (p.E + Es);
		double cosTheta = 1.0;
		if (p.E > 500e3) //use a simplified dipole distribution for E>500keV
		{
			cosTheta = 2.0*rng() - 1.0;
			if (rng() > 0.75)
			{
				if (cosTheta < 0.0) cosTheta = -pow(-cosTheta, 1.0 / 3.0);
				else cosTheta = pow(cosTheta, 1.0 / 3.0);
			}
			cosTheta = (cosTheta + beta) / (1.0 + beta*cosTheta);
		}
		else //use  fitting expression (3.169)
		{
			int IEE = 0;
			if (beta >BET[5]) IEE = 5;
			else if (beta < BET[0]) IEE = 0;
			else //search for the region where beta located
			{
				IEE = 0;
				int IEE1 = 5;
				int IEEM;
				do
				{
					IEEM = (IEE + IEE1) / 2;
					if (beta>BET[IEEM]) IEE = IEEM;
					else IEE1 = IEEM;
				} while (IEE1 - IEE>1);
			}

			//calculate A and betap first
			double rk = W / p.E * 20;
			int ik = min((int)rk, 19);
			double qL = BP1[IEE][ik][0] + beta*(BP1[IEE][ik][1] + beta*(BP1[IEE][ik][2] + beta*BP1[IEE][ik][3]));
			double qR = BP1[IEE][ik + 1][0] + beta*(BP1[IEE][ik + 1][1] + beta*(BP1[IEE][ik + 1][2] + beta*BP1[IEE][ik + 1][3]));
			double q1 = qL + (rk - ik)*(qR - qL);

			qL = BP2[IEE][ik][0] + beta*(BP2[IEE][ik][1] + beta*(BP2[IEE][ik][2] + beta*BP2[IEE][ik][3]));
			qR = BP2[IEE][ik + 1][0] + beta*(BP2[IEE][ik + 1][1] + beta*(BP2[IEE][ik + 1][2] + beta*BP2[IEE][ik + 1][3]));
			double q2 = qL + (rk - ik)*(qR - qL);

			double A = min(exp(q1) / beta, 1.0);
			double betap = min(max(beta*(1.0 + q2 / beta), 0.0), 0.999999999);
			if (rng() < A)
			{
				do
				{
					cosTheta = 2.0*rng() - 1.0;
				} while (2 * rng() > 1 + cosTheta*cosTheta);
			}
			else
			{
				do
				{
					cosTheta = 2.0*rng() - 1.0;
				} while (rng() > 1 - cosTheta*cosTheta);
			}
			cosTheta = (cosTheta + betap) / (1 + betap*cosTheta);
		}

		//here we got the scattering angle of photon
		Particle pnew = p;
		pnew.changeByCos(cosTheta, 2*PI*rng());;
		pnew.E = W;
		pnew.type = photon;
		secStack.push(pnew);
	}

	p.E -= W;
	if (p.E > Eabs[p.type]) return W;
	//lost all energy
	W = p.E + W;
	p.E = 0;
	if (p.type == positron)//positrion lost more energy and need call PANR()
	{
		W += TEs;
		PANR();
	}
	return W;
}

double Material::ESI()
{
	int it = omp_get_thread_num(); //thread index
	Particle &p = tp[it].p;
	int &IE = tp[it].IE; //left index of interval
	double &XIE = tp[it].XIE; //portion in that interval
	PENRNG &rng = tp[it].rng;

	//Ionization of inner shells by impact of electrons.
	//Output arguments, both start from 1 :
	//IZZ ....atomic number of the element where ionization has ocurred.
	//ISH ....atomic electron shell that has been  ionized.
	double PACSI[120];
	int IZSI[120], ISHSI[120];
	int kt = 0;
	PACSI[0] = 0;
	int IZZ, ISH, INDC;
	double SITOT = 0, XSI = 0;

	//prepare accumulation probability
	for (int j = 0; j < elemNum; ++j)
	{
		IZZ = IZ[j] - 1;
		INDC = IESIF[IZZ] - 1;
		for (int ish = 0; ish < NSESI[IZZ]; ++ish)
		{
			double Wcut = EB[IZZ][ish];
			if (Wcut>Ecutr && Wcut < p.E)
			{
				XSI = exp(XESI[INDC + IE][ish] + XIE*(XESI[INDC + IE + 1][ish] - XESI[INDC + IE][ish]));
				if (XSI > 1.1e-35)
				{
					SITOT += XSI*weight[j];
					IZSI[kt] = IZZ + 1;
					ISHSI[kt] = ish + 1;
					++kt;
					PACSI[kt] = SITOT;
				}
			}
		}
	}

	//sampling IZZ and ISH pair
	if (kt == 0) return 0.0;//do nothing

	double rx = PACSI[kt] * rng();
	int is = 0, im;
	do
	{
		im = (is + kt) / 2;
		if (rx > PACSI[im]) is = im;
		else kt = im;
	} while (kt - is>1);

	IZZ = IZSI[is];
	ISH = ISHSI[is];

	relax(IZZ, ISH);

	return 0.0; //no energy loss
}

double Material::PSI()
{
	int it = omp_get_thread_num(); //thread index
	Particle &p = tp[it].p;
	int &IE = tp[it].IE; //left index of interval
	double &XIE = tp[it].XIE; //portion in that interval
	PENRNG &rng = tp[it].rng;

	//Ionization of inner shells by impact of electrons.
	//Output arguments, both start from 1 :
	//IZZ ....atomic number of the element where ionization has ocurred.
	//ISH ....atomic electron shell that has been  ionized.
	double PACSI[120];
	int IZSI[120], ISHSI[120];
	int kt = 0;
	PACSI[0] = 0;
	int IZZ, ISH, INDC;
	double SITOT = 0, XSI = 0;

	//prepare accumulation probability
	for (int j = 0; j < elemNum; ++j)
	{
		IZZ = IZ[j] - 1;
		INDC = IPSIF[IZZ] - 1;
		for (int ish = 0; ish < NSPSI[IZZ]; ++ish)
		{
			double Wcut = EB[IZZ][ish];
			if (Wcut>Ecutr && Wcut < p.E)
			{
				XSI = exp(XPSI[INDC + IE][ish] + XIE*(XPSI[INDC + IE + 1][ish] - XPSI[INDC + IE][ish]));
				if (XSI > 1.1e-35)
				{
					SITOT += XSI*weight[j];
					IZSI[kt] = IZZ + 1;
					ISHSI[kt] = ish + 1;
					++kt;
					PACSI[kt] = SITOT;
				}
			}
		}
	}

	//sampling IZZ and ISH pair
	if (kt == 0) return 0.0;//do nothing

	double rx = PACSI[kt] * rng();
	int is = 0, im;
	do
	{
		im = (is + kt) / 2;
		if (rx > PACSI[im]) is = im;
		else kt = im;
	} while (kt - is > 1);

	IZZ = IZSI[is];
	ISH = ISHSI[is];

	relax(IZZ, ISH);

	return 0.0; //no energy loss
}

double Material::PAN()
{
	int it = omp_get_thread_num(); //thread index
	Particle &p = tp[it].p;
	SecondaryStack &secStack = tp[it].secStack;
	PENRNG &rng = tp[it].rng;

	double DE = 0;
	double E1, cosTheta1, E2, cosTheta2;
	if (p.E < Eabs[positron]) //Slow positrons (assumed at rest)
	{
		E1 = 0.5*(p.E + TEs);
		E2 = E1;
		cosTheta1 = -1.0 + 2.0*rng();
		cosTheta2 = -cosTheta1;
	}
	else
	{
		//Annihilation in flight(two photons with energy and directions
		//determined from the DCS and energy - momentum conservation).
		double gamma = 1.0 + max(p.E, 1.0) / Es;
		double GAM21 = sqrt(gamma*gamma - 1.0);
		double ANI = 1.0 + gamma;
		double CHIMIN = 1.0 / (ANI + GAM21);
		double RCHI = (1.0 - CHIMIN) / CHIMIN;
		double GT0 = ANI*ANI - 2.0;
		double CHI, GREJ;
		do
		{
			CHI = CHIMIN*pow(RCHI, rng());
			GREJ = ANI*ANI*(1.0 - CHI) + gamma + gamma - 1.0 / CHI;
		} while (rng()*GT0 > GREJ);

		double DET = p.E + TEs;
		E1 = CHI*DET;
		cosTheta1 = (ANI - 1.0 / CHI) / GAM21;
		double CHIP = 1.0 - CHI;
		E2 = DET - E1;
		cosTheta2 = (ANI - 1.0 / CHIP) / GAM21;
	}

	//need to store the photons
	double phi = 2 * PI*rng();
	if (E1 > Eabs[photon])
	{
		Particle pnew = p;
		pnew.changeByCos(cosTheta1, phi);
		pnew.E = E1;
		pnew.type = photon;
		secStack.push(pnew);
	}
	if (E2 > Eabs[photon])
	{
		Particle pnew = p;
		pnew.changeByCos(cosTheta2, phi + PI);
		pnew.E = E2;
		pnew.type = photon;
		secStack.push(pnew);
	}
	DE = p.E + TEs;
	p.E = 0.0;

	return DE;
}

//other process
void Material::PANR()
{
	int it = omp_get_thread_num(); //thread index
	Particle &p = tp[it].p;
	SecondaryStack &secStack = tp[it].secStack;
	PENRNG &rng = tp[it].rng;

	//simulate the static Annihilation of electron and positron
	Particle p1 = p; //copy position
	p1.type = photon;
	double cosTheta = -1 + 2 * rng();
	p1.changeByCos(cosTheta, 2*PI*rng());;
	p1.E = Es;

	Particle p2 = p1; //recoiled particle with opposite velocity
	p2.u = -p1.u;
	p2.v = -p1.v;
	p2.w = -p1.w;

	secStack.push(p1);
	secStack.push(p2);
}

void Material::relax(int IZ, int IS)
{
	int it = omp_get_thread_num();//thread index
	Particle &p = tp[it].p;
	SecondaryStack &secStack = tp[it].secStack;
	PENRNG &rng = tp[it].rng;

	if (IZ < 6 || IS > 9) return;
	--IZ;
	--IS;
	//If the shell ionization energy is less than Ecutr, the cascade is not followed.
	if (EB[IZ][IS] < Ecutr) return;
	const int stackDepth = 256;
	int isv[stackDepth]; //vacant shell index stack
	int nv = 0; // stack index, pointing to the first element
	isv[0] = IS; //index of initial vacant shell
	int isp = 0;
	int KF, KL, K1, KS;
	double RN, RX;
	int is1, is2;
	ParticleType ptype;// type of generated particle
	do
	{
		//fetch one vacancy from the stack
		isp = isv[nv];
		KF = IFIRST[IZ][isp] - 1;
		KL = ILAST[IZ][isp] - 1;
		--nv;
		//select one transition
		if (KL > KF)
		{
			RN = rng()*(KL - KF + 1); // KL-KF+1 is the point number
			K1 = (int)RN;
			RX = RN - K1; //the decimal part
			if (RX > FR[KF + K1]) KS = IAL[KF + K1] - 1;
			else KS = KF + K1;
		}
		else KS = KF;
		//KS is the index 
		is1 = IS1[KS];
		is2 = IS2[KS];
		if (is2 == 0) // radiative 
		{
			ptype = photon;
			if (is1 < 10 && EB[IZ][is1 - 1] > Ecutr)
			{
				++nv;
				isv[nv] = is1 - 1;
			}
		}
		else // Non-radiative 
		{
			ptype = electron;
			if (is1 < 10 && EB[IZ][is1 - 1] > Ecutr)
			{
				++nv;
				isv[nv] = is1 - 1; // note -1 is necessary since isv starts from 0
			}
			if (is2 < 10 && EB[IZ][is2 - 1] > Ecutr)
			{
				++nv;
				isv[nv] = is2 - 1; // note -1 is necessary since isv starts from 0
			}
		}

		//store the emitted particle
		if (ET[KS]>Eabs[ptype])
		{
			Particle pnew = p; //Position and material index copied
			pnew.type = ptype;
			pnew.w = -1.0 + 2.0*rng();
			double sinTheta = sqrt(1 - pnew.w*pnew.w);
			double phi = 2 * PI* rng();
			pnew.u = sinTheta*cos(phi);
			pnew.v = sinTheta*sin(phi);
			pnew.E = ET[KS];

			secStack.push(pnew); //store second particle
		}

	} while (nv >= 0);
}
/************************End: implemtents of class Material **************************/

void loadPenelopeCommonData()
{
	FILE *fp = NULL;
	if (true)
	{
		fp = fopen("Materials/common.dat", "rb");
		if (fp == NULL) exitApp("Cannot open the file common.dat!");
		fread(&ELow, sizeof(double), 1, fp);
		fread(&EUp, sizeof(double), 1, fp);
		fread(&FLCOH, sizeof(double), 1, fp);
		fread(EPH, sizeof(double), 8000, fp);
		for (int i = 0; i < 10; ++i)
		{
			for (int j = 0; j < 8000; ++j) fread(&XPH[j][i], sizeof(double), 1, fp);
		}
		fread(IPHF, sizeof(int), 99, fp);
		fread(IPHL, sizeof(int), 99, fp);
		fread(NPHS, sizeof(int), 99, fp);
		fread(ET, sizeof(double), 15000, fp);
		fread(FR, sizeof(double), 15000, fp);
		fread(IAL, sizeof(int), 15000, fp);
		fread(IS1, sizeof(int), 15000, fp);
		fread(IS2, sizeof(int), 15000, fp);

		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 99; ++j) fread(&IFIRST[j][i], sizeof(int), 1, fp);
		}

		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 99; ++j) fread(&ILAST[j][i], sizeof(int), 1, fp);
		}

		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 6000; ++j) fread(&XESI[j][i], sizeof(double), 1, fp);
		}
		fread(IESIF, sizeof(int), 99, fp);
		fread(NSESI, sizeof(int), 99, fp);
		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 6000; ++j) fread(&XPSI[j][i], sizeof(double), 1, fp);
		}
		fread(IPSIF, sizeof(int), 99, fp);
		fread(NSPSI, sizeof(int), 99, fp);
		for (int i = 0; i < 30; ++i)
		{
			for (int j = 0; j < 99; ++j) fread(&EB[j][i], sizeof(double), 1, fp);
		}
	}
	else
	{
		fp = fopen("common.dat", "r");
		if (fp == NULL) exitApp("Cannot open the file common.dat!");

		readline(fp);
		readline(fp);
		fscanf(fp, "%lf%lf%lf\n", &ELow, &EUp, &FLCOH);
		readline(fp);
		for (int i = 0; i < 8000; ++i) fscanf(fp, "%lf", &EPH[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 10; ++i)
		{
			for (int j = 0; j < 8000; ++j) fscanf(fp, "%lf", &XPH[j][i]); fgetc(fp);
		}
		readline(fp);
		for (int i = 0; i < 99; ++i) fscanf(fp, "%d", &IPHF[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 99; ++i) fscanf(fp, "%d", &IPHL[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 99; ++i) fscanf(fp, "%d", &NPHS[i]); fgetc(fp);

		readline(fp);
		for (int i = 0; i < 15000; ++i) fscanf(fp, "%lf", &ET[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 15000; ++i) fscanf(fp, "%lf", &FR[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 15000; ++i) fscanf(fp, "%d", &IAL[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 15000; ++i) fscanf(fp, "%d", &IS1[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 15000; ++i) fscanf(fp, "%d", &IS2[i]); fgetc(fp);

		readline(fp);
		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 99; ++j) fscanf(fp, "%d", &IFIRST[j][i]); fgetc(fp);
		}
		readline(fp);
		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 99; ++j) fscanf(fp, "%d", &ILAST[j][i]); fgetc(fp);
		}

		readline(fp);
		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 6000; ++j) fscanf(fp, "%lf", &XESI[j][i]); fgetc(fp);
		}
		readline(fp);
		for (int i = 0; i < 99; ++i) fscanf(fp, "%d", &IESIF[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 99; ++i) fscanf(fp, "%d", &NSESI[i]); fgetc(fp);

		readline(fp);
		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 6000; ++j) fscanf(fp, "%lf", &XPSI[j][i]); fgetc(fp);
		}
		readline(fp);
		for (int i = 0; i < 99; ++i) fscanf(fp, "%d", &IPSIF[i]); fgetc(fp);
		readline(fp);
		for (int i = 0; i < 99; ++i) fscanf(fp, "%d", &NSPSI[i]); fgetc(fp);

		readline(fp);
		for (int i = 0; i < 30; ++i)
		{
			for (int j = 0; j < 99; ++j) fscanf(fp, "%lf", &EB[j][i]); fgetc(fp);
		}
	}

	fclose(fp);

	//init BET[6], used in EBR()
	double E[6] = { 1.0e3, 5.0e3, 1.0e4, 5.0e4, 1.0e5, 5.0e5 };
	for (int i = 0; i < 6; ++i) BET[i] = sqrt(E[i] * (E[i] + TEs)) / (E[i] + Es);
	//init energy grid
	ELowL = log(ELow);
	double EUpL = log(EUp);
	DLE = (EUpL - ELowL) / (NEGrid - 1);
	for (int i = 0; i < NEGrid; ++i)
	{
		ELtab[i] = ELowL + i*DLE; //log(E) Energy table, uniform grid
		Etab[i] = exp(ELtab[i]);
	}

	//init accelerating energy grid
	LDE = (EUp - ELow) / (LNEGrid - 1);



}

int pID = 0;

void initPenelope(ConfigFile *cf, Material*& mat)
{
    //determine if we need to prepare the data tables
    string useLast;
    cf->getValue("useLastMats", useLast);
    bool useLastMats = useLast.compare("yes") == 0 ? true : false;
    bool fileExist = existFile("Materials/common.dat");
    ConfigFile *subCF = NULL;
    cf->resetSearchIndex();
    int _NMAT = 1;
    for (_NMAT = 1;; ++_NMAT)
    {
        subCF = cf->getBlock("MAT", true);
        if (NULL == subCF)
        {
            if (1 == _NMAT)
            {
                exitApp("cannot load any material parameters!"); //cannot even find the first material
            }
            else
            {
                --_NMAT;
                break;//at least loaded one material
            }
        }
        char ic[20];
        sprintf(ic, "Materials/mat%d.dat", _NMAT);
        fileExist = existFile(ic);
        int IDNumber = 0;
        subCF->getValue("ID number", IDNumber);
    }

    if (useLastMats&&fileExist)
    {
        Log("Penelope tries to load material files generated last time...\n");
        Log.warningColor();
        Log("Make sure that EMax and material parameters weren't changed! \n");
        Log.normalColor();
    }
    else //generate the data file
    {
        if (pID == 0) // Only the main process takes care of data table preparation
        {
            RunTimeCounter rc;
            Log("Penelope will prepare material data tables in advance...\n");
            //<<load the lib functions
#ifdef WIN32
            //SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, FALSE); //try to unload the control-C handler
			HINSTANCE hlib = ::LoadLibrary("PINIT.dll");
#else
            void* hlib = dlopen("./PINIT.so", RTLD_LAZY);
#endif
            if (hlib == NULL) exitApp("Cannot open PINIT.dll or PINIT.so");
            typedef void (CallStyle *SETPARA)(int*, double*, double*, double*, double*, double*, double*, double*);
            typedef void (CallStyle *PINIT)(double*, int*, int*, int*);
            SETPARA setPara = (SETPARA)getLibAddress(hlib, "setpara");
            if (setPara == NULL) exitApp("cannot load function setpara");
            PINIT initMat = (PINIT)getLibAddress(hlib, "pinit");
            if (initMat == NULL) exitApp("cannot load function pinit");
            //>>

            double EMax = SourceHead_Emax();
            cf->getValue("EMax", EMax); //may overwrite EMax


            ConfigFile *subCF = NULL;
            _NMAT = 1;
            double Eabs_e, Eabs_g, Eabs_p, C1, C2, Wcc, Wcr;
            int IDNumber;

            //read parameters of each material
            FILE* fp = NULL;
            fp = fopen("temp_filename.txt", "w");
            if (NULL == fp) exitApp("cannot create temp_filename.txt");
            fprintf(fp, "Materials/\n"); //tell the FORTRAN code in which directory it should export

            if (!fs::exists(fs::path("Materials"))) fs::create_directories(fs::path("Materials"));
            string dbpath;
            cf->getValue("material database", dbpath);
            //search for the definition of beams
            cf->resetSearchIndex();
            for (_NMAT = 1;; ++_NMAT)
            {
                subCF = cf->getBlock("MAT", true);
                if (NULL == subCF)
                {
                    if (1 == _NMAT)
                    {
                        exitApp("cannot load any material parameters!"); //cannot even find the first material
                    }
                    else
                    {
                        --_NMAT;
                        break;//at least loaded one material
                    }
                }
                subCF->getValue("Eabs_e", Eabs_e);
                subCF->getValue("Eabs_g", Eabs_g);
                subCF->getValue("Eabs_p", Eabs_p);
                subCF->getValue("C1", C1);
                subCF->getValue("C2", C2);
                subCF->getValue("Wcc", Wcc);
                subCF->getValue("Wcr", Wcr);
                subCF->getValue("ID number", IDNumber);

                if (Eabs_e == 0 || Eabs_g == 0 || Eabs_p == 0 || Wcc == 0 || Wcr == 0) //set the default cutoff energy
                {
                    Eabs_e = Eabs_p = 1e-2*EMax;
                    Eabs_g = 1e-3*EMax;
                    Wcc = Eabs_e;
                    Wcr = Eabs_g;
                }

                //search for the material file
                char ic[10];
                sprintf(ic, "%d", IDNumber);
                string filter(ic);
                filter += " ";
                string match = filter;
                matchFile("Materials", match);
                if (match.empty()) //need to exact from the database
                {
                    if (1 <= IDNumber && IDNumber <= 280)
                    {
                        string cmd = "7z e -oMaterials \"";
                        cmd += dbpath;
                        cmd += "\" \"";
                        cmd += ic;
                        cmd += " *\"";
                        system(cmd.c_str());

                        int numTry = 0;
                        do
                        {
                            ++numTry;
                            // wait some time for the exaction to be finished
                            std::this_thread::sleep_for(std::chrono::milliseconds(500));
                            match = filter;
                            matchFile("Materials", match);
                            if (numTry>5) exitApp("Failed to exact the material files from PENELOPE database. Please check the database file's path");
                        } while (match.empty());
                    }
                    else exitApp("Cannot find the designated material file in the database or in the Materials directory");
                }
                fprintf(fp, "%s\n", match.c_str()); //only main process write the file
                setPara(&_NMAT, &Eabs_e, &Eabs_g, &Eabs_p, &C1, &C2, &Wcc, &Wcr);
            }
            fclose(fp); // close the file

            int info = 0; // output information level
            int useBinaryFile = 1;
            initMat(&EMax, &_NMAT, &info, &useBinaryFile); // generate the data files
            if (hlib) freeLib(hlib);
            fs::remove(fs::path("temp_filename.txt"));

            Log("\nIt costs %f seconds to prepare Penelope materials. ", rc.stop());
        }
    }

#ifdef USE_MPI
    //let's make sure the data files have been created by synchronizing all processes using MPI_Bcast
	char ic[10];
	MPI_Bcast(ic, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

    //load the material data to CPU memory
    mat = new Material[_NMAT];
    for (int i = 0; i < _NMAT; ++i) mat[i].load(i + 1);

    loadPenelopeCommonData();
}


void executeJob(const char* configFileName, MPS& configMacro) //execute one job according to the config file
{
	RunTimeCounter totalTime;

	double DSMax = 0.1;//unit cm
	
	Material* mat = NULL;

	ConfigFile cf(configFileName); //parse the total config file
	//find out where is the config file located
	string cdir(configFileName);
	size_t pos = cdir.find_last_of("\\/");
	cdir.erase(++pos, string::npos);
	cf.macroReplace(string("$cdir$"), cdir);
	cf.macroReplace(configMacro);

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

	ConfigFile *subCF = NULL;
	int NThread = get_thread_num();

	subCF = cf.getBlock("SOURCEHEAD");
	SourceHead_Init(subCF);
	subCF->getValue("NThread", NThread); //the number in SourceHead file will take over the default one

	subCF = cf.getBlock("PENELOPE");
	subCF->getValue("DSMax", DSMax);
	initPenelope(subCF, mat);

	subCF = cf.getBlock("PHANTOM");
	Phantom tank;
	tank.loadPhantom(subCF);

	double fSIMU = 1;
	cf.getValue("NSIMU", fSIMU);

	tp = new ThreadParameters[NThread];
	/***************************************************************************************/
	Log("\nCalculating dose with %d thread(s), please wait patiently...\n\n", NThread);
	RunTimeCounter rc; //counting for the total simulation time
#ifdef STEP_DEBUG_STATUS
	FILE* fp = fopen("step_debug.dat", "w");
	if (NULL == fp) exitApp("Cannot create debug output files!");
	int nStepDebug = 0;
#endif

#ifdef ACCU_DEBUG
	int nIncident = 0;
#endif
	
#ifdef USE_OPENMP
#pragma omp parallel num_threads(NThread)
#endif
	{
		double fNSIMU = fSIMU/NThread;
		int it = omp_get_thread_num(); //thread index starting from 0
		double DS, DE;
		bool cross;
		Particle pars[400];//the max number of photon coming from the ViewRay head per time
		Particle &p = tp[it].p;
		SecondaryStack &secStack = tp[it].secStack;
		KnockMode &Mode = tp[it].Mode;
		PRNG rng; //RNG designed only for the source head
		rng.init(1234 + it); //initiate the RNG for source head
		tp[it].rng.init(1234, it); //initiate the RNG for the kernel simulation

		RunTimeCounter percentSpeed;
		double hist = 0;
		double percentHist = fNSIMU*0.01;
		int percent = 1;
		for (; hist <= fNSIMU; ++hist) //percentage need to be modified?????
		{
			int nSample = SourceHead_Sample(&rng, pars);
			for (int ipp = 0; ipp < nSample; ++ipp)
			{
				if (tank.lineIn(pars[ipp])) //the particle can enter the phantom
				{
					p = pars[ipp];
					if (p.E<1.00000001*ELow || p.E>0.99999999*EUp) exitApp("E out of range in function start()");
					Mode = SOFT;

					while (true)
					{
						//here we know the particle is inside the object
						DS = mat[0].jump(DSMax, it); //jump some distance
						cross = tank.move(p, DS); //move in the phantom, and judge if it crossed the boundary
						if (!cross) //still inside the object
						{
// #ifdef STEP_DEBUG
// 							//fprintf(fp, "v=(%.9f,%.9f,%.9f), x=(%.9f,%.9f,%.9f), E = %.9f\n", p.u, p.v, p.w, p.x + tank.DX*p.ivx, p.y + tank.DY*p.ivy, p.z + tank.DZ*p.ivz, p.E*1e-3);
// 							double temp[7];
// 							temp[0] = p.u;
// 							temp[1] = p.v;
// 							temp[2] = p.w;
// 							temp[3] = p.x + tank.DX*p.ivx;
// 							temp[4] = p.y + tank.DY*p.ivy;
// 							temp[5] = p.z + tank.DZ*p.ivz;
// 							temp[6] = p.E;
// 							fwrite(temp, sizeof(double), 7, fp);
// 							nStepDebug++;
// 							if (nStepDebug >= STEP_NUM)
// 							{
// 								fclose(fp);
// 								exitApp("Step debug finished!");
// 							}
// 
// #endif
							DE = mat[0].knock(it); //knock and loss some energy

							if (DE != 0)
							{
								tank.deposit(p, DE); //record the energy loss in the phantom lattice
#ifdef STEP_DEBUG_STATUS
								fprintf(fp, "%7d\t%.8g\n", (int)p.iabsv, DE);
								nStepDebug++;
								if (nStepDebug >= STEP_NUM)
								{
									fclose(fp);
									exitApp("Step debug finished!");
								}
#endif		
							}

							if (p.E == 0) cross = true; //energy too low, same as escape
						}

						if (cross) //we need to stop simulating current particle,
						{
							if (secStack.empty()) break; //got no particle at all, so we stop the while loop
							//still secondary particles, fetch the top one, and continue the while loop (starting a new jump)
							p = secStack.top();
							secStack.pop();
#ifdef STEP_DEBUG_STATUS
							fprintf(fp, "%7d\t%.8g\n", (int)p.iabsv, -p.E);
							nStepDebug++;
							if (nStepDebug >= STEP_NUM)
							{
								fclose(fp);
								exitApp("Step debug finished!");
							}
#endif

							tank.deposit(p, -p.E);

							if (p.E<1.00000001*ELow || p.E>0.99999999*EUp) exitApp("E out of range in function start()");
							Mode = SOFT;
							continue;
						}
						//do nothing, continue current particle's simulation
					}
				}
#ifdef ACCU_DEBUG
				++nIncident;
				if (nIncident >= INCIDENT_NUM) break;
#endif
			}

			if (0 == it && hist > percent*percentHist)
			{
				Log("Thread 0: processed ------------------------ %3d%%,   speed = %d h/s\n", percent,
					int(hist/ percentSpeed.stop(true)));
				++percent;
			}
#ifdef ACCU_DEBUG
			if (nIncident >= INCIDENT_NUM)
			{
				fSIMU = hist+1;
				break;
			}
#endif
		}
	}

	double norm = 1.0889e15 * SourceHead_BeamOnTime();
	
	tank.dose;
	tank.addDose(fSIMU, norm);

	Log("Thread 0: processed ------------------------ 100%%, <speed> = %d h/s\n", int(fSIMU / NThread / rc.stop(true)));
	Log("\nActual history number on this node is %.0f due to thread division\n", fSIMU);
	Log("\nThe simulation costs %f seconds, ie %f minutes\n", rc.stop(true), rc.stop(true) / 60);

	/*****************************************************************************************************/


	string outname;
	cf.getValue("output file name", outname);
	outname += ".dose";
	tank.output(outname.c_str());
	Log("\nThe total run time from entering main() is %f s\n\n\n\n", totalTime.stop());

	/*******************clean up work for this run***********************/
	SourceHead_Delete();
	delete[] tp;
	tp = NULL;
}

