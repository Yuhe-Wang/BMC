#ifdef WIN32
//windows specific, hide the cmd window
#pragma comment(linker, "/subsystem:windows /entry:mainCRTStartup")
#endif
#include "../Tools/Tools.h"
#include "../SourceHead/SourceHead.h"

double LeafAdjustment = 0;
void exitApp(const char *inf)
{
	printf("fatal error: %s", inf);
	getchar();
	exit(-1);
}

class SEGMENT
{
public:
	SEGMENT();
	SEGMENT(ConfigFile *cf);
	~SEGMENT()
	{
		if (primAT) delete[] primAT;
		if (scatAT) delete[] scatAT;
		if (primMask) delete[] primMask;
		if (scatMask) delete[] scatMask;
	}
	bool load(ConfigFile *cf);
	bool load(const string& vrs); //load segment from ViewRay's format
	bool setPos(vector<double> &values);
	
	//static constant about the MLC
	static const int NLeaf;
	static const double LeafWidth;
	static const double DX;

	vector<pair<double, double> >  leafPositions;  //!< List of left and right leaf positions in cm in the isocenter plane for all leaf pairs
	double                         onTime;         //!< Time (or relative fluence) for this segment

	unsigned short*    primAT;        //!< Alias sampling table for primary photons (used to sample x- and y-position in the plane above the MLC)
	unsigned short*    scatAT;        //!< Same as _primAT but for scattered particles
	char*              primMask;      //!< The primary particles mask 
	char*              scatMask;      //!< The scattered particles mask
	double             nSample;       //!< Number of times to sample
	double             pPrim;         //!< Probability for primary.
	int                beamID;        //!< which beam it belongs to. The ID starts from zero
	double             fluence;       //!< relative fluence during the treatment, used to sample the segment
	// Mask = 0 means within open+small margin
	// Mask = 1 means between small and large margin
	// Mask = 2 means outside large margin
};

const int     SEGMENT::NLeaf = 30;
const double  SEGMENT::LeafWidth = 1.05; //unit cm
const double  SEGMENT::DX = NLeaf*LeafWidth; //unit cm

SEGMENT::SEGMENT() :primAT(NULL), scatAT(NULL), primMask(NULL), scatMask(NULL)
{
	leafPositions.resize(NLeaf);
	for (int i = 0; i < NLeaf; ++i) //set all the leaves in close state
	{
		leafPositions[i] = make_pair(0, 0);
	}
}
SEGMENT::SEGMENT(ConfigFile *cf) :primAT(NULL), scatAT(NULL), primMask(NULL), scatMask(NULL)
{
	leafPositions.resize(NLeaf);
	for (int i = 0; i < NLeaf; ++i) //set all the leaves in close state
	{
		leafPositions[i] = make_pair(0, 0);
	}
	load(cf);
}
bool SEGMENT::load(ConfigFile *cf) //true means loaded successfully
{
	vector<double> values;
	if (!cf->getValue("on time", onTime)) return false;
	if (!cf->getValue("leaf position", values)) return false;
	bool status = setPos(values);
	if (!status) return status;
	while (cf->getValue("leaf position", values, true)) //if we can find more overlap definition
	{
		status = setPos(values);
	}
	return status;
}
bool SEGMENT::setPos(vector<double> &values)
{
	if (values.size() != 4) exitApp("The number of leaf position parameters isn't 4!");
	int start = int(values[0]);
	int end = int(values[1]);
	if (start<0 || end>NLeaf - 1) return false;
	if (values[3] - values[2] + LeafAdjustment >= 0) //may apply an overall adjustment on leaf distance
	{
		values[2] -= LeafAdjustment*0.5;
		values[3] += LeafAdjustment*0.5;
	}
	pair<double, double> pos = make_pair(values[2], values[3]);
	for (int i = start; i <= end; ++i) leafPositions[i] = pos;
	return true;
}

class BEAM //describing the beam configuration
{
public:
	BEAM(){};
	BEAM(ConfigFile *cf){ load(cf); }
	bool load(ConfigFile *cf);
	bool load(const string& vrb);//load beam configuration from ViewRay's plan
	double beamOnTime()
	{
		double sum = 0;
		for (size_t i = 0; i < segments.size(); ++i)
		{
			sum += segments[i].onTime;
		}
		return sum;
	}
	vector<SEGMENT> segments;     //!< The list of segments that belong to this beam
	double  gantryAngle;          //!< The gantry angle in degree
	double  sinPhi, cosPhi;       //!< sin and cos of the gantry angle to accelerate calculation
	// 	double  Xiso;                 //!< The x-position of the isocenter in cm
	// 	double  Yiso;                 //!< The y-position of the isocenter in cm
	// 	double  Ziso;                 //!< The z-position of the isocenter in cm
	int     headIndex;
};
bool BEAM::load(ConfigFile *cf)
{
	segments.resize(0);
	vector<double> values;
	if (!cf->getValue("gantry angle", gantryAngle)) return false;
	cosPhi = cos(gantryAngle*PI / 180);
	sinPhi = sin(gantryAngle*PI / 180);
	if (!cf->getValue("head index", headIndex)) return false;

	ConfigFile *seg_cf;
	cf->resetSearchIndex();
	for (int NSeg = 1;; ++NSeg)
	{
		seg_cf = cf->getBlock("segment", true);
		if (NULL == seg_cf)
		{
			if (1 == NSeg) return false; //cannot even find the first material
			else break; //at least loaded one material
		}
		segments.push_back(SEGMENT(seg_cf));
	}

	return true;
}


int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		printf("This program converts the beam configuration in text format(.py or PlanOverview.txt) to ViewRay's binary format(.beams)\n\n");
		printf("Usage: %s oldConfigFileName newConfigFileName\n\n", argv[0]);
		printf("You can also execute by dragging the text configFile to this program's icon, and the default output name would be mfQAData.beams\n\n");
		printf("Press enter to exit...");
		getchar();
		return 0;
	}

	ConfigFile cff;
	if(!cff.parse(argv[1]))
	{
		//try to parse the PlanOverview.txt
		const char* output = argc >= 3 ? argv[2] : "mfQAData.beams";
		if (SourceHead_ConvertBeams(argv[1], output)) return 0;
		printf("cannot parse beams file!");
		getchar();
		return -1;
	}
	ConfigFile* cf = cff.getBlock("SOURCEHEAD");
	if (!cf->getValue("Leaf spacing adjustment", LeafAdjustment)) LeafAdjustment = -0.07; //cm
	vector<BEAM> beams;
	vector<double> values;
	if (!cf->getValue("isocenter", values) || values.size() != 3) return false;
	double Isocenter_X = values[0];
	double Isocenter_Y = values[1];
	double Isocenter_Z = values[2];
	//need to calculate the isocenter if customized phantom is used
	bool b_fromFile = true;
	ConfigFile* pcf = cff.getBlock("PHANTOM");
	if (pcf)
	{
		string fname;
		if (!pcf->getValue("phantomName", fname)) b_fromFile = false; //no phantom file specified
		if (!b_fromFile)
		{
			pcf->getValue("DX DY DZ", values);
			double DY = values[1];
			pcf->resetSearchIndex();
			pcf->getValue("matrix", values, true);
			if (values.size() < 7) exitApp("Incorrect custom phantom definition");
			int NY = int(values[1]);
			Isocenter_X = -values[3];
			Isocenter_Y = -values[4];
			Isocenter_Z = -values[5];
			double SSD = 0;
			if (pcf->getValue("SSD", SSD)) Isocenter_Y = SAD - SSD - NY*DY / 2;
			vector<double> centerXZ;
			if (pcf->getValue("centerXZ", centerXZ))
			{
				Isocenter_X = -centerXZ[0];
				Isocenter_Z = -centerXZ[1];
			}
		}
	}
	
	//search for the definition of beams
	ConfigFile *beamCF;
	cf->resetSearchIndex();
	for (int NBeam = 1;; ++NBeam)
	{
		beamCF = cf->getBlock("Beam", true);
		if (NULL == beamCF)
		{
			if (1 == NBeam) exitApp("cannot load any beam parameter!"); //cannot even find the first material
			else  break;//at least loaded one beam
		}
		beams.resize(NBeam);
		beams[NBeam - 1].load(beamCF);
	}

	//output the beams in ViewRay's format
	
	//make sure current directory is where the config file locates
	fs::path p = fs::path(argv[1]).parent_path() / fs::path("mfQAData.beams");
	FILE* fp = fopen(p.string().c_str(), "wb");
	if (NULL == fp) return -1;

	unsigned int Major_Header = 1, Minor_Header = 0, Num_Of_Beams, Head, Num_Of_Segments;
	double GantryAngle, Beam_On_Time, Leaf_Position_Left, Leaf_Position_Right;
	Num_Of_Beams = beams.size();
	fwrite(&Major_Header, sizeof(unsigned int), 1, fp);
	fwrite(&Minor_Header, sizeof(unsigned int), 1, fp);
	fwrite(&Num_Of_Beams, sizeof(unsigned int), 1, fp);
	for (unsigned int i = 0; i < Num_Of_Beams; ++i) //for each beam
	{
		GantryAngle = beams[i].gantryAngle;
		Head = beams[i].headIndex;
		Num_Of_Segments = beams[i].segments.size();
		fwrite(&GantryAngle, sizeof(double), 1, fp);
		fwrite(&Isocenter_X, sizeof(double), 1, fp);
		fwrite(&Isocenter_Y, sizeof(double), 1, fp);
		fwrite(&Isocenter_Z, sizeof(double), 1, fp);
		fwrite(&Head, sizeof(unsigned int), 1, fp);
		fwrite(&Num_Of_Segments, sizeof(unsigned int), 1, fp);
		for (unsigned int j = 0; j < Num_Of_Segments; ++j)
		{
			Beam_On_Time = beams[i].segments[j].onTime;
			fwrite(&Beam_On_Time, sizeof(double), 1, fp);
			for (int k = 0; k < 30; ++k)
			{
				Leaf_Position_Left = beams[i].segments[j].leafPositions[k].first;
				Leaf_Position_Right = beams[i].segments[j].leafPositions[k].second;
				fwrite(&Leaf_Position_Left, sizeof(double), 1, fp);
				fwrite(&Leaf_Position_Right, sizeof(double), 1, fp);
			}
		}
	}

	fclose(fp);
	return 0;

}