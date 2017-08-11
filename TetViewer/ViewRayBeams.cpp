#include "PreCompile.h"
#include "ViewRayBeams.h"

const int SEGMENT::NLeaf = 30;

bool SEGMENT::loadViewRaySegment(const string& vrs) //load segment from ViewRay's format
{
	onTime = -1;
	size_t start = 0;
	for (int i = 0; i < 3; ++i)	start = vrs.find("\n", start + 1);
	sscanf(vrs.c_str() + start, "     Total Beam-On Time (s): %lf", &onTime);
	if (onTime < 0)
	{
		printf("Error: cannot read the beam on time for the segment.\n");
		return false;
	}
	for (int i = 0; i < 3; ++i)	start = vrs.find("\n", start + 1);
	double left = -1, right = -1;
	int sr = -100;
	for (int i = 0; i < NLeaf; ++i)
	{
		sscanf(vrs.c_str() + start, "      %*d      %lf                %lf\n", &left, &right);
		start = vrs.find("\n", start + 1); //update the search header index
		leafPositions[i] = make_pair(left, right);
	}
	return true;
}

bool BEAM::loadViewRayBeam(const string& vrb)//load beam configuration from ViewRay's plan
{
	//find beam index
	int angle = 0;
	sscanf(vrb.c_str(), "Beam %d:\n    Angle: %d¡ã", &headIndex, &angle);
	if (angle < 0 || angle>360)
	{
		printf("Error: incorrect beam angle!\n");
		return false;
	}
	gantryAngle = angle;

	segments.resize(0);//make sure previous segments are released
	//search the definition of each segment
	vector<string> vrSeg;
	size_t spos = vrb.find("Segment Id");
	size_t epos;
	do
	{
		epos = vrb.find("Segment Id", spos + 1);
		vrSeg.push_back(vrb.substr(spos, epos - spos));
		spos = epos;
	} while (epos != string::npos);

	segments.resize(vrSeg.size());
	for (size_t i = 0; i < vrSeg.size(); ++i) segments[i].loadViewRaySegment(vrSeg[i]);

	return true;
}

bool ViewRayBeams::getViewRayBeams(const string& file) // get beam definition from ViewRay's plan
{
	//let's find the file extension
	string ext = file.substr(file.find_last_of("."), string::npos);
	if (ext == ".txt") // use ViewRay plan overview, which may change due to ViewRay's update
	{
		FILE *fp = fopen(file.c_str(), "r");
		if (NULL == fp)
		{
			printf("Cannot open %s. Will try customized beam configuration.\nPress enter to continue...\n", file.c_str());
			getchar();
			return false;
		}
		string vrPlan;
		char ch = 0;
		do
		{
			ch = fgetc(fp);
			vrPlan += ch;
		} while (EOF != ch);
		fclose(fp);

		size_t spos = 0, epos = 0;

		//find the ISO center
// 		spos = vrPlan.find("Isocenter Coordinate (cm)");
// 		sscanf(vrPlan.c_str() + spos, "Isocenter Coordinate (cm): (%lf,%lf,%lf)", &isoX, &isoY, &isoZ);

		//find the prescription dose and fraction number
		spos = vrPlan.find("Total Prescription Dose");
		spos = vrPlan.find(":", spos) + 1;
		sscanf(vrPlan.c_str() + spos, "%lf", &stdValue);

		spos = vrPlan.find("Prescription Dose Per Fraction");
		spos = vrPlan.find(":", spos) + 1;
		double fractionDose;
		sscanf(vrPlan.c_str() + spos, "%lf", &fractionDose);
		treatmentFraction = int(stdValue / fractionDose);

		//find the ISO center
		spos = vrPlan.find("Isocenter Coordinate (cm)");
		sscanf(vrPlan.c_str() + spos, "Isocenter Coordinate (cm): (%lf,%lf,%lf)", &isoX, &isoY, &isoZ);
		//divide the setting to groups
		vector<string> vrGroup;
		spos = vrPlan.find("Planning Beam Angle Group");
		do
		{
			epos = vrPlan.find("Planning Beam Angle Group", spos + 1);
			vrGroup.push_back(vrPlan.substr(spos, epos - spos));
			spos = epos;
		} while (epos != string::npos);

		vector<string> vrBeam;
		//divide each group to beams
		for (size_t i = 0; i < vrGroup.size(); ++i)
		{
			size_t b1 = vrGroup[i].find("Beam 1");
			size_t b2 = vrGroup[i].find("Beam 2");
			size_t b3 = vrGroup[i].find("Beam 3");
			string sb1 = vrGroup[i].substr(b1, b2 - b1);
			string sb2 = vrGroup[i].substr(b2, b3 - b2);
			string sb3 = vrGroup[i].substr(b3, string::npos - b3);
			if (string::npos == sb1.find("Angle: N/A")) vrBeam.push_back(sb1);
			if (string::npos == sb2.find("Angle: N/A")) vrBeam.push_back(sb2);
			if (string::npos == sb3.find("Angle: N/A")) vrBeam.push_back(sb3);
		}

		beams.resize(vrBeam.size());
		for (size_t i = 0; i < beams.size(); ++i)
		{
			beams[i].loadViewRayBeam(vrBeam[i]);
		}
		getFlux();
		//need to change the coordinates of ISO center. This may vary depending on the verison of ViewRay software
		double isoT = isoY;
		isoY = -isoZ;
		isoZ = isoT;
		return true;
	}
	else if (ext == ".beams")
	{
		FILE *fp = fopen(file.c_str(), "rb");
		if (NULL == fp)
		{
			printf("Cannot open %s. Will try customized beam configuration.\nPress enter to continue...\n", file.c_str());
			getchar();
			return false;
		}

		unsigned int Major_Header, Minor_Header, Num_Of_Beams, Head, Num_Of_Segments;
		double GantryAngle, Isocenter_X, Isocenter_Y, Isocenter_Z, Beam_On_Time, Leaf_Position_Left, Leaf_Position_Right;
		fread(&Major_Header, sizeof(unsigned int), 1, fp);
		fread(&Minor_Header, sizeof(unsigned int), 1, fp);
		fread(&Num_Of_Beams, sizeof(unsigned int), 1, fp);
		beams.resize(Num_Of_Beams);
		for (unsigned int i = 0; i < Num_Of_Beams; ++i) //for each beam
		{
			fread(&GantryAngle, sizeof(double), 1, fp);
			fread(&Isocenter_X, sizeof(double), 1, fp);
			fread(&Isocenter_Y, sizeof(double), 1, fp);
			fread(&Isocenter_Z, sizeof(double), 1, fp);
			fread(&Head, sizeof(unsigned int), 1, fp);
			fread(&Num_Of_Segments, sizeof(unsigned int), 1, fp);
			beams[i].gantryAngle = GantryAngle;
			beams[i].headIndex = Head;
			beams[i].segments.resize(Num_Of_Segments);
			for (unsigned int j = 0; j < Num_Of_Segments; ++j)
			{
				fread(&Beam_On_Time, sizeof(double), 1, fp);
				beams[i].segments[j].onTime = Beam_On_Time;
				for (int k = 0; k < SEGMENT::NLeaf; ++k)
				{
					fread(&Leaf_Position_Left, sizeof(double), 1, fp);
					fread(&Leaf_Position_Right, sizeof(double), 1, fp);
					beams[i].segments[j].leafPositions[k] = make_pair(Leaf_Position_Left, Leaf_Position_Right);
				}
			}
		}
		isoX = Isocenter_X;
		isoY = Isocenter_Y;
		isoZ = Isocenter_Z;
		fclose(fp);
		getFlux();
		return true;
	}

	return false;
}