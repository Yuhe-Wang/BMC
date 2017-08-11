#ifndef _VIEWRAYBEAMS_H_
#define _VIEWRAYBEAMS_H_

class SEGMENT
{
public:
	SEGMENT()
	{
		leafPositions.resize(NLeaf);
		for (int i = 0; i < NLeaf; ++i) //set all the leaves in close state
		{
			leafPositions[i] = make_pair(0, 0);
		}
	}
	double getFlux()
	{
		double sum = 0;
		for (int i = 0; i < NLeaf; ++i)
		{
			double width = leafPositions[i].second - leafPositions[i].first;
			if (width < 0) width = 0;
			sum += width;
		}
		return sum*1.05*onTime;
	}
	bool loadViewRaySegment(const string& vrs); //load segment from ViewRay's format
	vector<pair<double, double> >  leafPositions;  //!< List of left and right leaf positions in cm in the isocenter plane for all leaf pairs
	double                         onTime;         //!< Time (or relative fluence) for this segment

	static const int NLeaf;
};

class BEAM //describing the beam configuration
{
public:
	bool loadViewRayBeam(const string& vrb);//load beam configuration from ViewRay's plan
	double beamOnTime()
	{
		double sum = 0;
		for (size_t i = 0; i < segments.size(); ++i)
		{
			sum += segments[i].onTime;
		}
		return sum;
	}
	double getFlux()
	{
		double sum = 0;
		for (size_t i = 0; i < segments.size(); ++i)
		{
			sum += segments[i].getFlux();
		}
		flux = sum;
		return sum;
	}
	vector<SEGMENT> segments;     //!< The list of segments that belong to this beam
	double  gantryAngle;          //!< The gantry angle in degree
	unsigned int headIndex;
	double flux;
};

class ViewRayBeams
{
public:
	//method
	ViewRayBeams() : isoX(0), isoY(0), isoZ(0), stdValue(0), treatmentFraction(1){}

	bool init(ConfigFile *cf);
	vector<BEAM> beams;
	double beamOnTime()
	{
		double sum = 0;
		for (size_t i = 0; i < beams.size(); ++i)
		{
			sum += beams[i].beamOnTime();
		}
		return sum;
	}
	void getFlux()
	{
		int nb = beams.size();
		for (int i = 0; i < nb; ++i) beams[i].getFlux();
	}
	void getMaxMinFlux(double& maxv, double& minv)
	{
		if (beams.size() > 0)
		{
			maxv = minv = beams[0].flux;
			int nb = beams.size();
			for (int i = 0; i < nb; ++i)
			{
				if (beams[i].flux>maxv) maxv = beams[i].flux;
				if (beams[i].flux<minv) minv = beams[i].flux;
			}
		}
		else maxv = minv = 0;
	}
	void getISO(double& x, double& y, double& z) { x = isoX, y = isoY, z = isoZ; }
	bool getViewRayBeams(const string& file);//get beam config from ViewRay's plan file
private:

	//data
	double isoX, isoY, isoZ;
	double stdValue;
	int treatmentFraction;
};

#endif