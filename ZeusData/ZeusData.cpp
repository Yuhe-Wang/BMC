#include "ZeusData.h"
#include "zeusCrossSection.h"
#include "zeusSpline.h"


bool ZeusData_load(const char* dataFolder) {

	string aDataFolder(dataFolder);
	//
	// Load the material data
	//
	if (Zeus::CrossSectionData::Material::loadData(aDataFolder)) return false;
	Zeus::CrossSectionData::Material _material;

	//
	// Load the various pre-computed data files 
	//
	if (Zeus::CrossSectionData::BremmConstants::loadData(aDataFolder)) return false;
	if (Zeus::CrossSectionData::LambdaBr::loadData(aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::LambdaCo::loadData(aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::LambdaMoller::loadData(aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::LambdaPhoton::loadData(aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::LambdaPair::loadData(aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::Kerma::loadData(aDataFolder, _material.numMaterials())) return false;
	//if ( Zeus::CrossSectionData::QSurface::loadData( true, aDataFolder, _material.numMaterials() ) ) return false;
	if (Zeus::CrossSectionData::QSurface::loadData(false, aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::RestrictedStoppingPower::loadData(aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::ScreeningParameter::loadData(aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::Range::loadData(aDataFolder, _material.numMaterials())) return false;
	if (Zeus::CrossSectionData::InverseRange::loadData(aDataFolder, _material.numMaterials())) return false;

	//
	// Prepare LambdaTot, InvLambdaTot and BremsProbability tables
	//
	int nmat = _material.numMaterials();
	vector<double*> lamTot(nmat);
	vector<double*> pbrem(nmat);
	vector<double*> invLamTot(nmat);
	double Eabs = _material.getEabs();
	double Emax = _material.getEmax();
	int nbin = 1024;
	double dE = (Emax - Eabs) / (nbin - 1);
	vector<double> eminArray(nmat, Eabs);
	vector<double> emaxArray(nmat, Emax);
	vector<double> lamMinArray(nmat, 0);
	vector<double> lamMaxArray(nmat);
	double *eners = new double[nbin];
	double *spla = new double[nbin], *splb = new double[nbin], *splc = new double[nbin], *spld = new double[nbin];
	for (int imat = 0; imat < nmat; ++imat) {
		//InformationLog("\n******** Working on material %d\n",imat);
		double *ltot = new double[nbin];
		double *pbr = new double[nbin];
		double *invlam = new double[nbin];
		double sumLam = 0;
		//InformationLog("+++ lambda\n");
		for (int i = 0; i < nbin - 1; ++i) {
			double E = Eabs + (0.5 + i)*dE;
			double L = Zeus::CrossSectionData::RestrictedStoppingPower::stpwr(E, imat);
			double lammo = 1 / Zeus::CrossSectionData::LambdaMoller::lmbdamo(E, imat);
			double lambr = Zeus::CrossSectionData::LambdaBr::lmbdabr(E, imat);
			double lamtot = lammo + lambr;
			double Li = dE / L;
			sumLam += lamtot*Li;
			//InformationLog("%g     %g  %g  %g\n",1e-6*(E+0.5*dE),sumLam,lammo,lambr);
			ltot[i + 1] = sumLam;
		}
		ltot[0] = 0;
		for (int i = 0; i < nbin; ++i) {
			double E = Eabs + dE*i;
			eners[i] = E;
			double lammo = 1 / Zeus::CrossSectionData::LambdaMoller::lmbdamo(E, imat);
			double lambr = Zeus::CrossSectionData::LambdaBr::lmbdabr(E, imat);
			double lamtot = lammo + lambr;
			pbr[i] = lamtot > 0 ? lambr / lamtot : 1;
		}
		lamTot[imat] = ltot;
		pbrem[imat] = pbr;
		invLamTot[imat] = invlam;

		Zeus::Spline::prepareSpline(ltot, eners, spla, splb, splc, spld, 0, 0, nbin);
		lamMaxArray[imat] = ltot[nbin - 1];
		double delta = ltot[nbin - 1] / (nbin - 1);
		invlam[0] = Eabs;
		//InformationLog("+++ inverse lambda\n");
		for (int i = 1; i < nbin; ++i) {
			double lam = delta*i;
			if (lam > 0.999999*ltot[nbin - 1]) lam = 0.999999*ltot[nbin - 1];
			invlam[i] = Zeus::Spline::interpolate(lam, nbin, ltot, spla, splb, splc, spld);
			//InformationLog("%g     %g\n",delta*i,invlam[i]);
		}

	}

	Zeus::CrossSectionData::LambdaTot::prepareData(eminArray, emaxArray, nbin, lamTot);
	Zeus::CrossSectionData::InvLambdaTot::prepareData(lamMinArray, lamMaxArray, nbin, invLamTot);
	Zeus::CrossSectionData::BremsProbability::prepareData(eminArray, emaxArray, nbin, pbrem);

	for (int imat = 0; imat < nmat; ++imat) {
		delete[] lamTot[imat];
		delete[] pbrem[imat];
		delete[] invLamTot[imat];
	}

	delete[] spla;
	delete[] splb;
	delete[] splc;
	delete[] spld;

	return true;

};

void ZeusData_unload() {

	Zeus::CrossSectionData::Material::unloadData();
	Zeus::CrossSectionData::BremmConstants::unloadData();
	Zeus::CrossSectionData::LambdaBr::unloadData();
	Zeus::CrossSectionData::LambdaCo::unloadData();
	Zeus::CrossSectionData::LambdaMoller::unloadData();
	Zeus::CrossSectionData::LambdaPhoton::unloadData();
	Zeus::CrossSectionData::LambdaPair::unloadData();
	Zeus::CrossSectionData::Kerma::unloadData();
	Zeus::CrossSectionData::QSurface::unloadData();
	Zeus::CrossSectionData::RestrictedStoppingPower::unloadData();
	Zeus::CrossSectionData::ScreeningParameter::unloadData();
	Zeus::CrossSectionData::Range::unloadData();
	Zeus::CrossSectionData::InverseRange::unloadData();
	Zeus::CrossSectionData::LambdaTot::unloadData();
	Zeus::CrossSectionData::InvLambdaTot::unloadData();
	Zeus::CrossSectionData::BremsProbability::unloadData();
}

void ZeusData_lambdaWoodcock(int aNumElements, const float *dens, const char *mat, int nvoxels, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b)
{
	Zeus::CrossSectionData::LambdaPhoton lmdphoton;
	lmdphoton.SetupWoodcock(aNumElements, dens, mat, nvoxels);
	lmdphoton.OutputWoodcock(aIndex, bIndex, a, b);
}

void ZeusData_BremConst(int matid, double f0[3], double& bcb)
{
	auto p = Zeus::CrossSectionData::BremmConstants::getData();
	f0[0] = p->fo(matid)[0];
	f0[1] = p->fo(matid)[1];
	f0[2] = p->fo(matid)[2];
	bcb = p->bcb()[matid];
}

void ZeusData_lambdaPhoton(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d)
{
	Zeus::CrossSectionData::LambdaPhoton::getData()->outputSpline(matid, aIndex, bIndex, a, b, c, d);
}

void ZeusData_lambdaCo(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d)
{
	Zeus::CrossSectionData::LambdaCo::getData()->outputSpline(matid, aIndex, bIndex, a, b, c, d);
}

void ZeusData_lmbdapp(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d)
{
	Zeus::CrossSectionData::LambdaPair::getData()->outputSpline(matid, aIndex, bIndex, a, b, c, d);
}

void ZeusData_ScreeningParameter(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d)
{
	Zeus::CrossSectionData::ScreeningParameter::getData()->outputSpline(matid, aIndex, bIndex, a, b, c, d);
}

void ZeusData_QSurface(int matid, double& fq, double& ieq0, vector<vector<double>>& qs)
{
	Zeus::CrossSectionData::QSurface::getData()->Output(matid, fq, ieq0, qs);
}

void ZeusData_Range(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b)
{
	Zeus::CrossSectionData::Range::getData()->outputLinear(matid, aIndex, bIndex, a, b);
}

void ZeusData_InverseRange(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b)
{
	Zeus::CrossSectionData::InverseRange::getData()->outputLinear(matid, aIndex, bIndex, a, b);
}