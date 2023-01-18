
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
//Contains public member functions that allow for easy access to the binning information of various useful quantities for the CaFe experiment.
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

#ifndef Bins_H
#define Bins_H
using namespace std;

class Bins
{
	public:
	//static Double_t Variable(int i)	{Double_t Variable[] = {nbins_MF, xmin_MF, xmax_MF, nbins_SRC, xmin_SRC, xmax_SRC}; return Variable[i];}
	
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//------------------------------------------------------------------------PID Histograms Bins-----------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static Double_t coin(int i)			{ Double_t coin[] = { 200, -20, 20, 200, -20, 20 };			return coin[i];}		//coincidence time [ns]
	//----HMS DETECTORS (FOR PID / TRACKING EFF. CHECKS)----
	static Double_t hcer(int i)			{ Double_t hcer[] = { 100, 0.001, 15, 100, 0.001, 15 };		return hcer[i];}		//cherenkov NPE sum
	static Double_t hcal(int i)			{ Double_t hcal[] = { 100, 0.001, 1.5, 100, 0.001, 1.5 };	return hcal[i];}		//calorimeter etotnorm / etottracknorm
	static Double_t hbeta(int i)		{ Double_t hbeta[] = { 100, 0.4, 1.6, 100, 0.4, 1.6 };		return hbeta[i];}		//calculated beta
	//----SHMS DETECTORS (FOR PID / TRACKING EFF. CHECKS)----
	static Double_t pngcer(int i)		{ Double_t pngcer[] = { 100, 0.001, 25, 100, 0.001, 25 };	return pngcer[i];}		//noble gas cherenkov npe sum
	static Double_t phgcer(int i)		{ Double_t phgcer[] = { 100, 0.001, 25, 100, 0.001, 25 };	return phgcer[i];}		//heavy gas cherenkob npe sum
	static Double_t pcal(int i)			{ Double_t pcal[] = { 100, 0.001, 1.5, 100, 0.001, 1.5 };	return pcal[i];}		//calorimeter etotnorm / etottracknorm
	static Double_t pbeta(int i)		{ Double_t pbeta[] = { 100, 0.4, 1.6, 100, 0.4, 1.6 };		return pbeta[i];}		//calculated beta

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-------------------------------------------------------------Primary Kinematics (electron kinematics)-------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static Double_t the(int i)			{ Double_t the[] = { 100, 0, 20, 100, 0, 20 };				return the[i]; }		//electron arm central angle [Deg]
	static Double_t W(int i)			{ Double_t W[] = { 100, 0, 2, 100, 0, 2 };					return W[i]; }			//invariant mass [GeV]
	static Double_t W2(int i)			{ Double_t W2[] = { 100, 0, 4, 100, 0, 4 };					return W2[i]; }			//invariant mass squared [(GeV)^2]
	static Double_t Q2(int i)			{ Double_t Q2[] = { 100, 1, 3, 100, 1, 3 };					return Q2[i];}			//4-momentum transfer [GeV]
	static Double_t X(int i)			{ Double_t X[] = { 100, 0, 2, 100, 0, 2 };					return X[i]; }			//Bjorken-X
	static Double_t nu(int i)			{ Double_t nu[] = { 100, 0.5, 2.5, 100, 0.5, 2.5 };			return nu[i]; }			//energy transfer [GeV]
	static Double_t q(int i)			{ Double_t q[] = { 100, 1, 3, 100, 1, 3 };					return q[i]; }			//3-momentum transfer [GeV]
	static Double_t qx(int i)			{ Double_t qx[] = { 100, -2, 2, 100, -2, 2 };				return qx[i]; }			//3-momentum transfer (x-component) [GeV]
	static Double_t qy(int i)			{ Double_t qy[] = { 100, -2, 2, 100, -2, 2 };				return qy[i]; }			//3-momentum transfer (y-component) [GeV]
	static Double_t qz(int i)			{ Double_t qz[] = { 100, 0, 3, 100, 0, 3 };					return qz[i]; }			//3-momentum transfer (z-component) [GeV]
	static Double_t thq(int i)			{ Double_t thq[] = { 100, 0, 80, 100, 0, 80 };				return thq[i]; }		//in-plane angle between q-vector and +z (beam) [Deg]
	static Double_t phq(int i)			{ Double_t phq[] = { 100, -200, 200, 100, -200, 200 };		return phq[i]; }		//out-of-plane angle between q-vector and +z (beam) [Deg]

	
	static Double_t epsilon(int i)		{ Double_t epsilon[] = { 100, -0.1, 1.5, 100, -0.1, 1.5 };	return epsilon[i];}		//virtual photon polarization transfer
	static Double_t omega(int i)		{ Double_t omega[] = { 100, 0, 6, 100, 0, 6 };				return omega[i];}		//????????????
	static Double_t kf(int i)			{ Double_t kf[] = { 100, 6, 11, 100, 8, 11 };				return kf[i];}			//final electron arm momentum [GeV]
	static Double_t ki(int i)			{ Double_t ki[] = { 100, 0, 11, 100, 0, 11 };				return ki[i]; }			//initial electron arm momentum [GeV]

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-------------------------------------------------------------------Secondary Kinematics (Hadron)------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static Double_t Em(int i)			{ Double_t Em[] = { 100, -0.1, 0.5, 100, -0.1, 0.55 };		return Em[i]; }			//Missing Energy [GeV]
	static Double_t Em_nuc(int i)		{ Double_t Em_nuc[] = { 100, -0.1, 0.5, 100, -0.1, 0.55 };	return Em_nuc[i]; }		//Nuclear Missing Energy [GeV]
	static Double_t Pm(int i)			{ Double_t Pm[] = { 100, -0.1, 2, 100, -0.1, 2 };			return Pm[i]; }			//Missing (Recoil) Momentum [GeV/c]
	static Double_t Pmx_lab(int i)		{ Double_t Pmx_lab[] = { 100, -1, 1, 100, -1, 1 };			return Pmx_lab[i]; }	//Missing Momentum (x-component, lab) [GeV]
	static Double_t Pmy_lab(int i)		{ Double_t Pmy_lab[] = { 100, -1, 1, 100, -1, 1 };			return Pmy_lab[i]; }	//Missing Momentum (y-component, lab) [GeV]
	static Double_t Pmz_lab(int i)		{ Double_t Pmz_lab[] = { 100, -0.1, 2, 100, -0.1, 2 };		return Pmz_lab[i]; }	//Missing Momentum (z-component, lab) [GeV]
	static Double_t Pmx_q(int i)		{ Double_t Pmx_q[] = { 100, -2, 2, 100, -2, 2 };			return Pmx_q[i]; }		//Missing Momentum (x-component, q-frame: +z_lab rotated to +q) [GeV]
	static Double_t Pmy_q(int i)		{ Double_t Pmy_q[] = { 100, -2, 2, 100, -2, 2 };			return Pmy_q[i]; }		//Missing Momentum (y-component, q-frame: +z_lab rotated to +q) [GeV]
	static Double_t Pmz_q(int i)		{ Double_t Pmz_q[] = { 100, -2, 2, 100, -2, 2 };			return Pmz_q[i]; }		//Missing Momentum (z-component, q-frame: +z_lab rotated to +q) [GeV]
	static Double_t Tx(int i)			{ Double_t Tx[] = { 100, -0.1, 1.5, 100, -0.1, 1 };			return Tx[i]; }			//Kinetic Energy (Detected Hadron, X) [GeV]
	static Double_t Tr(int i)			{ Double_t Tr[] = { 100, -0.1, 0.5, 100, -0.1, 0.5 };		return Tr[i]; }			//Kinetic Energy (Undetected Recoil System, R) [GeV]
	static Double_t MM(int i)			{ Double_t MM[] = { 200, -0.2, 1.2, 200, -0.2, 1.2 };		return MM[i]; }			//Missing Mass (Undetected Recoil System Mass) [GeV]
	static Double_t thxq(int i)			{ Double_t thxq[] = { 100, -10, 40, 100, -10, 40 };			return thxq[i]; }		//in-plane angle between detected particle and q [Deg]
	static Double_t thrq(int i)			{ Double_t thrq[] = { 100, -20, 200, 100, -20, 200 };		return thrq[i]; }		//In-plane angle between the recoil system and q [Deg]
	static Double_t phxq(int i)			{ Double_t phxq[] = { 100, -200, 200, 100, -200, 200 };		return phxq[i]; }		//Out-of-plane angle between detected particle and q [Deg]
	static Double_t phrq(int i)			{ Double_t phrq[] = { 100, -200, 200, 100, -200, 200 };		return phrq[i]; }		//Out-of-plane angle between recoil system and q [Deg]
	static Double_t xangle(int i)		{ Double_t xangle[] = { 100, 30, 120, 100, 60, 90 };		return xangle[i]; }		//Angle of detected particle with scattered electron [Deg]
	static Double_t Tx_cm(int i)		{ Double_t Tx_cm[] = { 100, -1, 3, 100, -1, 3 };			return Tx_cm[i]; }		//Kinetic Energy (Detected Hadron, in CM frame) [GeV]
	static Double_t Tr_cm(int i)		{ Double_t Tr_cm[] = { 100, -1, 8, 100, -1, 8 };			return Tr_cm[i]; }		//Kinetic Energy (Recoil System, in CM frame) [GeV]
	static Double_t thxq_cm(int i)		{ Double_t thxq_cm[] = { 100, -1, 5, 100, -1, 5 };			return thxq_cm[i]; }	//In-plane angle between detected particle and q in CM frame [Deg]
	static Double_t thrq_cm(int i)		{ Double_t thrq_cm[] = { 100, -1, 5, 100, -1, 5 };			return thrq_cm[i]; }	//In-plane angle between the recoil system and q in CM frame [Deg]
	static Double_t phxq_cm(int i)		{ Double_t phxq_cm[] = { 100, -5, 5, 100, -5, 5 };			return phxq_cm[i]; }	//Out-of-plane angle between detected particle and q in CM frame [Deg]
	static Double_t phrq_cm(int i)		{ Double_t phrq_cm[] = { 100, -5, 5, 100, -5, 5 };			return phrq_cm[i]; }	//Out-of-plane angle between recoil system and q in CM frame [Deg]
	static Double_t Ttot_cm(int i)		{ Double_t Ttot_cm[] = { 100, -1, 8, 100, -1, 8 };			return Ttot_cm[i]; }	//total kinetic energy in CM frame [GeV]
	static Double_t MandelS(int i)		{ Double_t MandelS[] = { 100, 0, 10, 100, 0, 10 };			return MandelS[i]; }	//Mandelstam s for secondary vertex
	static Double_t MandelT(int i)		{ Double_t MandelT[] = { 100, 0, 10, 100, 0, 10 };			return MandelT[i]; }	//Mandelstam t for secondary vertex
	static Double_t MandelU(int i)		{ Double_t MandelU[] = { 100, 0, 10, 100, 0, 10 };			return MandelU[i]; }	//Mandelstam u for secondary vertex
	//----- -----
	static Double_t MM2(int i)			{ Double_t MM2[] = { 100, -1.5, 1.5, 100, -1.5, 1.5 };		return MM2[i]; }		//Missing Mass Squared [(GeV)^2]
	static Double_t thx(int i)			{ Double_t thx[] = {100, 40, 55, 100, 60, 70};				return thx[i]; }		// Hadron arm central angle [Deg]
	static Double_t Pf(int i)			{ Double_t Pf[] = { 100, 1.5, 2.2, 100, 1, 1.8 };			return Pf[i]; }			//final hadron arm momentum [GeV/c]
	//----- -----
	static Double_t Erecoil(int i)		{ Double_t Erecoil[] = { 100, -1, 1, 100, -1, 1 };			return Erecoil[i];}		//Total Energy of Recoil System [GeV/c]
	static Double_t MMk(int i)			{ Double_t MMk[] = { 100, -1, 1, 100, -1, 1 };				return MMk[i];}			//Missing Mass (kaon) [GeV]
	static Double_t MMpi(int i)			{ Double_t MMpi[] = { 100, -1, 1, 100, -1, 1 };				return MMpi[i];}		//Missing mass (pion) [GeV]
	static Double_t Mrecoil(int i)		{ Double_t Mrecoil[] = { 100, -0.2, 1.2, 100, -0.2, 1.2 };	return Mrecoil[i];}		//Invariant Mass of Recoil System [GeV]
	static Double_t px_cm(int i)		{ Double_t px_cm[] = { 100, -1, 100, 100, -1, 100 };		return px_cm[i];}
	//----- -----
	static Double_t Ep(int i)			{ Double_t Ep[] = { 100, -0.1, 2, 100, -0.1, 2 };			return Ep[i]; }


	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------------Acceptance Histogram Bins---------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//----Electron Arm Focal Plane-----
	static Double_t exfp(int i)			{ Double_t exfp[] = { 100, -15, 40, 100, -15, 40 };			return exfp[i];}		//[cm]
	static Double_t expfp(int i)		{ Double_t expfp[] = { 100, -9, 9, 100, -9, 9 };			return expfp[i];}		//[Rad]
	static Double_t eyfp(int i)			{ Double_t eyfp[] = { 100, -25, 15, 100, -25, 15 };			return eyfp[i];}		//[cm]
	static Double_t eypfp(int i)		{ Double_t eypfp[] = { 100, -4, 4, 100, -4, 4 };			return eypfp[i];}		//[Rad]
	//----Electron Arm Reconstructed-----
	static Double_t eytar(int i)		{ Double_t eytar[] = { 100, -6, 6, 100, -6, 6 };			return eytar[i];}		//[cm]
	static Double_t eyptar(int i)		{ Double_t eyptar[] = { 100, -9, 9, 100, -9, 9 };			return eyptar[i];}		//[Rad]
	static Double_t exptar(int i)		{ Double_t exptar[] = { 100, -9, 9, 100, -9, 9 };			return exptar[i];}		//[Rad]
	static Double_t edelta(int i)		{ Double_t edelta[] = { 100, -30, 30, 100, -30, 30 };		return edelta[i];}		//[percent]
	//----Hadron Arm Focal Plane-----
	static Double_t hxfp(int i)			{ Double_t hxfp[] = { 100, -50, 50, 100, -50, 50 };			return hxfp[i];}		//[cm]
	static Double_t hxpfp(int i)		{ Double_t hxpfp[] = { 100, -12, 12, 100, -12, 12 };		return hxpfp[i];}		//[Rad]
	static Double_t hyfp(int i)			{ Double_t hyfp[] = { 100, -35, 35, 100, -35, 35 };			return hyfp[i];}		//[cm]
	static Double_t hypfp(int i)		{ Double_t hypfp[] = { 100, -12, 12, 100, -12, 12 };		return hypfp[i];}		//[Rad]
	//----Hadron Arm Reconstructed-----
	static Double_t hytar(int i)		{ Double_t hytar[] = { 100, -6, 6, 100, -6, 6 };			return hytar[i];}		//[cm]
	static Double_t hyptar(int i)		{ Double_t hyptar[] = { 100, -12, 12, 100, -12, 12 };		return hyptar[i];}		//[Rad]
	static Double_t hxptar(int i)		{ Double_t hxptar[] = { 100, -12, 12, 100, -12, 12 };		return hxptar[i];}		//[Rad]
	static Double_t hdelta(int i)		{ Double_t hdelta[] = { 100, -30, 30, 100, -30, 30 };		return hdelta[i];}		//[percent]
	//----Target Quantities----
	static Double_t tarx(int i)			{ Double_t tarx[] = { 100, -0.3, 0.3, 100, -0.3, 0.3 };		return tarx[i];}		//[cm]
	static Double_t tary(int i)			{ Double_t tary[] = { 100, -0.3, 0.3, 100, -0.3, 0.3 };		return tary[i];}		//[cm]
	static Double_t tarz(int i)			{ Double_t tarz[] = { 100, -15, 15, 100, -15, 15 };			return tarz[i];}		//[cm]
	static Double_t ztar_diff(int i)	{ Double_t ztar_diff[] = { 100, -10, 10, 100, -10, 10 };	return ztar_diff[i];}	//[percent]
	//----Collimator Quantities----
	static Double_t hXColl(int i)		{ Double_t hXColl[] = { 100, -15, 15, 100, -15, 15 };		return hXColl[i];}		//[cm]
	static Double_t hYColl(int i)		{ Double_t hYColl[] = { 100, -15, 15, 100, -15, 15 };		return hYColl[i];}		//[cm]
	static Double_t eXColl(int i)		{ Double_t eXColl[] = { 100, -15, 15, 100, -15, 15 };		return eXColl[i];}		//[cm]
	static Double_t eYColl(int i)		{ Double_t eYColl[] = { 100, -15, 15, 100, -15, 15 };		return eYColl[i];}		//[cm]

	private:
};

#endif