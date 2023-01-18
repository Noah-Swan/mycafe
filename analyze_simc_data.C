#include "utils/parse_utils.h"
#include "utils/hist_utils.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TImageDump.h"

#include "Cuts.h"
#include "Histo.h"
using namespace std;



//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//-----------------------------------------------------------------------------Constants----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
Double_t pi = 3.141592654;
Double_t dtr = pi/180.;
Double_t MP = 0.938272; //Proton Mass GeV
Double_t MD = 1.87561; //GeV
Double_t MN = 0.939566; //Neutron Mass GeV
Double_t me = 0.000510998; //GeV


//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//-----------------------------------------------------------------------------Functions----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
void plotstuff(TH1F*, const char*, const char*, TFile*);
void be_and_af(TH1F*, TH1F*, TH1F*, double, double, double, double, double, double, double, double, const char*, const char*, bool, TFile*);
void be_and_af2d(TH2F*, TH2F*, TH2F*, const char*, const char*, TFile*);


//*****************************************************************************************************************************************************************************
//Begin Main
//*****************************************************************************************************************************************************************************
void analyze_simc_data() 
{
	gStyle->SetOptStat(1);//have stat window
	gStyle->SetOptTitle(1);//have title
	gStyle->SetLabelSize(0.05, "xyz");
	gStyle->SetTitleSize(0.05, "xyz");
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadTopMargin(0.08);
	gStyle->SetPadLeftMargin(0.17);
	gStyle->SetPadRightMargin(0.2);
	gStyle->SetImageScaling(3.);




	//--------------------------------------	
	//Get Initial Information
	//--------------------------------------
	cout << "Enter Target: D2 or C12." << endl;
	string user_target;
	cin >> user_target;
	const char* target = user_target.c_str();

	cout << "Enter run type: SRC or MF." << endl;
	string user_runtype;
	cin >> user_runtype;
	const char* runtype = user_runtype.c_str();

	cout << "Enter rad type: norad or rad." << endl;
	string user_radtype;
	cin >> user_radtype;
	const char* radtype = user_radtype.c_str();

	int binning = -1;
	if (user_runtype == "SRC") { binning = 3; }
	else if (user_runtype == "MF") { binning = 0; }
	else { cout << "Invalid run type." << endl; return; }

	TString InputFileName;
	TFile* inROOT;
	TString OutputFileName;
	TFile* outROOT;
	TString infile;

	if (user_runtype == "MF")
	{
		if (user_target == "C12")
		{
			if (user_radtype == "norad")
			{
				inROOT = new TFile("worksim/cafe_c12_MF_norad.root", "READ");
				OutputFileName = "outHist_C12_MF_norad.root";
				infile = "infiles/cafe_c12_MF_norad.data";
			}
			else if (user_radtype == "rad")
			{
				inROOT = new TFile("worksim/cafe_c12_MF_rad.root", "READ");
				OutputFileName = "outHist_C12_MF_rad.root";
				infile = "infiles/cafe_c12_MF_rad.data";
			}
			else { cout << "Couldn't find file." << endl; }
		}
		else { cout << "Couldn't find file." << endl; }
	}
	else if (user_runtype == "SRC")
	{
		if (user_target == "D2")
		{
			if (user_radtype == "norad")
			{
				inROOT = new TFile("worksim/cafe_d2_SRC_norad.root", "READ");
				OutputFileName = "outHist_d2_MF_norad.root";
				infile = "infiles/cafe_d2_SRC_norad.data";
			}
			else if (user_radtype == "rad")
			{
				inROOT = new TFile("worksim/cafe_d2_SRC_rad.root", "READ");
				OutputFileName = "outHist_d2_MF_rad.root";
				infile = "infiles/cafe_d2_SRC_rad.data";
			}
			else { cout << "Couldn't find file." << endl; }
		}
		else { cout << "Couldn't find file." << endl; }
	}
	else { cout << "Couldn't find file." << endl; }







	int max_loops = 1000000; //maximum number of events to process
	int numLoops = 0;
	
	



	//Get the data tree
	TTree *inputtree;
	inputtree = (TTree*)inROOT->Get("SNT");

	Long64_t nentries;
	nentries = inputtree->GetEntries();
	cout << "There are " << nentries << " events." << endl;    

	//e- arm angle (deg)
	Double_t e_angle = (stod(split(split(FindString("spec%e%theta", infile.Data())[0], '!')[0], '=')[1]));
	//p arm angle (deg)
	Double_t h_angle = (stod(split(split(FindString("spec%p%theta", infile.Data())[0], '!')[0], '=')[1]));

	//SIMC Specific TTree Variable Names
	Double_t Normfac;		inputtree->SetBranchAddress("Normfac", &Normfac);               //normalization factor, defined as : normfac = luminosity * accepted / generated, luminosity = EXPER charge / targetfac
	Double_t Weight;		inputtree->SetBranchAddress("Weight",  &Weight);                //This Weight has the cross section in it
	Double_t Jacobian_corr;	inputtree->SetBranchAddress("Jacobian_corr", &Jacobian_corr);
	Double_t prob_abs;		inputtree->SetBranchAddress("probabs", &prob_abs);				// Probability of absorption of particle in the HMS Collimator
	
	//SIMC Collimator
	Double_t htarx_corr;
	Double_t etarx_corr;
  
	

	//--------------------------------------------------------
	//---------DECLARE HISTOGRAMS-----------------------------
	//--------------------------------------------------------
	//Create TLists to store categorical histograms
	TList* HList = new TList();

	TList* HList2 = new TList();



	//******************************************************************************************************************************************************************
	//-------------------------------------------------Primary Kinematics (electron kinematics) (USED BY DATA AND SIMC)-------------------------------------------------
	//******************************************************************************************************************************************************************
	
	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	Double_t th_e;				TH1F* H1_th_e =		Histo::mk_the_sim(inputtree, th_e, "H1_th_e", HList, binning);		//Electron scattering angle
	Double_t W;					TH1F* H1_W =		Histo::mk_W_sim(inputtree, W, "H1_W", HList, binning);				//Invariant mass
	Double_t Q2;				TH1F* H1_Q2 =		Histo::mk_Q2_sim(inputtree, Q2, "H1_Q2", HList, binning);			//Four-momentum trasfer
	Double_t nu;				TH1F* H1_nu =		Histo::mk_nu_sim(inputtree, nu, "H1_nu", HList, binning);			//Energy Transfer
	Double_t q;					TH1F* H1_q =		Histo::mk_q_sim(inputtree, q, "H1_q", HList, binning);				//Magnitude of the 3-vector q
	//ph_q //Out of plane angle between beamline and q
	Double_t Kf;				TH1F* H1_kf =		Histo::mk_kf_sim(inputtree, Kf, "H1_kf", HList, binning);

	TH1F* H1_th_e_acc =			Histo::mk_the_sim(inputtree, th_e, "H1_th_e_acc", HList, binning);						//Electron scattering angle
	TH1F* H1_W_acc =			Histo::mk_W_sim(inputtree, W, "H1_W_acc", HList, binning);								//Invariant mass
	TH1F* H1_Q2_acc =			Histo::mk_Q2_sim(inputtree, Q2, "H1_Q2_acc", HList, binning);							//Four-momentum trasfer
	TH1F* H1_nu_acc =			Histo::mk_nu_sim(inputtree, nu, "H1_nu_acc", HList, binning);							//Energy Transfer
	TH1F* H1_q_acc =			Histo::mk_q_sim(inputtree, q, "H1_q_acc", HList, binning);								//Magnitude of the 3-vector q
	TH1F* H1_kf_acc =			Histo::mk_kf_sim(inputtree, Kf, "H1_kf_acc", HList, binning);

	TH1F* H1_th_e_acc_kin =		Histo::mk_the_sim(inputtree, th_e, "H1_th_e_acc_kin", HList, binning);					//Electron scattering angle
	TH1F* H1_W_acc_kin =		Histo::mk_W_sim(inputtree, W, "H1_W_acc_kin", HList, binning);							//Invariant mass
	TH1F* H1_Q2_acc_kin =		Histo::mk_Q2_sim(inputtree, Q2, "H1_Q2_acc_kin", HList, binning);						//Four-momentum trasfer
	TH1F* H1_nu_acc_kin =		Histo::mk_nu_sim(inputtree, nu, "H1_nu_acc_kin", HList, binning);						//Energy Transfer
	TH1F* H1_q_acc_kin =		Histo::mk_q_sim(inputtree, q, "H1_q_acc_kin", HList, binning);							//Magnitude of the 3-vector q
	TH1F* H1_kf_acc_kin =		Histo::mk_kf_sim(inputtree, Kf, "H1_kf_acc_kin", HList, binning);

	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	Double_t x_bj;				TH1F* H1_xbj = Histo::mk_xbj_sim(inputtree, x_bj, "H1_xbj", HList, binning);			//B-jorken X  scaling variable
	Double_t th_q;				TH1F* H1_th_q = Histo::mk_thq_sim(inputtree, th_q, "H1_th_q", HList, binning);			//Angle between q and +z (hall coord. system)
	Double_t Ki;				TH1F* H1_ki = Histo::mk_ki_sim(inputtree, Ki, "H1_ki", HList, binning);

	TH1F* H1_xbj_acc =			Histo::mk_xbj_sim(inputtree, x_bj, "H1_xbj_acc", HList, binning);						//B-jorken X  scaling variable
	TH1F* H1_th_q_acc =			Histo::mk_thq_sim(inputtree, th_q, "H1_th_q_acc", HList, binning);						//Angle between q and +z (hall coord. system)
	TH1F* H1_ki_acc =			Histo::mk_ki_sim(inputtree, Ki, "H1_ki_acc", HList, binning);

	TH1F* H1_xbj_acc_kin =		Histo::mk_xbj_sim(inputtree, x_bj, "H1_xbj_acc_kin", HList, binning);					//B-jorken X  scaling variable
	TH1F* H1_th_q_acc_kin =		Histo::mk_thq_sim(inputtree, th_q, "H1_th_q_acc_kin", HList, binning);					//Angle between q and +z (hall coord. system)
	TH1F* H1_ki_acc_kin =		Histo::mk_ki_sim(inputtree, Ki, "H1_ki_acc_kin", HList, binning);



	//******************************************************************************************************************************************************************
	//--------------------------------------------Secondary (Hadron) Kinematics (recoil and missing are used interchageably)--------------------------------------------
	//******************************************************************************************************************************************************************

	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	//Standard Missing Energy for H(e,e'p)
	Double_t Em;			TH1F* H1_emiss = Histo::mk_Em_sim(inputtree, Em, "H1_emiss", HList, binning);
	Double_t pmiss;			TH1F* H1_pmiss = Histo::mk_Pm_sim(inputtree, pmiss, "H1_pmiss", HList, binning);
	Double_t th_pq;			TH1F* H1_th_pq = Histo::mk_thpq_sim(inputtree, th_pq, "H1_th_pq", HList, binning);
	Double_t th_rq;			TH1F* H1_th_rq = Histo::mk_thrq_sim(inputtree, th_rq, "H1_th_rq", HList, binning);
	//cth_rq
	Double_t ph_pq;			TH1F* H1_ph_pq = Histo::mk_phpq_sim(inputtree, ph_pq, "H1_ph_pq", HList, binning);
	Double_t Pf;			TH1F* H1_Pf = Histo::mk_Pf_sim(inputtree, Pf, "H1_Pf", HList, binning);
	Double_t th_p;			TH1F* H1_th_p = Histo::mk_thp_sim(inputtree, th_p, "H1_th_p", HList, binning);

	TH1F* H1_emiss_acc = Histo::mk_Em_sim(inputtree, Em, "H1_emiss_acc", HList, binning);
	TH1F* H1_pmiss_acc = Histo::mk_Pm_sim(inputtree, pmiss, "H1_pmiss_acc", HList, binning);
	TH1F* H1_th_pq_acc = Histo::mk_thpq_sim(inputtree, th_pq, "H1_th_pq_acc", HList, binning);
	TH1F* H1_th_rq_acc = Histo::mk_thrq_sim(inputtree, th_rq, "H1_th_rq_acc", HList, binning);
	TH1F* H1_ph_pq_acc = Histo::mk_phpq_sim(inputtree, ph_pq, "H1_ph_pq_acc", HList, binning);
	TH1F* H1_Pf_acc = Histo::mk_Pf_sim(inputtree, Pf, "H1_Pf_acc", HList, binning);
	TH1F* H1_th_p_acc = Histo::mk_thp_sim(inputtree, th_p, "H1_th_p_acc", HList, binning);

	TH1F* H1_emiss_acc_kin = Histo::mk_Em_sim(inputtree, Em, "H1_emiss_acc_kin", HList, binning);
	TH1F* H1_pmiss_acc_kin = Histo::mk_Pm_sim(inputtree, pmiss, "H1_pmiss_acc_kin", HList, binning);
	TH1F* H1_th_pq_acc_kin = Histo::mk_thpq_sim(inputtree, th_pq, "H1_th_pq_acc_kin", HList, binning);
	TH1F* H1_th_rq_acc_kin = Histo::mk_thrq_sim(inputtree, th_rq, "H1_th_rq_acc_kin", HList, binning);
	TH1F* H1_ph_pq_acc_kin = Histo::mk_phpq_sim(inputtree, ph_pq, "H1_ph_pq_acc_kin", HList, binning);
	TH1F* H1_Pf_acc_kin = Histo::mk_Pf_sim(inputtree, Pf, "H1_Pf_acc_kin", HList, binning);
	TH1F* H1_th_p_acc_kin = Histo::mk_thp_sim(inputtree, th_p, "H1_th_p_acc_kin", HList, binning);

	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	Double_t MM;			TH1F* H1_mmiss = Histo::mk_MM_sim(inputtree, MM, "H1_mmiss", HList, binning);
	Double_t MM2;			TH1F* H1_mmiss2 = Histo::mk_MM2_sim(inputtree, MM2, "H1_mmiss2", HList, binning);
	Double_t Ep;			TH1F* H1_Ep = Histo::mk_Ep_sim(inputtree, Ep, "H1_Ep", HList, binning);

	TH1F* H1_mmiss_acc = Histo::mk_MM_sim(inputtree, MM, "H1_mmiss_acc", HList, binning);
	TH1F* H1_mmiss2_acc = Histo::mk_MM2_sim(inputtree, MM2, "H1_mmiss2_acc", HList, binning);
	TH1F* H1_Ep_acc = Histo::mk_Ep_sim(inputtree, Ep, "H1_Ep_acc", HList, binning);	//Double_t Ep;            //final proton energy (needs to be calculated)

	TH1F* H1_mmiss_acc_kin = Histo::mk_MM_sim(inputtree, MM, "H1_mmiss_acc_kin", HList, binning);
	TH1F* H1_mmiss2_acc_kin = Histo::mk_MM2_sim(inputtree, MM2, "H1_mmiss2_acc_kin", HList, binning);
	TH1F* H1_Ep_acc_kin = Histo::mk_Ep_sim(inputtree, Ep, "H1_Ep_acc_kin", HList, binning);







	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------------Acceptance Histogram Bins---------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	
	//----------------------------------------
	//--------------- HMS Leafs --------------
	//----------------------------------------
	//----- Hadron Arm Focal Plane -----
	Double_t hxfp;			TH1F* H1_hxfp = Histo::mk_hxfp_sim(inputtree, hxfp, "H1_hxfp", HList, binning);
	Double_t hxpfp;			TH1F* H1_hxpfp = Histo::mk_hxpfp_sim(inputtree, hxpfp, "H1_hxpfp", HList, binning);
	Double_t hyfp;			TH1F* H1_hyfp = Histo::mk_hyfp_sim(inputtree, hyfp, "H1_hyfp", HList, binning);
	Double_t hypfp;			TH1F* H1_hypfp = Histo::mk_hypfp_sim(inputtree, hypfp, "H1_hypfp", HList, binning);

	TH1F* H1_hxfp_acc = Histo::mk_hxfp_sim(inputtree, hxfp, "H1_hxfp_acc", HList, binning);
	TH1F* H1_hxpfp_acc = Histo::mk_hxpfp_sim(inputtree, hxpfp, "H1_hxpfp_acc", HList, binning);
	TH1F* H1_hyfp_acc = Histo::mk_hyfp_sim(inputtree, hyfp, "H1_hyfp_acc", HList, binning);
	TH1F* H1_hypfp_acc = Histo::mk_hypfp_sim(inputtree, hypfp, "H1_hypfp_acc", HList, binning);

	TH1F* H1_hxfp_acc_kin = Histo::mk_hxfp_sim(inputtree, hxfp, "H1_hxfp_acc_kin", HList, binning);
	TH1F* H1_hxpfp_acc_kin = Histo::mk_hxpfp_sim(inputtree, hxpfp, "H1_hxpfp_acc_kin", HList, binning);
	TH1F* H1_hyfp_acc_kin = Histo::mk_hyfp_sim(inputtree, hyfp, "H1_hyfp_acc_kin", HList, binning);
	TH1F* H1_hypfp_acc_kin = Histo::mk_hypfp_sim(inputtree, hypfp, "H1_hypfp_acc_kin", HList, binning);

	//----- Hadron Arm Reconstructed Quantities -----
	Double_t hytar;			TH1F* H1_hytar = Histo::mk_hytar_sim(inputtree, hytar, "H1_hytar", HList, binning);
	Double_t hyptar;			TH1F* H1_hyptar = Histo::mk_hyptar_sim(inputtree, hyptar, "H1_hyptar", HList, binning);
	Double_t hxptar;			TH1F* H1_hxptar = Histo::mk_hxptar_sim(inputtree, hxptar, "H1_hxptar", HList, binning);
	Double_t hdelta;			TH1F* H1_hdelta = Histo::mk_hdelta_sim(inputtree, hdelta, "H1_hdelta", HList, binning);

	TH1F* H1_hytar_acc = Histo::mk_hytar_sim(inputtree, hytar, "H1_hytar_acc", HList, binning);
	TH1F* H1_hyptar_acc = Histo::mk_hyptar_sim(inputtree, hyptar, "H1_hyptar_acc", HList, binning);
	TH1F* H1_hxptar_acc = Histo::mk_hxptar_sim(inputtree, hxptar, "H1_hxptar_acc", HList, binning);
	TH1F* H1_hdelta_acc = Histo::mk_hdelta_sim(inputtree, hdelta, "H1_hdelta_acc", HList, binning);

	TH1F* H1_hytar_acc_kin = Histo::mk_hytar_sim(inputtree, hytar, "H1_hytar_acc_kin", HList, binning);
	TH1F* H1_hyptar_acc_kin = Histo::mk_hyptar_sim(inputtree, hyptar, "H1_hyptar_acc_kin", HList, binning);
	TH1F* H1_hxptar_acc_kin = Histo::mk_hxptar_sim(inputtree, hxptar, "H1_hxptar_acc_kin", HList, binning);
	TH1F* H1_hdelta_acc_kin = Histo::mk_hdelta_sim(inputtree, hdelta, "H1_hdelta_acc_kin", HList, binning);

	//----- Target Reconstruction (Hall Coord. System) -----
	Double_t htarx;			TH1F* H1_htarx = Histo::mk_htarx_sim(inputtree, htarx, "H1_htarx", HList, binning);
	Double_t htary;			TH1F* H1_htary = Histo::mk_htary_sim(inputtree, htary, "H1_htary", HList, binning);
	Double_t htarz;			TH1F* H1_htarz = Histo::mk_htarz_sim(inputtree, htarz, "H1_htarz", HList, binning);

	TH1F* H1_htarx_acc = Histo::mk_htarx_sim(inputtree, htarx, "H1_htarx_acc", HList, binning);
	TH1F* H1_htary_acc = Histo::mk_htary_sim(inputtree, htary, "H1_htary_acc", HList, binning);
	TH1F* H1_htarz_acc = Histo::mk_htarz_sim(inputtree, htarz, "H1_htarz_acc", HList, binning);

	TH1F* H1_htarx_acc_kin = Histo::mk_htarx_sim(inputtree, htarx, "H1_htarx_acc_kin", HList, binning);
	TH1F* H1_htary_acc_kin = Histo::mk_htary_sim(inputtree, htary, "H1_htary_acc_kin", HList, binning);
	TH1F* H1_htarz_acc_kin = Histo::mk_htarz_sim(inputtree, htarz, "H1_htarz_acc_kin", HList, binning);

	//----------------------------------------
	//----------- HMS Calculations -----------
	//----------------------------------------
	//----- HMS Collimator -----
	Double_t hXColl;			TH1F* H1_hXColl = Histo::mk_hXColl_sim(inputtree, hXColl, "H1_hXColl", HList, binning);
	TH1F* H1_hXColl_acc = Histo::mk_hXColl_sim(inputtree, hXColl, "H1_hXColl_acc", HList, binning);
	TH1F* H1_hXColl_acc_kin = Histo::mk_hXColl_sim(inputtree, hXColl, "H1_hXColl_acc_kin", HList, binning);

	Double_t hYColl;			TH1F* H1_hYColl = Histo::mk_hYColl_sim(inputtree, hYColl, "H1_hYColl", HList, binning);
	TH1F* H1_hYColl_acc = Histo::mk_hYColl_sim(inputtree, hYColl, "H1_hYColl_acc", HList, binning);
	TH1F* H1_hYColl_acc_kin = Histo::mk_hYColl_sim(inputtree, hYColl, "H1_hYColl_acc_kin", HList, binning);


	//----------------------------------------
	//-------------- SHMS Leafs --------------
	//----------------------------------------
	//----- Electron Arm Focal Plane -----
	Double_t exfp;			TH1F* H1_exfp = Histo::mk_exfp_sim(inputtree, exfp, "H1_exfp", HList, binning);
	Double_t expfp;			TH1F* H1_expfp = Histo::mk_expfp_sim(inputtree, expfp, "H1_expfp", HList, binning);
	Double_t eyfp;			TH1F* H1_eyfp = Histo::mk_eyfp_sim(inputtree, eyfp, "H1_eyfp", HList, binning);
	Double_t eypfp;			TH1F* H1_eypfp = Histo::mk_eypfp_sim(inputtree, eypfp, "H1_eypfp", HList, binning);

	TH1F* H1_exfp_acc = Histo::mk_exfp_sim(inputtree, exfp, "H1_exfp_acc", HList, binning);
	TH1F* H1_expfp_acc = Histo::mk_expfp_sim(inputtree, expfp, "H1_expfp_acc", HList, binning);
	TH1F* H1_eyfp_acc = Histo::mk_eyfp_sim(inputtree, eyfp, "H1_eyfp_acc", HList, binning);
	TH1F* H1_eypfp_acc = Histo::mk_eypfp_sim(inputtree, eypfp, "H1_eypfp_acc", HList, binning);

	TH1F* H1_exfp_acc_kin = Histo::mk_exfp_sim(inputtree, exfp, "H1_exfp_acc_kin", HList, binning);
	TH1F* H1_expfp_acc_kin = Histo::mk_expfp_sim(inputtree, expfp, "H1_expfp_acc_kin", HList, binning);
	TH1F* H1_eyfp_acc_kin = Histo::mk_eyfp_sim(inputtree, eyfp, "H1_eyfp_acc_kin", HList, binning);
	TH1F* H1_eypfp_acc_kin = Histo::mk_eypfp_sim(inputtree, eypfp, "H1_eypfp_acc_kin", HList, binning);

	//----- Electron Arm Reconstructed Quantities
	Double_t eytar;			TH1F* H1_eytar = Histo::mk_eytar_sim(inputtree, eytar, "H1_eytar", HList, binning);
	Double_t eyptar;			TH1F* H1_eyptar = Histo::mk_eyptar_sim(inputtree, eyptar, "H1_eyptar", HList, binning);
	Double_t exptar;			TH1F* H1_exptar = Histo::mk_exptar_sim(inputtree, exptar, "H1_exptar", HList, binning);
	Double_t edelta;			TH1F* H1_edelta = Histo::mk_edelta_sim(inputtree, edelta, "H1_edelta", HList, binning);

	TH1F* H1_eytar_acc = Histo::mk_eytar_sim(inputtree, eytar, "H1_eytar_acc", HList, binning);
	TH1F* H1_eyptar_acc = Histo::mk_eyptar_sim(inputtree, eyptar, "H1_eyptar_acc", HList, binning);
	TH1F* H1_exptar_acc = Histo::mk_exptar_sim(inputtree, exptar, "H1_exptar_acc", HList, binning);
	TH1F* H1_edelta_acc = Histo::mk_edelta_sim(inputtree, edelta, "H1_edelta_acc", HList, binning);

	TH1F* H1_eytar_acc_kin = Histo::mk_eytar_sim(inputtree, eytar, "H1_eytar_acc_kin", HList, binning);
	TH1F* H1_eyptar_acc_kin = Histo::mk_eyptar_sim(inputtree, eyptar, "H1_eyptar_acc_kin", HList, binning);
	TH1F* H1_exptar_acc_kin = Histo::mk_exptar_sim(inputtree, exptar, "H1_exptar_acc_kin", HList, binning);
	TH1F* H1_edelta_acc_kin = Histo::mk_edelta_sim(inputtree, edelta, "H1_edelta_acc_kin", HList, binning);

	//----- Target Reconstruction (Hall Coord. System) -----
	Double_t etary;			TH1F* H1_etary = Histo::mk_etary_sim(inputtree, etary, "H1_etary", HList, binning);
	Double_t etarz;			TH1F* H1_etarz = Histo::mk_etarz_sim(inputtree, etarz, "H1_etarz", HList, binning);

	TH1F* H1_etary_acc = Histo::mk_etary_sim(inputtree, etary, "H1_etary_acc", HList, binning);
	TH1F* H1_etarz_acc = Histo::mk_etarz_sim(inputtree, etarz, "H1_etarz_acc", HList, binning);

	TH1F* H1_etary_acc_kin = Histo::mk_etary_sim(inputtree, etary, "H1_etary_acc_kin", HList, binning);
	TH1F* H1_etarz_acc_kin = Histo::mk_etarz_sim(inputtree, etarz, "H1_etarz_acc_kin", HList, binning);


	//----------------------------------------
	//---------- SHMS Calculations -----------
	//----------------------------------------
	//----- SHMS Collimator -----
	Double_t eXColl;			TH1F* H1_eXColl = Histo::mk_eXColl_sim(inputtree, eXColl, "H1_eXColl", HList, binning);
	Double_t eYColl;			TH1F* H1_eYColl = Histo::mk_eYColl_sim(inputtree, eYColl, "H1_eYColl", HList, binning);

	TH1F* H1_eXColl_acc = Histo::mk_eXColl_sim(inputtree, eXColl, "H1_eXColl_acc", HList, binning);
	TH1F* H1_eYColl_acc = Histo::mk_eYColl_sim(inputtree, eYColl, "H1_eYColl_acc", HList, binning);

	TH1F* H1_eXColl_acc_kin = Histo::mk_eXColl_sim(inputtree, eXColl, "H1_eXColl_acc_kin", HList, binning);
	TH1F* H1_eYColl_acc_kin = Histo::mk_eYColl_sim(inputtree, eYColl, "H1_eYColl_acc_kin", HList, binning);

	//----- Calculated -----
	Double_t ztar_diff;			TH1F* H1_ztar_diff = Histo::mk_ztar_diff_sim(inputtree, ztar_diff, "H1_ztar_diff", HList, binning);
	TH1F* H1_ztar_diff_acc = Histo::mk_ztar_diff_sim(inputtree, ztar_diff, "H1_ztar_diff_acc", HList, binning);
	TH1F* H1_ztar_diff_acc_kin = Histo::mk_ztar_diff_sim(inputtree, ztar_diff, "H1_ztar_diff_acc_kin", HList, binning);



	//--------------------------------------------------------
	//--------------------- 2D HISTOGRAMS --------------------
	//--------------------------------------------------------
	
	//----- Kinematic -----
	TH2F* H2_thrq_pmiss = Histo::mk_thrq_pmiss(inputtree, "H2_thrq_pmiss", HList2, binning);
	//cth_rq
	//w_pmiss
	TH2F* H2_xbj_Q2 = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2", HList2, binning);
	TH2F* H2_Em_Pm = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm", HList2, binning);
	TH2F* H2_the_kf = Histo::mk_the_kf(inputtree, "H2_the_kf", HList2, binning);
	TH2F* H2_thp_Pf = Histo::mk_thp_Pf(inputtree, "H2_thp_Pf", HList2, binning);

	TH2F* H2_thrq_pmiss_acc = Histo::mk_thrq_pmiss(inputtree, "H2_thrq_pmiss_acc", HList2, binning);
	TH2F* H2_xbj_Q2_acc = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_acc", HList2, binning);
	TH2F* H2_Em_Pm_acc = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_acc", HList2, binning);
	TH2F* H2_the_kf_acc = Histo::mk_the_kf(inputtree, "H2_the_kf_acc", HList2, binning);
	TH2F* H2_thp_Pf_acc = Histo::mk_thp_Pf(inputtree, "H2_thp_Pf_acc", HList2, binning);

	TH2F* H2_thrq_pmiss_acc_kin = Histo::mk_thrq_pmiss(inputtree, "H2_thrq_pmiss_acc_kin", HList2, binning);
	TH2F* H2_xbj_Q2_acc_kin = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_acc_kin", HList2, binning);
	TH2F* H2_Em_Pm_acc_kin = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_acc_kin", HList2, binning);
	TH2F* H2_the_kf_acc_kin = Histo::mk_the_kf(inputtree, "H2_the_kf_acc_kin", HList2, binning);
	TH2F* H2_thp_Pf_acc_kin = Histo::mk_thp_Pf(inputtree, "H2_thp_Pf_acc_kin", HList2, binning);
	
	//----- Acceptance -----
	TH2F* H2_hXColl_hYColl = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl", HList2, binning);
	TH2F* H2_eXColl_eYColl = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl", HList2, binning);
	TH2F* H2_hxfp_hyfp = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp", HList2, binning);
	TH2F* H2_exfp_eyfp = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp", HList2, binning);
	TH2F* H2_hxptar_exptar = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar", HList2, binning);
	TH2F* H2_hyptar_eyptar = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar", HList2, binning);
	TH2F* H2_hdelta_edelta = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta", HList2, binning);

	TH2F* H2_hXColl_hYColl_acc = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_acc", HList2, binning);
	TH2F* H2_eXColl_eYColl_acc = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_acc", HList2, binning);
	TH2F* H2_hxfp_hyfp_acc = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_acc", HList2, binning);
	TH2F* H2_exfp_eyfp_acc = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_acc", HList2, binning);
	TH2F* H2_hxptar_exptar_acc = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_acc", HList2, binning);
	TH2F* H2_hyptar_eyptar_acc = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_acc", HList2, binning);
	TH2F* H2_hdelta_edelta_acc = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_acc", HList2, binning);

	TH2F* H2_hXColl_hYColl_acc_kin = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_acc_kin", HList2, binning);
	TH2F* H2_eXColl_eYColl_acc_kin = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_acc_kin", HList2, binning);
	TH2F* H2_hxfp_hyfp_acc_kin = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_acc_kin", HList2, binning);
	TH2F* H2_exfp_eyfp_acc_kin = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_acc_kin", HList2, binning);
	TH2F* H2_hxptar_exptar_acc_kin = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_acc_kin", HList2, binning);
	TH2F* H2_hyptar_eyptar_acc_kin = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_acc_kin", HList2, binning);
	TH2F* H2_hdelta_edelta_acc_kin = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_acc_kin", HList2, binning);



	
	









	//-------------------------------
	//  DEFINE FULL WEIGHT VARIABLES
	//-------------------------------
	// STEP1: Determine the charge factor:
	// definition: total charge deposited on target over a time period
	// SIMC input files are set to 'events / 1mC'

	// Charge factor is the total integrated charge assuming a beam current and run time
	Double_t Ib = 80;		//beam current in (uA) microAmps (micro-Coulombs / sec),   1 mC = 1000 uC
	Double_t time = 1.0;	//estimated time (in hours) a run takes (start - end) of run
	Double_t charge_factor = Ib * time * 3600. / 1000.;  // units in mC

	//target boiling slopes for Hydrofen and Deuterium (during commissioning)
	//Double_t LH2_slope = 0.00063396;	//Double_t LD2_slope = 0.00080029;  //NORM. yield loss / uA

	//more realistic target boiling slopes based on improved measurements later on
	Double_t LH2_slope = 0.0002;
	Double_t LD2_slope = 0.00025;

	// STEP2: Estimate Efficiencies (use efficiencies from commissioning experiment)
	// coin. rates were ~ 2.5 Hz in commissioning,
	//Double_t e_trk      = 0.964;	//Double_t h_trk      = 0.988;	//Double_t daq_lt     = 0.98;   //(it was 0.926 during commissioning due to large logic windows ~100 ns in HMS, but now is smaller)

	//for heep checks
	Double_t e_trk      = 0.99;
	Double_t h_trk      = 0.99;
	Double_t daq_lt     = 0.99;

	Double_t proton_abs = 1.0; // 0.9534;  let assume no proton absorption thru material (since for heep singles, only electron thru SHMS, and does NOT get absorbed) 

	Double_t eff_factor;
	//assuming check singles
	Double_t tgt_boil   = 1. - LD2_slope * Ib;
	eff_factor = 0.5; //e_trk * h_trk * daq_lt * tgt_boil * proton_abs;
	
	Double_t FullWeight;

	//Transparency Factors
	Double_t TF_D2 = 0.9;
	Double_t TF_C12 = 0.6;//0.56
	Double_t TF_Be9 = 0.6;
	Double_t TF_B10 = 0.6;
	Double_t TF_B11 = 0.6;
	Double_t TF_Ca40 = 0.4;
	Double_t TF_Ca48 = 0.4;
	Double_t TF_Fe54 = 0.4;



	//Areal densities (g/cm^2)
	Double_t AD_D2 = 1.67;
	Double_t AD_Be9 = 0.978;
	Double_t AD_B10 = 0.5722;
	Double_t AD_B11 = 0.6344;
	Double_t AD_C12 = 0.5244;
	Double_t AD_Ca40 = 0.8;
	Double_t AD_Ca48 = 0.8;
	Double_t AD_Fe54 = 0.4152;

	//Momentum Distribution (MD) Ratios of MD_A / MD_D2
	Double_t a2_D2 = 1;
	Double_t a2_Be9 = 3.9;
	Double_t a2_B10 = 4;
	Double_t a2_B11 = 4;
	Double_t a2_C12 = 4.5;
	Double_t a2_Ca40 = 4.5;
	Double_t a2_Ca48 = 4.5;
	Double_t a2_Fe54 = 5.2;

	//Definining mean field scale factors
	Double_t MF_C12_scale_fac = TF_C12 * eff_factor * charge_factor;
	Double_t MF_D2_scale_fac = (TF_D2 / TF_C12) * MF_C12_scale_fac * (AD_D2 / AD_C12);
	Double_t MF_Be9_scale_fac = (TF_Be9 / TF_C12) * MF_C12_scale_fac * (AD_Be9 / AD_C12);
	Double_t MF_B10_scale_fac = (TF_B10 / TF_C12) * MF_C12_scale_fac * (AD_B10 / AD_C12);
	Double_t MF_B11_scale_fac = (TF_B11 / TF_C12) * MF_C12_scale_fac * (AD_B11 / AD_C12);
	Double_t MF_Ca40_scale_fac = (TF_Ca40 / TF_C12) * MF_C12_scale_fac * (AD_Ca40 / AD_C12);
	Double_t MF_Ca48_scale_fac = (TF_Ca48 / TF_C12) * MF_C12_scale_fac * (AD_Ca48 / AD_C12);
	Double_t MF_Fe54_scale_fac = (TF_Fe54 / TF_C12) * MF_C12_scale_fac * (AD_Fe54 / AD_C12);

	//Definining SRC scale factors
	Double_t SRC_D2_scale_fac = TF_D2 * eff_factor * charge_factor;
	Double_t SRC_Be9_scale_fac = (TF_Be9 / TF_D2) * SRC_D2_scale_fac * a2_Be9 * (AD_Be9 / AD_D2);
	Double_t SRC_B10_scale_fac = (TF_B10 / TF_D2) * SRC_D2_scale_fac * a2_B10 * (AD_B10 / AD_D2);
	Double_t SRC_B11_scale_fac = (TF_B11 / TF_D2) * SRC_D2_scale_fac * a2_B11 * (AD_B11 / AD_D2);
	Double_t SRC_C12_scale_fac = (TF_C12 / TF_D2) * SRC_D2_scale_fac * a2_C12 * (AD_C12 / AD_D2);
	Double_t SRC_Ca40_scale_fac = (TF_Ca40 / TF_D2) * SRC_D2_scale_fac * a2_Ca40 * (AD_Ca40 / AD_D2);
	Double_t SRC_Ca48_scale_fac = (TF_Ca48 / TF_D2) * SRC_D2_scale_fac * a2_Ca48 * (AD_Ca48 / AD_D2);
	Double_t SRC_Fe54_scale_fac = (TF_Fe54 / TF_D2) * SRC_D2_scale_fac * a2_Fe54 * (AD_Fe54 / AD_D2);



	//**************************************************
	//	Begin Main Loop
	//**************************************************
	for (int i = 0; i < nentries; i++)
	{
		inputtree->GetEntry(i);	//Get the ith entry from the SNT TTree

		//-----Define Additional Kinematic Variables--------
		x_bj = Q2 / (2.*MP*nu);
		th_q = th_rq + th_p;
		ztar_diff =  htarz - etarz;

		Pf = Pf/1000.; //final proton momentum (GeV/c)
		Ep = sqrt(MP*MP + Pf*Pf);

		//Missing Mass
		if(user_runtype == "MF")//generalize to heavier nuclei
		{
			MM2 = Em*Em - pmiss*pmiss;
			MM = sqrt(MM2);
		}
		else//only works for D2
		{
			MM = sqrt( pow(nu+MD-Ep,2) - pmiss*pmiss ); //recoil mass (neutron missing mass)
			MM2 = MM * MM;
		}
		
		//SIMC Collimator (definition based on HCANA collimator)
		htarx_corr = htarx - hxptar*htarz*cos(h_angle*dtr);
		etarx_corr = htarx - exptar*etarz*cos(e_angle*dtr);

		//Define Collimator (same as in HCANA)
		hXColl = htarx_corr + hxptar*168.;	//in cm
		hYColl = hytar + hyptar*168.;
		eXColl = etarx_corr + exptar*253.;
		eYColl = eytar + eyptar*253.-(0.019+40.*.01*0.052)*edelta+(0.00019+40*.01*.00052)*edelta*edelta; //correct for HB horizontal bend	  



		if(user_runtype == "MF")
		{
			FullWeight = (Normfac * Weight) / nentries * MF_C12_scale_fac;
		}
		else if(user_runtype == "SRC")
		{
			FullWeight = (Normfac * Weight) / nentries * SRC_D2_scale_fac;
		}
		
	
		bool accp_hms = false;
		bool accp_shms = false;
		bool accp_z = false;
		bool kinheepsing = false;
		bool kinheepcoin = false;
		bool kin = false;
		bool hmscoll = true;
		bool shmscoll = true;

		
		accp_hms = Cuts::Accp_HMS(hdelta, hxptar / dtr, hyptar / dtr, hmscoll);
		accp_shms = Cuts::Accp_SHMS(edelta, exptar / dtr, eyptar / dtr, shmscoll);
		accp_z = Cuts::Accp_Z(ztar_diff);
		kinheepsing = Cuts::kinHeepSing(Q2, W, x_bj);
		kinheepcoin = Cuts::kinHeepCoin(Q2, W, x_bj, Em, MM);
		if (user_runtype == "SRC")
		{
			kin = Cuts::kinSRC(Q2, pmiss, x_bj, (th_rq / dtr), Em, Em, target);
		}
		else if (user_runtype == "MF")
		{
			kin = Cuts::kinMF(Q2, pmiss, Em, Em, (th_rq / dtr), target);
		}
		else
		{
			cout << "Invalid type." << endl;
		}






		H1_kf->				Fill(Kf/1000., FullWeight);		//Primary (electron) Kinematics
		H1_th_e->			Fill(th_e/dtr, FullWeight);
		H1_Q2->				Fill(Q2, FullWeight);
		H1_xbj->			Fill(x_bj, FullWeight);
		H1_nu->				Fill(nu, FullWeight);
		H1_q->				Fill(q, FullWeight);
		H1_th_q->			Fill(th_q/dtr, FullWeight);
		H1_W->				Fill(W, FullWeight);
		H1_Pf->				Fill(Pf, FullWeight);				//Secondary (Hadron) Kinematics
		H1_th_p->				Fill(th_p/dtr, FullWeight);
		H1_emiss->				Fill(Em, FullWeight);
		H1_pmiss->				Fill(pmiss, FullWeight);     
		H1_th_pq->			Fill(th_rq/dtr, FullWeight);
		H1_th_rq->			Fill(th_rq/dtr, FullWeight);
		H1_mmiss->				Fill(MM, FullWeight);
		H1_mmiss2->			Fill(MM2, FullWeight);
		H1_htarx->			Fill(htarx, FullWeight);			//Target Reconstruction (Hall Coord. System)
		H1_htary->			Fill(htary, FullWeight);
		H1_htarz->			Fill(htarz, FullWeight);
		//H1_etarx->			Fill(htarx, FullWeight);
		H1_etary->			Fill(etary, FullWeight);
		H1_etarz->			Fill(etarz, FullWeight);
		H1_ztar_diff->		Fill(ztar_diff, FullWeight);
		H1_hytar->			Fill(hytar, FullWeight);			//Hadron arm Reconstructed Quantities (xtar, ytar, xptar, yptar, delta)
		H1_hxptar->			Fill(hxptar / dtr, FullWeight);
		H1_hyptar->			Fill(hyptar / dtr, FullWeight);
		H1_hdelta->			Fill(hdelta, FullWeight);
		H1_hxfp->			Fill(hxfp, FullWeight);			//Hadron arm Focal Plane Quantities
		H1_hyfp->			Fill(hyfp, FullWeight);
		H1_hxpfp->			Fill(hxpfp / dtr, FullWeight);
		H1_hypfp->			Fill(hypfp / dtr, FullWeight);
		H1_hXColl->			Fill(hXColl, FullWeight);
		H1_hYColl->			Fill(hYColl, FullWeight);

		H1_eytar->			Fill(eytar, FullWeight);			//Electron Arm Reconstructed Quantities (xtar, ytar, xptar, yptar, delta)
		H1_exptar->			Fill(exptar / dtr, FullWeight);
		H1_eyptar->			Fill(eyptar / dtr, FullWeight);
		H1_edelta->			Fill(edelta, FullWeight);
		H1_exfp->			Fill(exfp, FullWeight);			//Electron Arm Focal Plane Quantities
		H1_eyfp->			Fill(eyfp, FullWeight);
		H1_expfp->			Fill(expfp / dtr, FullWeight);
		H1_eypfp->			Fill(eypfp / dtr, FullWeight);
		H1_eXColl->			Fill(eXColl, FullWeight);
		H1_eYColl->			Fill(eYColl, FullWeight);

		//------- Fill 2D Histos -------
		H2_hxfp_hyfp->				Fill(hyfp, hxfp, FullWeight);			// Xfp vs Yfp
		H2_exfp_eyfp->				Fill(eyfp, exfp, FullWeight);
		H2_hXColl_hYColl->			Fill(hYColl, hXColl, FullWeight);		//2D Collimator Histos
		H2_eXColl_eYColl->			Fill(eYColl, eXColl, FullWeight);
		H2_hxptar_exptar->			Fill(exptar / dtr, hxptar / dtr, FullWeight);		//2D HMS v. SHMS Acceptance Correlations
		H2_hyptar_eyptar->			Fill(eyptar / dtr, hyptar / dtr, FullWeight);
		H2_hdelta_edelta->			Fill(edelta, hdelta, FullWeight);
		H2_thrq_pmiss->				Fill(pmiss, th_rq/dtr, FullWeight);		//This is for the 2D cross section Pm vs. thnq binned in thnq
		H2_thrq_pmiss->				Fill(pmiss, th_rq/dtr, FullWeight);
		H2_the_kf->					Fill(Kf/1000, th_e/dtr, FullWeight);
		H2_thp_Pf->					Fill(Pf, th_p/dtr, FullWeight);
		H2_xbj_Q2->					Fill(Q2, x_bj, FullWeight);
		H2_Em_Pm->					Fill(pmiss, Em, FullWeight);
		
		




		if (accp_hms && accp_shms && accp_z)
		{
			H1_kf_acc->Fill(Kf / 1000., FullWeight);		//Primary (electron) Kinematics
			H1_th_e_acc->Fill(th_e / dtr, FullWeight);
			H1_Q2_acc->Fill(Q2, FullWeight);
			H1_xbj_acc->Fill(x_bj, FullWeight);
			H1_nu_acc->Fill(nu, FullWeight);
			H1_q_acc->Fill(q, FullWeight);
			H1_th_q_acc->Fill(th_q / dtr, FullWeight);
			H1_W_acc->Fill(W, FullWeight);
			H1_Pf_acc->Fill(Pf, FullWeight);				//Secondary (Hadron) Kinematics
			H1_th_p_acc->Fill(th_p / dtr, FullWeight);
			H1_emiss_acc->Fill(Em, FullWeight);
			H1_pmiss_acc->Fill(pmiss, FullWeight);
			H1_th_pq_acc->Fill(th_rq / dtr, FullWeight);
			H1_th_rq_acc->Fill(th_rq / dtr, FullWeight);
			H1_mmiss_acc->Fill(MM, FullWeight);
			H1_mmiss2_acc->Fill(MM2, FullWeight);
			H1_htarx_acc->Fill(htarx, FullWeight);			//Target Reconstruction (Hall Coord. System)
			H1_htary_acc->Fill(htary, FullWeight);
			H1_htarz_acc->Fill(htarz, FullWeight);
			//H1_etarx_acc->Fill(htarx, FullWeight);
			H1_etary_acc->Fill(etary, FullWeight);
			H1_etarz_acc->Fill(etarz, FullWeight);
			H1_ztar_diff_acc->Fill(ztar_diff, FullWeight);
			H1_hytar_acc->Fill(hytar, FullWeight);			//Hadron arm Reconstructed Quantities (xtar, ytar, xptar, yptar, delta)
			H1_hxptar_acc->Fill(hxptar / dtr, FullWeight);
			H1_hyptar_acc->Fill(hyptar / dtr, FullWeight);
			H1_hdelta_acc->Fill(hdelta, FullWeight);
			H1_hxfp_acc->Fill(hxfp, FullWeight);			//Hadron arm Focal Plane Quantities
			H1_hyfp_acc->Fill(hyfp, FullWeight);
			H1_hxpfp_acc->Fill(hxpfp / dtr, FullWeight);
			H1_hypfp_acc->Fill(hypfp / dtr, FullWeight);
			H1_hXColl_acc->Fill(hXColl, FullWeight);
			H1_hYColl_acc->Fill(hYColl, FullWeight);
			
			H1_eytar_acc->Fill(eytar, FullWeight);			//Electron Arm Reconstructed Quantities (xtar, ytar, xptar, yptar, delta)
			H1_exptar_acc->Fill(exptar / dtr, FullWeight);
			H1_eyptar_acc->Fill(eyptar / dtr, FullWeight);
			H1_edelta_acc->Fill(edelta, FullWeight);
			H1_exfp_acc->Fill(exfp, FullWeight);			//Electron Arm Focal Plane Quantities
			H1_eyfp_acc->Fill(eyfp, FullWeight);
			H1_expfp_acc->Fill(expfp / dtr, FullWeight);
			H1_eypfp_acc->Fill(eypfp / dtr, FullWeight);
			H1_eXColl_acc->Fill(eXColl, FullWeight);
			H1_eYColl_acc->Fill(eYColl, FullWeight);

			//------- Fill 2D Histos -------
			H2_hxfp_hyfp_acc->Fill(hyfp, hxfp, FullWeight);			// Xfp vs Yfp
			H2_exfp_eyfp_acc->Fill(eyfp, exfp, FullWeight);
			H2_hXColl_hYColl_acc->Fill(hYColl, hXColl, FullWeight);		//2D Collimator Histos
			H2_eXColl_eYColl_acc->Fill(eYColl, eXColl, FullWeight);
			H2_hxptar_exptar_acc->Fill(exptar / dtr, hxptar / dtr, FullWeight);		//2D HMS v. SHMS Acceptance Correlations
			H2_hyptar_eyptar_acc->Fill(eyptar / dtr, hyptar / dtr, FullWeight);
			H2_hdelta_edelta_acc->Fill(edelta, hdelta, FullWeight);
			H2_thrq_pmiss_acc->Fill(pmiss, th_rq / dtr, FullWeight);		//This is for the 2D cross section Pm vs. thnq binned in thnq
			H2_thrq_pmiss_acc->Fill(pmiss, th_rq / dtr, FullWeight);
			H2_the_kf_acc->Fill(Kf / 1000, th_e / dtr, FullWeight);
			H2_thp_Pf_acc->Fill(Pf, th_p / dtr, FullWeight);
			H2_xbj_Q2_acc->Fill(Q2, x_bj, FullWeight);
			H2_Em_Pm_acc->Fill(pmiss, Em, FullWeight);
		}




		


		if (accp_hms && accp_shms && accp_z && kin)
		{
			H1_kf_acc_kin->Fill(Kf / 1000., FullWeight);		//Primary (electron) Kinematics
			H1_th_e_acc_kin->Fill(th_e / dtr, FullWeight);
			H1_Q2_acc_kin->Fill(Q2, FullWeight);
			H1_xbj_acc_kin->Fill(x_bj, FullWeight);
			H1_nu_acc_kin->Fill(nu, FullWeight);
			H1_q_acc_kin->Fill(q, FullWeight);
			H1_th_q_acc_kin->Fill(th_q / dtr, FullWeight);
			H1_W_acc_kin->Fill(W, FullWeight);
			H1_Pf_acc_kin->Fill(Pf, FullWeight);				//Secondary (Hadron) Kinematics
			H1_th_p_acc_kin->Fill(th_p / dtr, FullWeight);
			H1_emiss_acc_kin->Fill(Em, FullWeight);
			H1_pmiss_acc_kin->Fill(pmiss, FullWeight);
			H1_th_pq_acc_kin->Fill(th_rq / dtr, FullWeight);
			H1_th_rq_acc_kin->Fill(th_rq / dtr, FullWeight);
			H1_mmiss_acc_kin->Fill(MM, FullWeight);
			H1_mmiss2_acc_kin->Fill(MM2, FullWeight);
			H1_htarx_acc_kin->Fill(htarx, FullWeight);			//Target Reconstruction (Hall Coord. System)
			H1_htary_acc_kin->Fill(htary, FullWeight);
			H1_htarz_acc_kin->Fill(htarz, FullWeight);
			//H1_etarx_acc_kin->Fill(htarx, FullWeight);
			H1_etary_acc_kin->Fill(etary, FullWeight);
			H1_etarz_acc_kin->Fill(etarz, FullWeight);
			H1_ztar_diff_acc_kin->Fill(ztar_diff, FullWeight);
			H1_hytar_acc_kin->Fill(hytar, FullWeight);			//Hadron arm Reconstructed Quantities (xtar, ytar, xptar, yptar, delta)
			H1_hxptar_acc_kin->Fill(hxptar / dtr, FullWeight);
			H1_hyptar_acc_kin->Fill(hyptar / dtr, FullWeight);
			H1_hdelta_acc_kin->Fill(hdelta, FullWeight);
			H1_hxfp_acc_kin->Fill(hxfp, FullWeight);			//Hadron arm Focal Plane Quantities
			H1_hyfp_acc_kin->Fill(hyfp, FullWeight);
			H1_hxpfp_acc_kin->Fill(hxpfp / dtr, FullWeight);
			H1_hypfp_acc_kin->Fill(hypfp / dtr, FullWeight);
			H1_hXColl_acc_kin->Fill(hXColl, FullWeight);
			H1_hYColl_acc_kin->Fill(hYColl, FullWeight);
			
			H1_eytar_acc_kin->Fill(eytar, FullWeight);			//Electron Arm Reconstructed Quantities (xtar, ytar, xptar, yptar, delta)
			H1_exptar_acc_kin->Fill(exptar / dtr, FullWeight);
			H1_eyptar_acc_kin->Fill(eyptar / dtr, FullWeight);
			H1_edelta_acc_kin->Fill(edelta, FullWeight);
			H1_exfp_acc_kin->Fill(exfp, FullWeight);			//Electron Arm Focal Plane Quantities
			H1_eyfp_acc_kin->Fill(eyfp, FullWeight);
			H1_expfp_acc_kin->Fill(expfp / dtr, FullWeight);
			H1_eypfp_acc_kin->Fill(eypfp / dtr, FullWeight);
			H1_eXColl_acc_kin->Fill(eXColl, FullWeight);
			H1_eYColl_acc_kin->Fill(eYColl, FullWeight);

			//------- Fill 2D Histos -------
			H2_hxfp_hyfp_acc_kin->Fill(hyfp, hxfp, FullWeight);			// Xfp vs Yfp
			H2_exfp_eyfp_acc_kin->Fill(eyfp, exfp, FullWeight);
			H2_hXColl_hYColl_acc_kin->Fill(hYColl, hXColl, FullWeight);		//2D Collimator Histos
			H2_eXColl_eYColl_acc_kin->Fill(eYColl, eXColl, FullWeight);
			H2_hxptar_exptar_acc_kin->Fill(exptar / dtr, hxptar / dtr, FullWeight);		//2D HMS v. SHMS Acceptance Correlations
			H2_hyptar_eyptar_acc_kin->Fill(eyptar / dtr, hyptar / dtr, FullWeight);
			H2_hdelta_edelta_acc_kin->Fill(edelta, hdelta, FullWeight);
			H2_thrq_pmiss_acc_kin->Fill(pmiss, th_rq / dtr, FullWeight);		//This is for the 2D cross section Pm vs. thnq binned in thnq			
			H2_thrq_pmiss_acc_kin->Fill(pmiss, th_rq / dtr, FullWeight);
			H2_the_kf_acc_kin->Fill(Kf / 1000, th_e / dtr, FullWeight);
			H2_thp_Pf_acc_kin->Fill(Pf, th_p / dtr, FullWeight);
			H2_xbj_Q2_acc_kin->Fill(Q2, x_bj, FullWeight);
			H2_Em_Pm_acc_kin->Fill(pmiss, Em, FullWeight);
		}


		if(numLoops > max_loops){ cout<<"Exceeded Max Loops"<<endl;	break;} //leave the loop after completing the loop max_loops times
		numLoops++;
	}//end for



	Double_t error;
	Double_t Rate;

	if(user_runtype == "MF")
	{
		Rate = H1_pmiss_acc_kin->IntegralAndError(H1_pmiss_acc_kin->FindBin(0), H1_pmiss_acc_kin->FindBin(0.25), error, "");
		cout << "The integral over P_{m} for MF D2 is: " << (Rate * MF_D2_scale_fac / MF_C12_scale_fac) << endl;
		cout << "The integral over P_{m} for MF Be9 is: " << (Rate * MF_Be9_scale_fac / MF_C12_scale_fac) << endl;
		cout << "The integral over P_{m} for MF B10 is: " <<(Rate * MF_B10_scale_fac / MF_C12_scale_fac) << endl;
		cout << "The integral over P_{m} for MF B11 is: " << (Rate * MF_B11_scale_fac / MF_C12_scale_fac) << endl;
		cout << "The integral over P_{m} for MF C12 is: " << Rate << endl;
		cout << "The integral over P_{m} for MF Ca40 is: " << (Rate * MF_Ca40_scale_fac / MF_C12_scale_fac) << endl;
		cout << "The integral over P_{m} for MF Ca48 is: " << (Rate * MF_Ca48_scale_fac / MF_C12_scale_fac) << endl;
		cout << "The integral over P_{m} for MF Fe54 is: " << (Rate * MF_Fe54_scale_fac / MF_C12_scale_fac) << endl;
	}
	else
	{
		Rate = H1_pmiss_acc_kin->IntegralAndError(H1_pmiss_acc_kin->FindBin(0.35), H1_pmiss_acc_kin->FindBin(2), error, "");
		cout << "The integral over P_{m} for SRC D2 is: " << Rate << endl;
		cout << "The integral over P_{m} for SRC Be9 is: " << (Rate * SRC_Be9_scale_fac / SRC_D2_scale_fac) << endl;
		cout << "The integral over P_{m} for SRC B10 is: " << (Rate * SRC_B10_scale_fac / SRC_D2_scale_fac) << endl;
		cout << "The integral over P_{m} for SRC B11 is: " << (Rate * SRC_B11_scale_fac / SRC_D2_scale_fac) << endl;
		cout << "The integral over P_{m} for SRC C12 is: " << (Rate * SRC_C12_scale_fac / SRC_D2_scale_fac) << endl;
		cout << "The integral over P_{m} for SRC Ca40 is: " << (Rate * SRC_Ca40_scale_fac / SRC_D2_scale_fac) << endl;
		cout << "The integral over P_{m} for SRC Ca48 is: " << (Rate * SRC_Ca48_scale_fac / SRC_D2_scale_fac) << endl;
		cout << "The integral over P_{m} for SRC Fe54 is: " << (Rate * SRC_Fe54_scale_fac / SRC_D2_scale_fac) << endl;
	}



	outROOT = new TFile(OutputFileName, "RECREATE");	//Create Output ROOTfile
	
	outROOT->mkdir("plots");								//Make directories to store histograms based on Kinematic
	outROOT->cd("plots");									//Write Kinematics histos to kin_plots directory
	HList->Write();	

	outROOT->mkdir("plots2");								//Make directories to store histograms based on Kinematic
	outROOT->cd("plots2");									//Write Kinematics histos to kin_plots directory
	HList2->Write();








	//--------------------------------------
	//PID Histograms Bins
	//--------------------------------------
	//----- HMS -----
	double y2f = 0;
	double y1f = 0;

	Double_t hdelta_min = Cuts::hdelta(0), hdelta_max = Cuts::hdelta(1);
	Double_t hxptar_min = Cuts::hxptar(0), hxptar_max = Cuts::hxptar(1);
	Double_t hyptar_min = Cuts::hyptar(0), hyptar_max = Cuts::hyptar(1);
	Double_t edelta_min = Cuts::edelta(0), edelta_max = Cuts::edelta(1);
	Double_t exptar_min = Cuts::exptar(0), exptar_max = Cuts::exptar(1);
	Double_t eyptar_min = Cuts::eyptar(0), eyptar_max = Cuts::eyptar(1);
	Double_t ztarDiff_min = Cuts::ztarDiff(0), ztarDiff_max = Cuts::ztarDiff(1);
	//----- Hadron Arm Focal Plane -----
	be_and_af(H1_hxfp, H1_hxfp_acc, H1_hxfp_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hxfp.png", "hxfp", false, outROOT);
	be_and_af(H1_hxpfp, H1_hxpfp_acc, H1_hxpfp_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hxpfp.png", "hxpfp", false, outROOT);
	be_and_af(H1_hyfp, H1_hyfp_acc, H1_hyfp_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hyfp.png", "hyfp", false, outROOT);
	be_and_af(H1_hypfp, H1_hypfp_acc, H1_hypfp_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hypfp.png", "hypfp", false, outROOT);
	//----- Hadron Arm Reconstructed Quantities -----
	be_and_af(H1_hytar, H1_hytar_acc, H1_hytar_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hytar.png", "hytar", false, outROOT);
	be_and_af(H1_hyptar, H1_hyptar_acc, H1_hyptar_acc_kin, hyptar_min, 0, hyptar_min, y1f, hyptar_max, 0, hyptar_max, y2f, "hyptar.png", "hyptar", false, outROOT);
	be_and_af(H1_hxptar, H1_hxptar_acc, H1_hxptar_acc_kin, hxptar_min, 0, hxptar_min, y1f, hxptar_max, 0, hxptar_max, y2f, "hxptar.png", "hxptar", false, outROOT);
	be_and_af(H1_hdelta, H1_hdelta_acc, H1_hdelta_acc_kin, hdelta_min, 0, hdelta_min, y1f, hdelta_max, 0, hdelta_max, y2f, "hdelta.png", "hdelta", true, outROOT);
	//----- Target Reconstruction (Hall Coord. System) -----
	be_and_af(H1_htarx, H1_htarx_acc, H1_htarx_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htarx.png", "htarx", false, outROOT);
	be_and_af(H1_htary, H1_htary_acc, H1_htary_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htary.png", "htary", false, outROOT);
	be_and_af(H1_htarz, H1_htarz_acc, H1_htarz_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htarz.png", "htarz", false, outROOT);
	//----- HMS Collimator -----
	be_and_af(H1_hXColl, H1_hXColl_acc, H1_hXColl_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hXColl.png", "hXColl", false, outROOT);
	be_and_af(H1_hYColl, H1_hYColl_acc, H1_hYColl_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hYColl.png", "hYColl", false, outROOT);

	//----- Hadron Arm Focal Plane -----
	be_and_af(H1_exfp, H1_exfp_acc, H1_exfp_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "exfp.png", "exfp", false, outROOT);
	be_and_af(H1_expfp, H1_expfp_acc, H1_expfp_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "expfp.png", "expfp", false, outROOT);
	be_and_af(H1_eyfp, H1_eyfp_acc, H1_eyfp_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eyfp.png", "eyfp", false, outROOT);
	be_and_af(H1_eypfp, H1_eypfp_acc, H1_eypfp_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eypfp.png", "eypfp", false, outROOT);
	//----- Hadron Arm Reconstructed Quantities -----
	be_and_af(H1_eytar, H1_eytar_acc, H1_eytar_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eytar.png", "eytar", false, outROOT);
	be_and_af(H1_eyptar, H1_eyptar_acc, H1_eyptar_acc_kin, eyptar_min, 0, eyptar_min, y1f, eyptar_max, 0, eyptar_max, y2f, "eyptar.png", "eyptar", false, outROOT);
	be_and_af(H1_exptar, H1_exptar_acc, H1_exptar_acc_kin, exptar_min, 0, exptar_min, y1f, exptar_max, 0, exptar_max, y2f, "exptar.png", "exptar", false, outROOT);
	be_and_af(H1_edelta, H1_edelta_acc, H1_edelta_acc_kin, edelta_min, 0, edelta_min, y1f, edelta_max, 0, edelta_max, y2f, "edelta.png", "edelta", true, outROOT);
	//----- Target Reconstruction (Hall Coord. System) -----
	//be_and_af(H1_etarx, H1_etarx, H1_etarx_acc, H1_etarx_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etarx.png", "etarx", false, outROOT);
	be_and_af(H1_etary, H1_etary_acc, H1_etary_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etary.png", "etary", false, outROOT);
	be_and_af(H1_etarz, H1_etarz_acc, H1_etarz_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etarz.png", "etarz", false, outROOT);
	//----- SHMS Collimator -----
	be_and_af(H1_eXColl, H1_eXColl_acc, H1_eXColl_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eXColl.png", "eXColl", false, outROOT);
	be_and_af(H1_eYColl, H1_eYColl_acc, H1_eYColl_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eYColl.png", "eYColl", false, outROOT);

	//----- -----
	be_and_af(H1_ztar_diff, H1_ztar_diff_acc, H1_ztar_diff_acc_kin, ztarDiff_min, 0, ztarDiff_min, y1f, ztarDiff_max, 0, ztarDiff_max, y2f, "ztar_diff.png", "ztar_diff", true, outROOT);



	//--------------------------------------
	//Primary Electron Kinematics
	//--------------------------------------

	Double_t Q2_min = -1000, Q2_max = -1000;
	Double_t W_min = -1000, W_max = -1000;
	Double_t xbj_min = -1000, xbj_max = -1000;
	Double_t Em_min = -1000, Em_max = -1000;
	Double_t MM_min = -1000, MM_max = -1000;
	Double_t Pm_min = -1000, Pm_max = -1000;
	Double_t Src_thrq_min = -1000, Src_thrq_max = -1000;
	Double_t Mf_thrq_min = -1000, Mf_thrq_max = -1000;

	/*
	if (runtype == Heep)
	{
		Q2_min = Cuts::heep_Q2(0), Q2_max = Cuts::heep_Q2(1);
		W_min = Cuts::heep_W(0), W_max = Cuts::heep_W(1);
		xbj_min = Cuts::heep_xbj(0), xbj_max = Cuts::heep_xbj(1);
		Em_min = Cuts::heep_Em(0), Em_max = Cuts::heep_Em(1);
		MM_min = Cuts::heep_MM(0), MM_max = Cuts::heep_MM(1);
	}
	*/
	if (user_runtype == "SRC")
	{
		Q2_min = Cuts::SRC_Q2(0), Q2_max = Cuts::SRC_Q2(1);
		Pm_min = Cuts::SRC_Pm(0), Pm_max = Cuts::SRC_Pm(1);
		Em_min = Cuts::d2SRC_Em(0), Em_max = Cuts::d2SRC_Em(1);
		Src_thrq_min = Cuts::SRC_thrq(0), Src_thrq_max = Cuts::SRC_thrq(1);
	}
	else if (user_runtype == "MF")
	{
		Q2_min = Cuts::SRC_Q2(0), Q2_max = Cuts::SRC_Q2(1);
		Pm_min = Cuts::MF_Pm(0), Pm_max = Cuts::MF_Pm(1);
		Mf_thrq_min = Cuts::MF_thrq(0), Mf_thrq_max = Cuts::MF_thrq(1);
		if (user_target == "D2")
		{
			Em_min = Cuts::d2MF_Em(0), Em_max = Cuts::d2MF_Em(1);
		}
		else
		{
			Em_min = Cuts::MF_Em(0), Em_max = Cuts::MF_Em(1);
		}
	}


	be_and_af(H1_Q2, H1_Q2_acc, H1_Q2_acc_kin, Q2_min, 0, Q2_min, y1f, Q2_max, 0, Q2_max, y2f, "Q2.png", "Q2", true, outROOT);
	be_and_af(H1_xbj, H1_xbj_acc, H1_xbj_acc_kin, xbj_min, 0, xbj_min, y1f, xbj_max, 0, xbj_max, y2f, "xbj.png", "xbj", true, outROOT);
	be_and_af(H1_W, H1_W_acc, H1_W_acc_kin, W_min, 0, W_min, y1f, W_max, 0, W_max, y2f, "W.png", "W", true, outROOT);
	be_and_af(H1_nu, H1_nu_acc, H1_nu_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "nu.png", "nu", false, outROOT);
	//be_and_af(H1_ph_q, H1_ph_q_acc, H1_ph_q_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_q.png", "ph_q", false, outROOT);
	be_and_af(H1_q, H1_q_acc, H1_q_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "q.png", "q", false, outROOT);
	be_and_af(H1_th_e, H1_th_e_acc, H1_th_e_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "the.png", "the", false, outROOT);
	be_and_af(H1_th_q, H1_th_q_acc, H1_th_q_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "thq.png", "thq", false, outROOT);



	//--------------------------------------
	//Secondary Hadron Kinematics
	//--------------------------------------
	be_and_af(H1_pmiss, H1_pmiss_acc, H1_pmiss_acc_kin, Pm_min, 0, Pm_min, y1f, Pm_max, 0, Pm_max, y2f, "pmiss.png", "pmiss", true, outROOT);
	be_and_af(H1_th_rq, H1_th_rq_acc, H1_th_rq_acc_kin, Src_thrq_min, 0, Src_thrq_min, y1f, Src_thrq_max, 0, Src_thrq_max, y2f, "th_rq.png", "th_rq", true, outROOT);
	be_and_af(H1_emiss, H1_emiss_acc, H1_emiss_acc_kin, Em_min, 0, Em_min, y1f, Em_max, 0, Em_max, y2f, "emiss.png", "emiss", true, outROOT);
	//be_and_af(H1_ph_rq, H1_ph_rq_acc, H1_ph_rq_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_rq.png", "ph_rq", false, outROOT);
	//be_and_af(H1_ph_pq, H1_ph_pq_acc, H1_ph_pq_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_pq.png", "ph_pq", false, outROOT);
	be_and_af(H1_th_pq, H1_th_pq_acc, H1_th_pq_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "th_pq.png", "th_pq", false, outROOT);
	be_and_af(H1_mmiss, H1_mmiss_acc, H1_mmiss_acc_kin, MM_min, 0, MM_min, y1f, MM_max, 0, MM_max, y2f, "mmiss.png", "mmiss", true, outROOT);
	//mm2
	be_and_af(H1_th_p, H1_th_p_acc, H1_th_p_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "th_p.png", "th_p", false, outROOT);

	be_and_af2d(H2_hxfp_hyfp, H2_hxfp_hyfp_acc, H2_hxfp_hyfp_acc_kin, "H2_hxfp_hyfp.png", "H2_hxfp_hyfp", outROOT);
	be_and_af2d(H2_exfp_eyfp, H2_exfp_eyfp_acc, H2_exfp_eyfp_acc_kin, "H2_exfp_eyfp.png", "H2_exfp_eyfp", outROOT);
	be_and_af2d(H2_hxptar_exptar, H2_hxptar_exptar_acc, H2_hxptar_exptar_acc_kin, "H2_hxptar_exptar.png", "H2_hxptar_exptar", outROOT);
	be_and_af2d(H2_hyptar_eyptar, H2_hyptar_eyptar_acc, H2_hyptar_eyptar_acc_kin, "H2_hyptar_eyptar.png", "H2_hyptar_eyptar", outROOT);
	be_and_af2d(H2_hdelta_edelta, H2_hdelta_edelta_acc, H2_hdelta_edelta_acc_kin, "H2_hdelta_edelta.png", "H2_hdelta_edelta", outROOT);
	be_and_af2d(H2_hXColl_hYColl, H2_hXColl_hYColl_acc, H2_hXColl_hYColl_acc_kin, "H2_hXColl_hYColl.png", "H2_hXColl_hYColl", outROOT);
	be_and_af2d( H2_eXColl_eYColl, H2_eXColl_eYColl_acc, H2_eXColl_eYColl_acc_kin, "H2_eXColl_eYColl.png", "H2_eXColl_eYColl", outROOT);
	be_and_af2d(H2_xbj_Q2, H2_xbj_Q2_acc, H2_xbj_Q2_acc_kin, "H2_xbj_Q2.png", "H2_xbj_Q2", outROOT);
	be_and_af2d(H2_the_kf, H2_the_kf_acc, H2_the_kf_acc_kin, "H2_the_kf.png", "H2_the_kf", outROOT);
	be_and_af2d(H2_Em_Pm, H2_Em_Pm_acc, H2_Em_Pm_acc_kin, "H2_Em_Pm.png", "H2_Em_Pm", outROOT);
	be_and_af2d(H2_thrq_pmiss, H2_thrq_pmiss_acc, H2_thrq_pmiss_acc_kin, "H2_thrq_pmiss.png", "H2_thrq_pmiss", outROOT);
	//be_and_af2d(H2_cthrq_pmiss, H2_cthrq_pmiss_acc, H2_cthrq_pmiss_acc_kin, "H2_cthrq_pmiss.png", "H2_cthrq_pmiss", outROOT);
	be_and_af2d(H2_thp_Pf, H2_thp_Pf_acc, H2_thp_Pf_acc_kin, "H2_thp_Pf.png", "H2_thp_Pf", outROOT);
	//be_and_af2d(H2_W_pmiss, H2_W_pmiss_acc, H2_W_pmiss_acc_kin, "H2_W_pmiss.png", "H2_W_pmiss", outROOT);








	outROOT->Close();										//Close File
}
//*****************************************************************************************************************************************************************************
//End Main
//*****************************************************************************************************************************************************************************



void plotstuff(TH1F* hist1, const char* name, const char* title, TFile* outRoot)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	hist1->Draw();

	outRoot->cd(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}



void be_and_af(TH1F* hist1, TH1F* hist2, TH1F* hist3, double x1i, double y1i, double x1f, double y1f, double x2i, double y2i, double x2f, double y2f, const char* name, const char* title, bool line, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->Divide(2, 2);

	c->cd(1); hist1->Draw("HIST");

	if (line)
	{
		y1f = hist1->GetMaximum();
		y2f = hist1->GetMaximum();
		TLine* line1 = new TLine(x1i, y1i, x1f, y1f); line1->SetLineColor(2); line1->SetLineWidth(2); line1->Draw("SAME");
		TLine* line2 = new TLine(x2i, y2i, x2f, y2f); line2->SetLineColor(2); line2->SetLineWidth(2); line2->Draw("SAME");
	}

	c->cd(2); hist2->Draw("HIST SAME");

	if (line)
	{
		y1f = hist2->GetMaximum();
		y2f = hist2->GetMaximum();
		TLine* line3 = new TLine(x1i, y1i, x1f, y1f); line3->SetLineColor(2); line3->SetLineWidth(2); line3->Draw("SAME");
		TLine* line4 = new TLine(x2i, y2i, x2f, y2f); line4->SetLineColor(2); line4->SetLineWidth(2); line4->Draw("SAME");
	}

	c->cd(4); hist3->Draw("HIST SAME");

	if (line)
	{
		y1f = hist3->GetMaximum();
		y2f = hist3->GetMaximum();
		TLine* line4 = new TLine(x1i, y1i, x1f, y1f); line4->SetLineColor(2); line4->SetLineWidth(2); line4->Draw("SAME");
		TLine* line5 = new TLine(x2i, y2i, x2f, y2f); line5->SetLineColor(2); line5->SetLineWidth(2); line5->Draw("SAME");
	}

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}



void be_and_af2d(TH2F* hist1, TH2F* hist2, TH2F* hist3, const char* name, const char* title, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->Divide(2, 2);

	c->cd(1); hist1->Draw("Col");
	c->cd(2); hist2->Draw("Col SAME");
	c->cd(4); hist3->Draw("Col SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}