//simulation uses 1mc
//subtract by weight instead of nuclei?

//calculate ki, Kf, & Pf
//already built cos(thrq)
//centralize header files with symbolic links

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

#include "../../../../header_files/Histo.h"
#include "../../../../header_files/Cuts.h"
using namespace std;


//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//-----------------------------------------------------------------------------Constants----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
Double_t pi = TMath::Pi();			//Pi
Double_t dtr = pi/180.;				//Conversion between Rad & Deg
Double_t MP = 0.938272;				//Proton Mass GeV
Double_t MD = 1.87561;				//GeV
Double_t MN = 0.939566;				//Neutron Mass GeV
Double_t me = 0.000510998;			//GeV


//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//-----------------------------------------------------------------------------Functions----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
void plotstuff(TH1F*, const char*, const char*, TFile*);
void be_and_af(TH1F*, TH1F*, TH1F*, TH1F*, double, double, double, double, double, double, double, double, const char*, const char*, bool, TFile*);
void be_and_af2d(TH2F*, TH2F*, TH2F*, TH2F*, const char*, const char*, TFile*);


//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//----------------------------------------------------------------------------Begin Main----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
void useskims()
{
	//--------------------------------------	
	//Histogram Settings
	//--------------------------------------
	gStyle->SetOptStat(1);				//(1) for a stats window. (0) for no stats window.
	gStyle->SetOptTitle(1);				//(1) for a title. (0) for no title.
	gStyle->SetLabelSize(0.05, "xyz");	//Text size for the axis labels.
	gStyle->SetTitleSize(0.05, "xyz");	//Text size for the title.
	gStyle->SetPadBottomMargin(0.15);	//Margin spacing on the bottom of the histogram.
	gStyle->SetPadTopMargin(0.08);		//Margin spacing on the top of the histogram.
	gStyle->SetPadLeftMargin(0.17);		//Margin spacing on the left of the histogram.
	gStyle->SetPadRightMargin(0.2);		//Margin spacing on the right of the histogram.
	gStyle->SetImageScaling(3.);		//????.


	//--------------------------------------	
	//Get Initial Information
	//--------------------------------------
	cout << "Enter Target: LD2, Be9, B10, B11, C12, Ca40, Ca48, or Fe54." << endl;
	string user_target;
	cin >> user_target;
	const char* target = user_target.c_str();

	cout << "Enter run type: SRC or MF." << endl;
	string user_runtype;
	cin >> user_runtype;
	const char* runtype = user_runtype.c_str();

	int binning = -1;
	if (user_runtype == "SRC") { binning = 3; }
	else if (user_runtype == "MF") { binning = 0; }
	else { cout << "Invalid run type." << endl; return; }






	

	//******************************************************************************************************************************************************************
	//---------------------------------------------------------------------------Manage FILES---------------------------------------------------------------------------
	//******************************************************************************************************************************************************************

	//Just here as a reminder that it is possible
	//TChain * inputtree = new TChain("skim");
	//inputtree->Add("qe_skim_run_11286.root");		
	//inputtree->Add("qe_skim_run_11287.root");	//And so on
	//TFile* inROOT = new TFile("../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/cafe_prod_Fe54_MF_16982_-1_skimmed.root", "READ");
	TFile* inROOT;

	if (user_runtype == "MF")
	{
		if (user_target == "LD2")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/d2_MF_skim_total.root", "READ");
		}
		else if (user_target == "Be9")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/be9_MF_skim_total.root", "READ");
		}
		else if (user_target == "B10")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/b10_MF_skim_total.root", "READ");
		}
		else if (user_target == "B11")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/b11_MF_skim_total.root", "READ");
		}
		else if (user_target == "C12")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/c12_MF_skim_total.root", "READ");
		}
		else if (user_target == "Ca40")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/ca40_MF_skim_total.root", "READ");
		}
		else if (user_target == "Ca48")
		{
			//inROOT = new TFile("../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/ca48_MF_skim_total.root", "READ");//
			//TChain* inROOT = new TChain("skim");
			//inROOT->Add("../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/cafe_prod_Ca48_MF_17093_-1_skimmed.root");//excluding 978 & 979 due to contamination
			//inROOT->Add("../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/cafe_prod_Ca48_MF_17094_-1_skimmed.root");
			//inROOT->Add("../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/cafe_prod_Ca48_MF_17096_-1_skimmed.root");
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/cafe_prod_Ca48_MF_17096_-1_skimmed.root", "READ");
		}
		else if (user_target == "Fe54")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/fe54_MF_skim_total.root", "READ");
		}
		else
		{
			cout << "Couldn't find file." << endl;
		}
	}
	else if (user_runtype == "SRC")
	{
		if (user_target == "LD2")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/cafe_prod_LD2_SRC_17134_-1_skimmed.root", "READ");//
		}
		else if (user_target == "Be9")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/be9_SRC_skim_total.root", "READ");
		}
		else if (user_target == "B10")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/b10_SRC_skim_total.root", "READ");
		}
		else if (user_target == "B11")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/b11_SRC_skim_total.root", "READ");
		}
		else if (user_target == "C12")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/c12_SRC_skim_total.root", "READ");
		}
		else if (user_target == "Ca40")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/ca40_SRC_skim_total.root", "READ");
		}
		else if (user_target == "Ca48")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/ca48_SRC_skim_total.root", "READ");
		}
		else if (user_target == "Fe54")
		{
			inROOT = new TFile("../../../../../../../../../../volatile/hallc/c-cafe-2022/OFFLINE/PASS1/CAFE_OUTPUT/skimmed_pass1/fe54_SRC_skim_total.root", "READ");
		}
		else
		{
			cout << "Couldn't find file." << endl;
		}
	}
	else
	{
		cout << "Couldn't find file." << endl;
	}
	

	TTree* inputtree = 0;
	//if (user_target == "Ca48" && user_runtype == "MF")
	//{
	//	cout << 4 << endl;
	//	inputtree = (TChain*)inROOT->Get("T");								//Get the data tree
	//	cout << 5 << endl;
	//}
	//else
	//{
	inputtree = (TTree*)inROOT->Get("T");								//Get the data tree
	//}

	//TFile newfile("good.root", "RECREATE");									//Created a new file to write reduced events to
	//auto good = inputtree->CloneTree(0);										//Clone the relevant tree


	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	//------------------------------------------------------------------------DECLARE HISTOGRAMS------------------------------------------------------------------------
	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	TList* PIDList = new TList();		//Create TList to store histograms
	TList* PKList = new TList();		//Create TList to store histograms
	TList* SKList = new TList();		//Create TList to store histograms
	TList* HAccList = new TList();		//Create TList to store histograms
	TList* SAccList = new TList();		//Create TList to store histograms
	TList* HList2 = new TList();		//Create TList to store histograms

	//double hhod_GoodSinHit; inputtree->SetBranchAddress("H.hod.goodscinhit", &hhod_GoodScinHit);
	//double hdc_ntrack; inputtree->SetBranchAddress("H.dc.ntrack", &hdc_ntrack);
	//double phod_GoodSinHit; inputtree->SetBranchAddress("P.hod.goodscinhit", &phod_GoodScinHit);
	//double pdc_ntrack; inputtree->SetBranchAddress("P.dc.ntrack", &pdc_ntrack);

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//------------------------------------------------------------------------PID Histograms Bins-----------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	bool hmscoll; inputtree->SetBranchAddress("hms_collimator_cut_flag", &hmscoll);
	bool shmscoll; inputtree->SetBranchAddress("shms_collimator_cut_flag", &shmscoll);

	//Find offset for ep_ctime
	//Double_t deadTime;  inputtree->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw", &deadTime);
	//Double_t ep_ctime;	inputtree->SetBranchAddress("CTime.epCoinTime_ROC2", &ep_ctime);
	//TH1F* H1_ctime_peak = new TH1F("H1_ctime_peak", "", 200, -100, 100);
	//for (int i = 0; i <= 50000; i++) 
	//{
	//	inputtree->GetEntry(i);
	//	H1_ctime_peak->Fill(ep_ctime);
	//}
	//int binmax_ctime = H1_ctime_peak->GetMaximumBin();
	//double ctime_offset = H1_ctime_peak->GetXaxis()->GetBinCenter(binmax_ctime);

	Double_t ep_ctime;
	TH1F* H1_ep_ctime = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime", PIDList, binning);						//Four-momentum trasfer
	TH1F* H1_ep_ctime_pid = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid", PIDList, binning);						//Four-momentum trasfer
	TH1F* H1_ep_ctime_pid_acc = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc", PIDList, binning);						//Four-momentum trasfer
	TH1F* H1_ep_ctime_pid_acc_kin = Histo::mk_epctime(inputtree, ep_ctime, "H1_ep_ctime_pid_acc_kin", PIDList, binning);						//Four-momentum trasfer

	//----- HMS -----
	Double_t hCerNpeSum;			TH1F* H1_hCerNpeSum =			Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum", PIDList, binning);					//Four-momentum trasfer
	Double_t hCalEtotNorm;			TH1F* H1_hCalEtotNorm =			Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm", PIDList, binning);			//Four-momentum trasfer
	Double_t hCalEtotTrkNorm;		TH1F* H1_hCalEtotTrkNorm =		Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm", PIDList, binning);	//Four-momentum trasfer
	Double_t hHodBetaNtrk;			TH1F* H1_hHodBetaNtrk =			Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk", PIDList, binning);			//Four-momentum trasfer
	//Double_t hHodBetaTrk;			TH1F* H1_hHodBetaTrk =			Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk", PIDList, binning);				//Four-momentum trasfer
	TH1F* H1_hCerNpeSum_pid =			Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum_pid", PIDList, binning);					//Four-momentum trasfer
	TH1F* H1_hCalEtotNorm_pid =		Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm_pid", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_hCalEtotTrkNorm_pid =	Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm_pid", PIDList, binning);	//Four-momentum trasfer
	TH1F* H1_hHodBetaNtrk_pid =		Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk_pid", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_hHodBetaTrk_pid =		Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk_pid", PIDList, binning);				//Four-momentum trasfer
	TH1F* H1_hCerNpeSum_pid_acc = Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum_pid_acc", PIDList, binning);					//Four-momentum trasfer
	TH1F* H1_hCalEtotNorm_pid_acc = Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm_pid_acc", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_hCalEtotTrkNorm_pid_acc = Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm_pid_acc", PIDList, binning);	//Four-momentum trasfer
	TH1F* H1_hHodBetaNtrk_pid_acc = Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk_pid_acc", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_hHodBetaTrk_pid_acc = Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk_pid_acc", PIDList, binning);				//Four-momentum trasfer
	TH1F* H1_hCerNpeSum_pid_acc_kin = Histo::mk_hCerNpeSum(inputtree, hCerNpeSum, "H1_hCerNpeSum_pid_acc_kin", PIDList, binning);					//Four-momentum trasfer
	TH1F* H1_hCalEtotNorm_pid_acc_kin = Histo::mk_hCalEtotNorm(inputtree, hCalEtotNorm, "H1_hCalEtotNorm_pid_acc_kin", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_hCalEtotTrkNorm_pid_acc_kin = Histo::mk_hCalEtotTrkNorm(inputtree, hCalEtotTrkNorm, "H1_hCalEtotTrkNorm_pid_acc_kin", PIDList, binning);	//Four-momentum trasfer
	TH1F* H1_hHodBetaNtrk_pid_acc_kin = Histo::mk_hHodBetaNtrk(inputtree, hHodBetaNtrk, "H1_hHodBetaNtrk_pid_acc_kin", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_hHodBetaTrk_pid_acc_kin = Histo::mk_hHodBetaTrk(inputtree, hHodBetaTrk, "H1_hHodBetaTrk_pid_acc_kin", PIDList, binning);				//Four-momentum trasfer

	//----- SHMS -----
	Double_t pNGCerNpeSum;			TH1F* H1_pNGCerNpeSum =			Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum", PIDList, binning);			//Four-momentum trasfer
	Double_t pHGCerNpeSum;			TH1F* H1_pHGCerNpeSum =			Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum", PIDList, binning);			//Four-momentum trasfer
	Double_t pCalEtotNorm;			TH1F* H1_pCalEtotNorm =			Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm", PIDList, binning);			//Four-momentum trasfer
	Double_t pCalEtotTrkNorm;		TH1F* H1_pCalEtotTrkNorm =		Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm", PIDList, binning);	//Four-momentum trasfer
	Double_t pHodBetaNtrk;			TH1F* H1_pHodBetaNtrk =			Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk", PIDList, binning);			//Four-momentum trasfer
	//Double_t pHodBetaTrk;			TH1F* H1_pHodBetaTrk =			Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk", PIDList, binning);				//Four-momentum trasfer
	TH1F* H1_pNGCerNpeSum_pid =		Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum_pid", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_pHGCerNpeSum_pid =		Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum_pid", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_pCalEtotNorm_pid =		Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm_pid", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_pCalEtotTrkNorm_pid =	Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_pid", PIDList, binning);	//Four-momentum trasfer
	TH1F* H1_pHodBetaNtrk_pid =		Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk_pid", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_pHodBetaTrk_pid =		Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk_pid", PIDList, binning);				//Four-momentum trasfer
	TH1F* H1_pNGCerNpeSum_pid_acc = Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum_pid_acc", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_pHGCerNpeSum_pid_acc = Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum_pid_acc", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_pCalEtotNorm_pid_acc = Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm_pid_acc", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_pCalEtotTrkNorm_pid_acc = Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_pid_acc", PIDList, binning);	//Four-momentum trasfer
	TH1F* H1_pHodBetaNtrk_pid_acc = Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk_pid_acc", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_pHodBetaTrk_pid_acc = Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk_pid_acc", PIDList, binning);				//Four-momentum trasfer
	TH1F* H1_pNGCerNpeSum_pid_acc_kin = Histo::mk_pNGCerNpeSum(inputtree, pNGCerNpeSum, "H1_pNGCerNpeSum_pid_acc_kin", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_pHGCerNpeSum_pid_acc_kin = Histo::mk_pHGCerNpeSum(inputtree, pHGCerNpeSum, "H1_pHGCerNpeSum_pid_acc_kin", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_pCalEtotNorm_pid_acc_kin = Histo::mk_pCalEtotNorm(inputtree, pCalEtotNorm, "H1_pCalEtotNorm_pid_acc_kin", PIDList, binning);			//Four-momentum trasfer
	TH1F* H1_pCalEtotTrkNorm_pid_acc_kin = Histo::mk_pCalEtotTrkNorm(inputtree, pCalEtotTrkNorm, "H1_pCalEtotTrkNorm_pid_acc_kin", PIDList, binning);	//Four-momentum trasfer
	TH1F* H1_pHodBetaNtrk_pid_acc_kin = Histo::mk_pHodBetaNtrk(inputtree, pHodBetaNtrk, "H1_pHodBetaNtrk_pid_acc_kin", PIDList, binning);			//Four-momentum trasfer
	//TH1F* H1_pHodBetaTrk_pid_acc_kin = Histo::mk_pHodBetaTrk(inputtree, pHodBetaTrk, "H1_pHodBetaTrk_pid_acc_kin", PIDList, binning);				//Four-momentum trasfer

	//******************************************************************************************************************************************************************
	//-------------------------------------------------Primary Kinematics (electron kinematics) (USED BY DATA AND SIMC)-------------------------------------------------
	//******************************************************************************************************************************************************************
	
	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	Double_t th_e;				TH1F* H1_th_e =						Histo::mk_the(inputtree, th_e, "H1_th_e", PKList, binning);			//Electron scattering angle
	Double_t W;					TH1F* H1_W =						Histo::mk_W(inputtree, W, "H1_W", PKList, binning);					//Invariant mass
	Double_t Q2;				TH1F* H1_Q2 =						Histo::mk_Q2(inputtree, Q2, "H1_Q2", PKList, binning);				//Four-momentum trasfer
	Double_t x_bj;				TH1F* H1_xbj =						Histo::mk_xbj(inputtree, x_bj, "H1_xbj", PKList, binning);			//B-jorken X  scaling variable
	Double_t nu;				TH1F* H1_nu =						Histo::mk_nu(inputtree, nu, "H1_nu", PKList, binning);				//Energy Transfer
	Double_t q;					TH1F* H1_q =						Histo::mk_q(inputtree, q, "H1_q", PKList, binning);					//Magnitude of the 3-vector q
	Double_t q_x;				TH1F* H1_qx =						Histo::mk_qx(inputtree, q_x, "H1_qx", PKList, binning);				//x-component of the energy transfer
	Double_t q_y;				TH1F* H1_qy =						Histo::mk_qy(inputtree, q_y, "H1_qy", PKList, binning);				//y-component of the energy transfer
	Double_t q_z;				TH1F* H1_qz =						Histo::mk_qz(inputtree, q_z, "H1_qz", PKList, binning);				//z-component of the energy transfer
	Double_t th_q;				TH1F* H1_th_q =						Histo::mk_thq(inputtree, th_q, "H1_th_q", PKList, binning);			//Angle between q and +z (hall coord. system)
	Double_t ph_q;				TH1F* H1_ph_q =						Histo::mk_phq(inputtree, ph_q, "H1_ph_q", PKList, binning);			//Out of plane angle between beamline and q
	//Double_t omega;				TH1F* H1_omega =					Histo::mk_omega(inputtree, omega, "H1_omega", PKList, binning);		//
	
	TH1F* H1_th_e_pid =			Histo::mk_the(inputtree, th_e, "H1_th_e_pid", PKList, binning);			//Electron scattering angle
	TH1F* H1_W_pid =			Histo::mk_W(inputtree, W, "H1_W_pid", PKList, binning);					//Invariant mass
	TH1F* H1_Q2_pid =			Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid", PKList, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid =			Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid", PKList, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid =			Histo::mk_nu(inputtree, nu, "H1_nu_pid", PKList, binning);				//Energy Transfer
	TH1F* H1_q_pid =			Histo::mk_q(inputtree, q, "H1_q_pid", PKList, binning);				//Magnitude of the 3-vector q
	TH1F* H1_qx_pid =			Histo::mk_qx(inputtree, q_x, "H1_qx_pid", PKList, binning);				//x-component of the energy transfer
	TH1F* H1_qy_pid =			Histo::mk_qy(inputtree, q_y, "H1_qy_pid", PKList, binning);				//y-component of the energy transfer
	TH1F* H1_qz_pid =			Histo::mk_qz(inputtree, q_z, "H1_qz_pid", PKList, binning);				//z-component of the energy transfer
	TH1F* H1_th_q_pid =			Histo::mk_thq(inputtree, th_q, "H1_th_q_pid", PKList, binning);			//Angle between q and +z (hall coord. system)
	TH1F* H1_ph_q_pid =			Histo::mk_phq(inputtree, ph_q, "H1_ph_q_pid", PKList, binning);			//Out of plane angle between beamline and q
	//TH1F* H1_omega_pid =		Histo::mk_omega(inputtree, omega, "H1_omega_pid", PKList, binning);		//

	TH1F* H1_th_e_pid_acc =		Histo::mk_the(inputtree, th_e, "H1_th_e_pid_acc", PKList, binning);			//Electron scattering angle
	TH1F* H1_W_pid_acc =		Histo::mk_W(inputtree, W, "H1_W_pid_acc", PKList, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc =		Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc", PKList, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc =		Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc", PKList, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc =		Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc", PKList, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc =		Histo::mk_q(inputtree, q, "H1_q_pid_acc", PKList, binning);				//Magnitude of the 3-vector q
	TH1F* H1_qx_pid_acc =		Histo::mk_qx(inputtree, q_x, "H1_qx_pid_acc", PKList, binning);				//x-component of the energy transfer
	TH1F* H1_qy_pid_acc =		Histo::mk_qy(inputtree, q_y, "H1_qy_pid_acc", PKList, binning);				//y-component of the energy transfer
	TH1F* H1_qz_pid_acc =		Histo::mk_qz(inputtree, q_z, "H1_qz_pid_acc", PKList, binning);				//z-component of the energy transfer
	TH1F* H1_th_q_pid_acc =		Histo::mk_thq(inputtree, th_q, "H1_th_q_pid_acc", PKList, binning);			//Angle between q and +z (hall coord. system)
	TH1F* H1_ph_q_pid_acc =		Histo::mk_phq(inputtree, ph_q, "H1_ph_q_pid_acc", PKList, binning);			//Out of plane angle between beamline and q
	//TH1F* H1_omega_pid_acc =	Histo::mk_omega(inputtree, omega, "H1_omega_pid_acc", PKList, binning);		//

	TH1F* H1_th_e_pid_acc_kin =	Histo::mk_the(inputtree, th_e, "H1_th_e_pid_acc_kin", PKList, binning);			//Electron scattering angle
	TH1F* H1_W_pid_acc_kin =	Histo::mk_W(inputtree, W, "H1_W_pid_acc_kin", PKList, binning);					//Invariant mass
	TH1F* H1_Q2_pid_acc_kin =	Histo::mk_Q2(inputtree, Q2, "H1_Q2_pid_acc_kin", PKList, binning);				//Four-momentum trasfer
	TH1F* H1_xbj_pid_acc_kin =	Histo::mk_xbj(inputtree, x_bj, "H1_xbj_pid_acc_kin", PKList, binning);			//B-jorken X  scaling variable
	TH1F* H1_nu_pid_acc_kin =	Histo::mk_nu(inputtree, nu, "H1_nu_pid_acc_kin", PKList, binning);				//Energy Transfer
	TH1F* H1_q_pid_acc_kin =	Histo::mk_q(inputtree, q, "H1_q_pid_acc_kin", PKList, binning);				//Magnitude of the 3-vector q
	TH1F* H1_qx_pid_acc_kin =	Histo::mk_qx(inputtree, q_x, "H1_qx_pid_acc_kin", PKList, binning);				//x-component of the energy transfer
	TH1F* H1_qy_pid_acc_kin =	Histo::mk_qy(inputtree, q_y, "H1_qy_pid_acc_kin", PKList, binning);				//y-component of the energy transfer
	TH1F* H1_qz_pid_acc_kin =	Histo::mk_qz(inputtree, q_z, "H1_qz_pid_acc_kin", PKList, binning);				//z-component of the energy transfer
	TH1F* H1_th_q_pid_acc_kin =	Histo::mk_thq(inputtree, th_q, "H1_th_q_pid_acc_kin", PKList, binning);			//Angle between q and +z (hall coord. system)
	TH1F* H1_ph_q_pid_acc_kin = Histo::mk_phq(inputtree, ph_q, "H1_ph_q_pid_acc_kin", PKList, binning);			//Out of plane angle between beamline and q
	//TH1F* H1_omega_pid_acc_kin = Histo::mk_omega(inputtree, omega, "H1_omega_pid_acc_kin", PKList, binning);	//

	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	//Double_t Ki;			TH1F* H1_Ki = Histo::mk_ki(inputtree, Ki, "H1_Ki", SKList);
	//Double_t Kf;			TH1F* H1_Kf = Histo::mk_kf(inputtree, Kf, "H1_Kf", SKList);
	

	//******************************************************************************************************************************************************************
	//--------------------------------------------Secondary (Hadron) Kinematics (recoil and missing are used interchageably)--------------------------------------------
	//******************************************************************************************************************************************************************
	
	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	//Double_t emiss;			TH1F* H1_emiss =			Histo::mk_Em(inputtree, emiss, "H1_emiss", SKList, binning);			//Standard Missing Energy for H(e,e'p)
	Double_t emiss_nuc;		TH1F* H1_emiss_nuc =		Histo::mk_Em_nuc(inputtree, emiss_nuc, "H1_emiss_nuc", SKList, binning);
	Double_t pmiss;			TH1F* H1_pmiss =			Histo::mk_Pm(inputtree, pmiss, "H1_pmiss", SKList, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	Double_t prec_x;		TH1F* H1_prec_x =			Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x", SKList, binning); 
	Double_t prec_y;		TH1F* H1_prec_y =			Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y", SKList, binning); 
	Double_t prec_z;		TH1F* H1_prec_z =			Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z", SKList, binning); 
	Double_t pmiss_x;		TH1F* H1_pmiss_x =			Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x", SKList, binning); 
	Double_t pmiss_y;		TH1F* H1_pmiss_y =			Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y", SKList, binning); 
	Double_t pmiss_z;		TH1F* H1_pmiss_z =			Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z", SKList, binning); 
	Double_t Tx;			TH1F* H1_Tx =				Histo::mk_Tx(inputtree, Tx, "H1_Tx", SKList, binning);
	Double_t Tr;			TH1F* H1_Tr =				Histo::mk_Tr(inputtree, Tr, "H1_Tr", SKList, binning);
	Double_t mmiss;			TH1F* H1_mmiss =			Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss", SKList, binning);
	Double_t th_pq;			TH1F* H1_th_pq =			Histo::mk_thpq(inputtree, th_pq, "H1_th_pq", SKList, binning); 		//detected particle in-plane angle w.r.to q-vector
	Double_t th_rq;			TH1F* H1_th_rq =			Histo::mk_thrq(inputtree, th_rq, "H1_th_rq", SKList, binning); 		//recoil particle in-plane angle w.r.to q-vector
	Double_t cth_rq;		TH1F* H1_cth_rq =			Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq", SKList, binning); 		//recoil particle in-plane angle w.r.to q-vector
	Double_t ph_pq;			TH1F* H1_ph_pq =			Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq", SKList, binning); 		//detected particle ???
	Double_t ph_rq;			TH1F* H1_ph_rq =			Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq", SKList, binning); 
	Double_t xangle;		TH1F* H1_xangle =			Histo::mk_xangle(inputtree, xangle, "H1_xangle", SKList, binning); 
	
	//TH1F* H1_emiss_pid =		Histo::mk_Em(inputtree, emiss, "H1_emiss_pid", SKList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_emiss_nuc_pid =	Histo::mk_Em_nuc(inputtree, emiss_nuc, "H1_emiss_nuc_pid", SKList, binning);
	TH1F* H1_pmiss_pid =		Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid", SKList, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_prec_x_pid =		Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x_pid", SKList, binning);
	TH1F* H1_prec_y_pid =		Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y_pid", SKList, binning);
	TH1F* H1_prec_z_pid =		Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z_pid", SKList, binning);
	TH1F* H1_pmiss_x_pid =		Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x_pid", SKList, binning);
	TH1F* H1_pmiss_y_pid =		Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y_pid", SKList, binning);
	TH1F* H1_pmiss_z_pid =		Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z_pid", SKList, binning);
	TH1F* H1_Tx_pid =			Histo::mk_Tx(inputtree, Tx, "H1_Tx_pid", SKList, binning);
	TH1F* H1_Tr_pid =			Histo::mk_Tr(inputtree, Tr, "H1_Tr_pid", SKList, binning);
	TH1F* H1_mmiss_pid =		Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss_pid", SKList, binning);
	TH1F* H1_th_pq_pid =		Histo::mk_thpq(inputtree, th_pq, "H1_th_pq_pid", SKList, binning); 		//detected particle in-plane angle w.r.to q-vector
	TH1F* H1_th_rq_pid =		Histo::mk_thrq(inputtree, th_rq, "H1_th_rq_pid", SKList, binning); 		//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_cth_rq_pid =		Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid", SKList, binning); 		//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_ph_pq_pid =		Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq_pid", SKList, binning); 		//detected particle ???
	TH1F* H1_ph_rq_pid =		Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq_pid", SKList, binning);
	TH1F* H1_xangle_pid =		Histo::mk_xangle(inputtree, xangle, "H1_xangle_pid", SKList, binning);
	
	//TH1F* H1_emiss_pid_acc = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc", SKList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_emiss_nuc_pid_acc = Histo::mk_Em_nuc(inputtree, emiss_nuc, "H1_emiss_nuc_pid_acc", SKList, binning);
	TH1F* H1_pmiss_pid_acc = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc", SKList, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_prec_x_pid_acc = Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x_pid_acc", SKList, binning);
	TH1F* H1_prec_y_pid_acc = Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y_pid_acc", SKList, binning);
	TH1F* H1_prec_z_pid_acc = Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z_pid_acc", SKList, binning);
	TH1F* H1_pmiss_x_pid_acc = Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x_pid_acc", SKList, binning);
	TH1F* H1_pmiss_y_pid_acc = Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y_pid_acc", SKList, binning);
	TH1F* H1_pmiss_z_pid_acc = Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z_pid_acc", SKList, binning);
	TH1F* H1_Tx_pid_acc = Histo::mk_Tx(inputtree, Tx, "H1_Tx_pid_acc", SKList, binning);
	TH1F* H1_Tr_pid_acc = Histo::mk_Tr(inputtree, Tr, "H1_Tr_pid_acc", SKList, binning);
	TH1F* H1_mmiss_pid_acc = Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss_pid_acc", SKList, binning);
	TH1F* H1_th_pq_pid_acc = Histo::mk_thpq(inputtree, th_pq, "H1_th_pq_pid_acc", SKList, binning); 		//detected particle in-plane angle w.r.to q-vector
	TH1F* H1_th_rq_pid_acc = Histo::mk_thrq(inputtree, th_rq, "H1_th_rq_pid_acc", SKList, binning); 		//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_cth_rq_pid_acc = Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid_acc", SKList, binning); 		//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_ph_pq_pid_acc = Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq_pid_acc", SKList, binning); 		//detected particle ???
	TH1F* H1_ph_rq_pid_acc = Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq_pid_acc", SKList, binning);
	TH1F* H1_xangle_pid_acc = Histo::mk_xangle(inputtree, xangle, "H1_xangle_pid_acc", SKList, binning);

	//TH1F* H1_emiss_pid_acc_kin = Histo::mk_Em(inputtree, emiss, "H1_emiss_pid_acc_kin", SKList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_emiss_nuc_pid_acc_kin = Histo::mk_Em_nuc(inputtree, emiss_nuc, "H1_emiss_nuc_pid_acc_kin", SKList, binning);
	TH1F* H1_pmiss_pid_acc_kin = Histo::mk_Pm(inputtree, pmiss, "H1_pmiss_pid_acc_kin", SKList, binning);			//Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
	TH1F* H1_prec_x_pid_acc_kin = Histo::mk_Pmx_lab(inputtree, prec_x, "H1_prec_x_pid_acc_kin", SKList, binning);
	TH1F* H1_prec_y_pid_acc_kin = Histo::mk_Pmy_lab(inputtree, prec_y, "H1_prec_y_pid_acc_kin", SKList, binning);
	TH1F* H1_prec_z_pid_acc_kin = Histo::mk_Pmz_lab(inputtree, prec_z, "H1_prec_z_pid_acc_kin", SKList, binning);
	TH1F* H1_pmiss_x_pid_acc_kin = Histo::mk_Pmx_q(inputtree, pmiss_x, "H1_pmiss_x_pid_acc_kin", SKList, binning);
	TH1F* H1_pmiss_y_pid_acc_kin = Histo::mk_Pmy_q(inputtree, pmiss_y, "H1_pmiss_y_pid_acc_kin", SKList, binning);
	TH1F* H1_pmiss_z_pid_acc_kin = Histo::mk_Pmz_q(inputtree, pmiss_z, "H1_pmiss_z_pid_acc_kin", SKList, binning);
	TH1F* H1_Tx_pid_acc_kin = Histo::mk_Tx(inputtree, Tx, "H1_Tx_pid_acc_kin", SKList, binning);
	TH1F* H1_Tr_pid_acc_kin = Histo::mk_Tr(inputtree, Tr, "H1_Tr_pid_acc_kin", SKList, binning);
	TH1F* H1_mmiss_pid_acc_kin = Histo::mk_Mrecoil(inputtree, mmiss, "H1_mmiss_pid_acc_kin", SKList, binning);
	TH1F* H1_th_pq_pid_acc_kin = Histo::mk_thpq(inputtree, th_pq, "H1_th_pq_pid_acc_kin", SKList, binning); 		//detected particle in-plane angle w.r.to q-vector
	TH1F* H1_th_rq_pid_acc_kin = Histo::mk_thrq(inputtree, th_rq, "H1_th_rq_pid_acc_kin", SKList, binning); 		//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_cth_rq_pid_acc_kin = Histo::mk_cthrq(inputtree, cth_rq, "H1_cth_rq_pid_acc_kin", SKList, binning); 		//recoil particle in-plane angle w.r.to q-vector
	TH1F* H1_ph_pq_pid_acc_kin = Histo::mk_phpq(inputtree, ph_pq, "H1_ph_pq_pid_acc_kin", SKList, binning); 		//detected particle ???
	TH1F* H1_ph_rq_pid_acc_kin = Histo::mk_phrq(inputtree, ph_rq, "H1_ph_rq_pid_acc_kin", SKList, binning);
	TH1F* H1_xangle_pid_acc_kin = Histo::mk_xangle(inputtree, xangle, "H1_xangle_pid_acc_kin", SKList, binning);

	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	//Double_t Pf;			TH1F* H1_Pf =				Histo::mk_Pf(inputtree, Pf, "H1_Pf", SKList, binning);
	
	Double_t th_p;			TH1F* H1_th_p = Histo::mk_thp(inputtree, th_p, "H1_th_p", SKList, binning);
	TH1F* H1_th_p_pid = Histo::mk_thp(inputtree, th_p, "H1_th_p_pid", SKList, binning);
	TH1F* H1_th_p_pid_acc = Histo::mk_thp(inputtree, th_p, "H1_th_p_pid_acc", SKList, binning);
	TH1F* H1_th_p_pid_acc_kin = Histo::mk_thp(inputtree, th_p, "H1_th_p_pid_acc_kin", SKList, binning);

	//mmiss2

	Double_t Em_src;		TH1F* H1_Em_src = Histo::mk_Em_src(inputtree, Em_src, "H1_Em_src", SKList, binning);

	

	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------------Acceptance Histogram Bins---------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//----------------------------------------
	//--------------- HMS Leafs --------------
	//----------------------------------------
	//----- Hadron Arm Focal Plane -----
	Double_t hxfp;			TH1F* H1_hxfp =						Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp", HAccList, binning); 			//Standard Missing Energy for H(e,e'p)
	Double_t hxpfp;			TH1F* H1_hxpfp =					Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	Double_t hyfp;			TH1F* H1_hyfp =						Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp", HAccList, binning); 			//Standard Missing Energy for H(e,e'p)
	Double_t hypfp;			TH1F* H1_hypfp =					Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hxfp_pid = Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp_pid", HAccList, binning); 			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hxpfp_pid = Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp_pid", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hyfp_pid = Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp_pid", HAccList, binning); 			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hypfp_pid = Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp_pid", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hxfp_pid_acc = Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp_pid_acc", HAccList, binning); 			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hxpfp_pid_acc = Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp_pid_acc", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hyfp_pid_acc = Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp_pid_acc", HAccList, binning); 			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hypfp_pid_acc = Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp_pid_acc", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hxfp_pid_acc_kin = Histo::mk_hxfp(inputtree, hxfp, "H1_hxfp_pid_acc_kin", HAccList, binning); 			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hxpfp_pid_acc_kin = Histo::mk_hxpfp(inputtree, hxpfp, "H1_hxpfp_pid_acc_kin", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hyfp_pid_acc_kin = Histo::mk_hyfp(inputtree, hyfp, "H1_hyfp_pid_acc_kin", HAccList, binning); 			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hypfp_pid_acc_kin = Histo::mk_hypfp(inputtree, hypfp, "H1_hypfp_pid_acc_kin", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)

	//----- Hadron Arm Reconstructed Quantities -----
	Double_t hytar;			TH1F* H1_hytar =					Histo::mk_hytar(inputtree, hytar, "H1_hytar", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	Double_t hyptar;		TH1F* H1_hyptar =					Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	Double_t hxptar;		TH1F* H1_hxptar =					Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	Double_t hdelta;		TH1F* H1_hdelta =					Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hytar_pid = Histo::mk_hytar(inputtree, hytar, "H1_hytar_pid", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hyptar_pid = Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar_pid", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hxptar_pid = Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar_pid", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hdelta_pid = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_pid", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hytar_pid_acc = Histo::mk_hytar(inputtree, hytar, "H1_hytar_pid_acc", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hyptar_pid_acc = Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar_pid_acc", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hxptar_pid_acc = Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar_pid_acc", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hdelta_pid_acc = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_pid_acc", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hytar_pid_acc_kin = Histo::mk_hytar(inputtree, hytar, "H1_hytar_pid_acc_kin", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hyptar_pid_acc_kin = Histo::mk_hyptar(inputtree, hyptar, "H1_hyptar_pid_acc_kin", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hxptar_pid_acc_kin = Histo::mk_hxptar(inputtree, hxptar, "H1_hxptar_pid_acc_kin", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hdelta_pid_acc_kin = Histo::mk_hdelta(inputtree, hdelta, "H1_hdelta_pid_acc_kin", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	
	//----- Target Reconstruction (Hall Coord. System) -----
	Double_t htarx;			TH1F* H1_htarx =					Histo::mk_htarx(inputtree, htarx, "H1_htarx", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	Double_t htary;			TH1F* H1_htary =					Histo::mk_htary(inputtree, htary, "H1_htary", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	Double_t htarz;			TH1F* H1_htarz =					Histo::mk_htarz(inputtree, htarz, "H1_htarz", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_htarx_pid = Histo::mk_htarx(inputtree, htarx, "H1_htarx_pid", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_htary_pid = Histo::mk_htary(inputtree, htary, "H1_htary_pid", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_htarz_pid = Histo::mk_htarz(inputtree, htarz, "H1_htarz_pid", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_htarx_pid_acc = Histo::mk_htarx(inputtree, htarx, "H1_htarx_pid_acc", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_htary_pid_acc = Histo::mk_htary(inputtree, htary, "H1_htary_pid_acc", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_htarz_pid_acc = Histo::mk_htarz(inputtree, htarz, "H1_htarz_pid_acc", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_htarx_pid_acc_kin = Histo::mk_htarx(inputtree, htarx, "H1_htarx_pid_acc_kin", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_htary_pid_acc_kin = Histo::mk_htary(inputtree, htary, "H1_htary_pid_acc_kin", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_htarz_pid_acc_kin = Histo::mk_htarz(inputtree, htarz, "H1_htarz_pid_acc_kin", HAccList, binning); 		//Standard Missing Energy for H(e,e'p)
	
	//----- HMS Collimator -----
	Double_t hXColl;		TH1F* H1_hXColl =					Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	Double_t hYColl;		TH1F* H1_hYColl =					Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hXColl_pid = Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl_pid", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hYColl_pid = Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl_pid", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hXColl_pid_acc = Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl_pid_acc", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hYColl_pid_acc = Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl_pid_acc", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_hXColl_pid_acc_kin = Histo::mk_hXColl(inputtree, hXColl, "H1_hXColl_pid_acc_kin", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)
	TH1F* H1_hYColl_pid_acc_kin = Histo::mk_hYColl(inputtree, hYColl, "H1_hYColl_pid_acc_kin", HAccList, binning); 	//Standard Missing Energy for H(e,e'p)


	//----------------------------------------
	//-------------- SHMS Leafs --------------
	//----------------------------------------
	//----- Electron Arm Focal Plane -----
	Double_t exfp;			TH1F* H1_exfp =						Histo::mk_exfp(inputtree, exfp, "H1_exfp", SAccList, binning);			//Standard Missing Energy for H(e,e'p)
	Double_t expfp;			TH1F* H1_expfp =					Histo::mk_expfp(inputtree, expfp, "H1_expfp", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	Double_t eyfp;			TH1F* H1_eyfp =						Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp", SAccList, binning);			//Standard Missing Energy for H(e,e'p)
	Double_t eypfp;			TH1F* H1_eypfp =					Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_exfp_pid =		Histo::mk_exfp(inputtree, exfp, "H1_exfp_pid", SAccList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_expfp_pid =		Histo::mk_expfp(inputtree, expfp, "H1_expfp_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eyfp_pid =		Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp_pid", SAccList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eypfp_pid =		Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_exfp_pid_acc = Histo::mk_exfp(inputtree, exfp, "H1_exfp_pid_acc", SAccList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_expfp_pid_acc = Histo::mk_expfp(inputtree, expfp, "H1_expfp_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eyfp_pid_acc = Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp_pid_acc", SAccList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eypfp_pid_acc = Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_exfp_pid_acc_kin = Histo::mk_exfp(inputtree, exfp, "H1_exfp_pid_acc_kin", SAccList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_expfp_pid_acc_kin = Histo::mk_expfp(inputtree, expfp, "H1_expfp_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eyfp_pid_acc_kin = Histo::mk_eyfp(inputtree, eyfp, "H1_eyfp_pid_acc_kin", SAccList, binning);			//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eypfp_pid_acc_kin = Histo::mk_eypfp(inputtree, eypfp, "H1_eypfp_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)

	//----- Electron Arm Reconstructed Quantities -----
	Double_t eytar;			TH1F* H1_eytar =					Histo::mk_eytar(inputtree, eytar, "H1_eytar", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	Double_t eyptar;		TH1F* H1_eyptar =					Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	Double_t exptar;		TH1F* H1_exptar =					Histo::mk_exptar(inputtree, exptar, "H1_exptar", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	Double_t edelta;		TH1F* H1_edelta =					Histo::mk_edelta(inputtree, edelta, "H1_edelta", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_eytar_pid =		Histo::mk_eytar(inputtree, eytar, "H1_eytar_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eyptar_pid =		Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_exptar_pid =		Histo::mk_exptar(inputtree, exptar, "H1_exptar_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_edelta_pid =		Histo::mk_edelta(inputtree, edelta, "H1_edelta_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_eytar_pid_acc = Histo::mk_eytar(inputtree, eytar, "H1_eytar_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eyptar_pid_acc = Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_exptar_pid_acc = Histo::mk_exptar(inputtree, exptar, "H1_exptar_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_edelta_pid_acc = Histo::mk_edelta(inputtree, edelta, "H1_edelta_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_eytar_pid_acc_kin = Histo::mk_eytar(inputtree, eytar, "H1_eytar_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eyptar_pid_acc_kin = Histo::mk_eyptar(inputtree, eyptar, "H1_eyptar_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_exptar_pid_acc_kin = Histo::mk_exptar(inputtree, exptar, "H1_exptar_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_edelta_pid_acc_kin = Histo::mk_edelta(inputtree, edelta, "H1_edelta_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	//----- Target Reconstruction (Hall Coord. System) -----
	Double_t etarx;			TH1F* H1_etarx =					Histo::mk_etarx(inputtree, etarx, "H1_etarx", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	Double_t etary;			TH1F* H1_etary =					Histo::mk_etary(inputtree, etary, "H1_etary", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	Double_t etarz;			TH1F* H1_etarz =					Histo::mk_etarz(inputtree, etarz, "H1_etarz", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_etarx_pid =		Histo::mk_etarx(inputtree, etarx, "H1_etarx_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_etary_pid =		Histo::mk_etary(inputtree, etary, "H1_etary_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_etarz_pid =		Histo::mk_etarz(inputtree, etarz, "H1_etarz_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_etarx_pid_acc = Histo::mk_etarx(inputtree, etarx, "H1_etarx_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_etary_pid_acc = Histo::mk_etary(inputtree, etary, "H1_etary_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_etarz_pid_acc = Histo::mk_etarz(inputtree, etarz, "H1_etarz_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_etarx_pid_acc_kin = Histo::mk_etarx(inputtree, etarx, "H1_etarx_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_etary_pid_acc_kin = Histo::mk_etary(inputtree, etary, "H1_etary_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_etarz_pid_acc_kin = Histo::mk_etarz(inputtree, etarz, "H1_etarz_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	//----- HMS Collimator -----
	Double_t eXColl;		TH1F* H1_eXColl =					Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	Double_t eYColl;		TH1F* H1_eYColl =					Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_eXColl_pid =		Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eYColl_pid =		Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_eXColl_pid_acc = Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eYColl_pid_acc = Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	TH1F* H1_eXColl_pid_acc_kin = Histo::mk_eXColl(inputtree, eXColl, "H1_eXColl_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_eYColl_pid_acc_kin = Histo::mk_eYColl(inputtree, eYColl, "H1_eYColl_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	
	//----- Calculated -----
	Double_t ztar_diff;		TH1F* H1_ztar_diff =				Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_ztar_diff_pid =	Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff_pid", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_ztar_diff_pid_acc = Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff_pid_acc", SAccList, binning);		//Standard Missing Energy for H(e,e'p)
	TH1F* H1_ztar_diff_pid_acc_kin = Histo::mk_ztar_diff(inputtree, ztar_diff, "H1_ztar_diff_pid_acc_kin", SAccList, binning);		//Standard Missing Energy for H(e,e'p)



	//--------------------------------------------------------
	//--------------------- 2D HISTOGRAMS --------------------
	//--------------------------------------------------------

	//----- Kinematic -----
	TH2F* H2_thrq_pmiss = Histo::mk_thrq_pmiss(inputtree, "H2_thrq_pmiss", HList2, binning);
	TH2F* H2_cthrq_pmiss = Histo::mk_cthrq_pmiss(inputtree, "H2_cthrq_pmiss", HList2, binning);
	TH2F* H2_W_pmiss = Histo::mk_W_pmiss(inputtree, "H2_W_pmiss", HList2, binning);
	TH2F* H2_xbj_Q2 = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2", HList2, binning);
	TH2F* H2_Em_Pm = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm", HList2, binning);
	//TH2F* H2_the_kf = Histo::mk_the_kf(inputtree, "H2_the_kf", HList2, binning);
	//TH2F* H2_thp_Pf = Histo::mk_thp_Pf(inputtree, "H2_thp_Pf", HList2, binning);

	TH2F* H2_thrq_pmiss_pid = Histo::mk_thrq_pmiss(inputtree, "H2_thrq_pmiss_pid", HList2, binning);
	TH2F* H2_cthrq_pmiss_pid = Histo::mk_cthrq_pmiss(inputtree, "H2_cthrq_pmiss_pid", HList2, binning);
	TH2F* H2_W_pmiss_pid = Histo::mk_W_pmiss(inputtree, "H2_W_pmiss_pid", HList2, binning);
	TH2F* H2_xbj_Q2_pid = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_pid", HList2, binning);
	TH2F* H2_Em_Pm_pid = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_pid", HList2, binning);
	
	TH2F* H2_thrq_pmiss_pid_acc = Histo::mk_thrq_pmiss(inputtree, "H2_thrq_pmiss_pid_acc", HList2, binning);
	TH2F* H2_cthrq_pmiss_pid_acc = Histo::mk_cthrq_pmiss(inputtree, "H2_cthrq_pmiss_pid_acc", HList2, binning);
	TH2F* H2_W_pmiss_pid_acc = Histo::mk_W_pmiss(inputtree, "H2_W_pmiss_pid_acc", HList2, binning);
	TH2F* H2_xbj_Q2_pid_acc = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_pid_acc", HList2, binning);
	TH2F* H2_Em_Pm_pid_acc = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_pid_acc", HList2, binning);
	
	TH2F* H2_thrq_pmiss_pid_acc_kin = Histo::mk_thrq_pmiss(inputtree, "H2_thrq_pmiss_pid_acc_kin", HList2, binning);
	TH2F* H2_cthrq_pmiss_pid_acc_kin = Histo::mk_cthrq_pmiss(inputtree, "H2_cthrq_pmiss_pid_acc_kin", HList2, binning);
	TH2F* H2_W_pmiss_pid_acc_kin = Histo::mk_W_pmiss(inputtree, "H2_W_pmiss_pid_acc_kin", HList2, binning);
	TH2F* H2_xbj_Q2_pid_acc_kin = Histo::mk_xbj_Q2(inputtree, "H2_xbj_Q2_pid_acc_kin", HList2, binning);
	TH2F* H2_Em_Pm_pid_acc_kin = Histo::mk_Em_Pm(inputtree, "H2_Em_Pm_pid_acc_kin", HList2, binning);
	
	//----- Acceptance -----
	TH2F* H2_hXColl_hYColl = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl", HList2, binning);
	TH2F* H2_eXColl_eYColl = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl", HList2, binning);
	TH2F* H2_hxfp_hyfp = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp", HList2, binning);
	TH2F* H2_exfp_eyfp = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp", HList2, binning);
	TH2F* H2_hxptar_exptar = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar", HList2, binning);
	TH2F* H2_hyptar_eyptar = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar", HList2, binning);
	TH2F* H2_hdelta_edelta = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta", HList2, binning);

	TH2F* H2_hXColl_hYColl_pid = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_pid", HList2, binning);
	TH2F* H2_eXColl_eYColl_pid = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_pid", HList2, binning);
	TH2F* H2_hxfp_hyfp_pid = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_pid", HList2, binning);
	TH2F* H2_exfp_eyfp_pid = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_pid", HList2, binning);
	TH2F* H2_hxptar_exptar_pid = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_pid", HList2, binning);
	TH2F* H2_hyptar_eyptar_pid = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_pid", HList2, binning);
	TH2F* H2_hdelta_edelta_pid = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_pid", HList2, binning);

	TH2F* H2_hXColl_hYColl_pid_acc = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_pid_acc", HList2, binning);
	TH2F* H2_eXColl_eYColl_pid_acc = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_pid_acc", HList2, binning);
	TH2F* H2_hxfp_hyfp_pid_acc = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_pid_acc", HList2, binning);
	TH2F* H2_exfp_eyfp_pid_acc = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_pid_acc", HList2, binning);
	TH2F* H2_hxptar_exptar_pid_acc = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_pid_acc", HList2, binning);
	TH2F* H2_hyptar_eyptar_pid_acc = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_pid_acc", HList2, binning);
	TH2F* H2_hdelta_edelta_pid_acc = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_pid_acc", HList2, binning);

	TH2F* H2_hXColl_hYColl_pid_acc_kin = Histo::mk_hXColl_hYColl(inputtree, "H2_hXColl_hYColl_pid_acc_kin", HList2, binning);
	TH2F* H2_eXColl_eYColl_pid_acc_kin = Histo::mk_eXColl_eYColl(inputtree, "H2_eXColl_eYColl_pid_acc_kin", HList2, binning);
	TH2F* H2_hxfp_hyfp_pid_acc_kin = Histo::mk_hxfp_hyfp(inputtree, "H2_hxfp_hyfp_pid_acc_kin", HList2, binning);
	TH2F* H2_exfp_eyfp_pid_acc_kin = Histo::mk_exfp_eyfp(inputtree, "H2_exfp_eyfp_pid_acc_kin", HList2, binning);
	TH2F* H2_hxptar_exptar_pid_acc_kin = Histo::mk_hxptar_exptar(inputtree, "H2_hxptar_exptar_pid_acc_kin", HList2, binning);
	TH2F* H2_hyptar_eyptar_pid_acc_kin = Histo::mk_hyptar_eyptar(inputtree, "H2_hyptar_eyptar_pid_acc_kin", HList2, binning);
	TH2F* H2_hdelta_edelta_pid_acc_kin = Histo::mk_hdelta_edelta(inputtree, "H2_hdelta_edelta_pid_acc_kin", HList2, binning);




	//bool noedtm = false;
	//bool epctime = false;
	bool pid_hms = false;
	bool pid_shms = false;
	bool accp_hms = false;
	bool accp_shms = false;
	bool accp_z = false;
	bool kinheepsing = false;
	bool kinheepcoin = false;
	bool kin = false;

	
	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	//-------------------------------------------------------------------------- Loop ---------------------------------------------------------------------------
	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	int numLoops = 0;							//Initialize event count to 0
	int max_loops = 50000000;					//Define maximum number of events to process
	Long64_t nentries;
	nentries = inputtree->GetEntries();		//Get the total number of entries
	cout << "There are " << nentries << " events." << endl;    

	for (int i = 0; i < nentries; i++)
	{
		inputtree->GetEntry(i);				//Get the ith entry from the T TTree

		Em_src = nu - Tx - (sqrt(MN * MN + pmiss * pmiss) - MN);
		ztar_diff = htarz - etarz;
		th_p = xangle - th_e;

		//epctime = Cuts::cTime(ep_ctime - ctime_offset);
		pid_hms = Cuts::PID_HMS(hCalEtotTrkNorm, hCerNpeSum);
		pid_shms = Cuts::PID_SHMS(pCalEtotTrkNorm, pNGCerNpeSum, pHGCerNpeSum);
		accp_hms = Cuts::Accp_HMS(hdelta, hxptar, hyptar, hmscoll);
		accp_shms = Cuts::Accp_SHMS(edelta, exptar, eyptar, shmscoll);
		accp_z = Cuts::Accp_Z(ztar_diff);
		kinheepsing = Cuts::kinHeepSing(Q2, W, x_bj);
		kinheepcoin = Cuts::kinHeepCoin(Q2, W, x_bj, emiss_nuc, mmiss);
		/*
		if (deadTime == 0)
		{
			noedtm = true;
		}
		else
		{
			noedtm = false;
		}
		*/
		if (user_runtype == "SRC")
		{
			kin = Cuts::kinSRC(Q2, pmiss, x_bj, (cth_rq / dtr), emiss_nuc, Em_src, target);
		}
		else if (user_runtype == "MF")
		{
			kin = Cuts::kinMF(Q2, pmiss, emiss_nuc, emiss_nuc, (cth_rq / dtr), target);
		}
		else
		{
			cout << "Invalid type." << endl;
		}
		

		//--------------------------------------
		//PID Histograms Bins
		//--------------------------------------
		//----- HMS -----
		H1_hCerNpeSum->Fill(hCerNpeSum);
		H1_hCalEtotNorm->Fill(hCalEtotNorm);
		H1_hCalEtotTrkNorm->Fill(hCalEtotTrkNorm);
		H1_hHodBetaNtrk->Fill(hHodBetaNtrk);
		//H1_hHodBetaTrk->Fill(hHodBetaTrk);
		H1_pNGCerNpeSum->Fill(pNGCerNpeSum);
		//H1_pHGCerNpeSum->Fill(pHGCerNpeSum);
		H1_pCalEtotNorm->Fill(pCalEtotNorm);
		H1_pCalEtotTrkNorm->Fill(pCalEtotTrkNorm);
		H1_pHodBetaNtrk->Fill(pHodBetaNtrk);
		//H1_pHodBetaTrk->Fill(pHodBetaTrk);
		//----- SHMS -----
		if (pid_hms && pid_shms)
		{
			H1_hCerNpeSum_pid->Fill(hCerNpeSum);
			H1_hCalEtotNorm_pid->Fill(hCalEtotNorm);
			H1_hCalEtotTrkNorm_pid->Fill(hCalEtotTrkNorm);
			H1_hHodBetaNtrk_pid->Fill(hHodBetaNtrk);
			//H1_hHodBetaTrk_pid->Fill(hHodBetaTrk);
			H1_pNGCerNpeSum_pid->Fill(pNGCerNpeSum);
			//H1_pHGCerNpeSum_pid->Fill(pHGCerNpeSum);
			H1_pCalEtotNorm_pid->Fill(pCalEtotNorm);
			H1_pCalEtotTrkNorm_pid->Fill(pCalEtotTrkNorm);
			H1_pHodBetaNtrk_pid->Fill(pHodBetaNtrk);
			//H1_pHodBetaTrk_pid->Fill(pHodBetaTrk);
		}


		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z)
		{
			H1_hCerNpeSum_pid_acc->Fill(hCerNpeSum);
			H1_hCalEtotNorm_pid_acc->Fill(hCalEtotNorm);
			H1_hCalEtotTrkNorm_pid_acc->Fill(hCalEtotTrkNorm);
			H1_hHodBetaNtrk_pid_acc->Fill(hHodBetaNtrk);
			//H1_hHodBetaTrk_pid_acc->Fill(hHodBetaTrk);
			H1_pNGCerNpeSum_pid_acc->Fill(pNGCerNpeSum);
			//H1_pHGCerNpeSum_pid_acc->Fill(pHGCerNpeSum);
			H1_pCalEtotNorm_pid_acc->Fill(pCalEtotNorm);
			H1_pCalEtotTrkNorm_pid_acc->Fill(pCalEtotTrkNorm);
			H1_pHodBetaNtrk_pid_acc->Fill(pHodBetaNtrk);
			//H1_pHodBetaTrk_pid_acc->Fill(pHodBetaTrk);
		}


		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z && kin)
		{
			H1_hCerNpeSum_pid_acc_kin->Fill(hCerNpeSum);
			H1_hCalEtotNorm_pid_acc_kin->Fill(hCalEtotNorm);
			H1_hCalEtotTrkNorm_pid_acc_kin->Fill(hCalEtotTrkNorm);
			H1_hHodBetaNtrk_pid_acc_kin->Fill(hHodBetaNtrk);
			//H1_hHodBetaTrk_pid_acc_kin->Fill(hHodBetaTrk);
			H1_pNGCerNpeSum_pid_acc_kin->Fill(pNGCerNpeSum);
			//H1_pHGCerNpeSum_pid_acc_kin->Fill(pHGCerNpeSum);
			H1_pCalEtotNorm_pid_acc_kin->Fill(pCalEtotNorm);
			H1_pCalEtotTrkNorm_pid_acc_kin->Fill(pCalEtotTrkNorm);
			H1_pHodBetaNtrk_pid_acc_kin->Fill(pHodBetaNtrk);
			//H1_pHodBetaTrk_pid_acc_kin->Fill(pHodBetaTrk);
		}
			
		/*
		H1_ep_ctime->Fill(ep_ctime - ctime_offset);

		if (pid_hms && pid_shms)
		{
			H1_ep_ctime_pid->Fill(ep_ctime - ctime_offset);
		}

		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z)
		{
			H1_ep_ctime_pid_acc->Fill(ep_ctime - ctime_offset);
		}

		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z && kin)
		{
			H1_ep_ctime_pid_acc_kin->Fill(ep_ctime - ctime_offset);
		}
		*/

		//--------------------------------------
		//Acceptance Histogram Bins
		//--------------------------------------
		//----- Hadron Arm Focal Plane -----
		H1_hxfp->Fill(hxfp);
		H1_hxpfp->Fill(hxpfp / dtr);
		H1_hyfp->Fill(hyfp);
		H1_hypfp->Fill(hypfp / dtr);
		//----- Hadron Arm Reconstructed Quantities -----
		H1_hytar->Fill(hytar);
		H1_hyptar->Fill(hyptar / dtr);
		H1_hxptar->Fill(hxptar / dtr);
		H1_hdelta->Fill(hdelta);
		//----- Target Reconstruction (Hall Coord. System) -----
		H1_htarx->Fill(htarx);
		H1_htary->Fill(htary);
		H1_htarz->Fill(htarz);
		//----- HMS Collimator -----
		H1_hXColl->Fill(hXColl);
		H1_hYColl->Fill(hYColl);
		//----- Hadron Arm Focal Plane -----
		H1_exfp->Fill(exfp);
		H1_expfp->Fill(expfp / dtr);
		H1_eyfp->Fill(eyfp);
		H1_eypfp->Fill(eypfp / dtr);
		//----- Hadron Arm Reconstructed Quantities -----
		H1_eytar->Fill(eytar);
		H1_eyptar->Fill(eyptar / dtr);
		H1_exptar->Fill(exptar / dtr);
		H1_edelta->Fill(edelta);
		//----- Target Reconstruction (Hall Coord. System) -----
		H1_etarx->Fill(etarx);
		H1_etary->Fill(etary);
		H1_etarz->Fill(etarz);
		//----- SHMS Collimator -----
		H1_eXColl->Fill(eXColl);
		H1_eYColl->Fill(eYColl);
		//----- -----
		H1_ztar_diff->Fill(ztar_diff);
		

		if (pid_hms && pid_shms)
		{
			//----- Hadron Arm Focal Plane -----
			H1_hxfp_pid->Fill(hxfp);
			H1_hxpfp_pid->Fill(hxpfp / dtr);
			H1_hyfp_pid->Fill(hyfp);
			H1_hypfp_pid->Fill(hypfp / dtr);
			//----- Hadron Arm Reconstructed Quantities -----
			H1_hytar_pid->Fill(hytar);
			H1_hyptar_pid->Fill(hyptar / dtr);
			H1_hxptar_pid->Fill(hxptar / dtr);
			H1_hdelta_pid->Fill(hdelta);
			//----- Target Reconstruction (Hall Coord. System) -----
			H1_htarx_pid->Fill(htarx);
			H1_htary_pid->Fill(htary);
			H1_htarz_pid->Fill(htarz);
			//----- HMS Collimator -----
			H1_hXColl_pid->Fill(hXColl);
			H1_hYColl_pid->Fill(hYColl);
			//----- Hadron Arm Focal Plane -----
			H1_exfp_pid->Fill(exfp);
			H1_expfp_pid->Fill(expfp / dtr);
			H1_eyfp_pid->Fill(eyfp);
			H1_eypfp_pid->Fill(eypfp / dtr);
			//----- Hadron Arm Reconstructed Quantities -----
			H1_eytar_pid->Fill(eytar);
			H1_eyptar_pid->Fill(eyptar / dtr);
			H1_exptar_pid->Fill(exptar / dtr);
			H1_edelta_pid->Fill(edelta);
			//----- Target Reconstruction (Hall Coord. System) -----
			H1_etarx_pid->Fill(etarx);
			H1_etary_pid->Fill(etary);
			H1_etarz_pid->Fill(etarz);
			//----- SHMS Collimator -----
			H1_eXColl_pid->Fill(eXColl);
			H1_eYColl_pid->Fill(eYColl);
			//----- -----
			H1_ztar_diff_pid->Fill(ztar_diff);
		}

		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z)
		{
			//----- Hadron Arm Focal Plane -----
			H1_hxfp_pid_acc->Fill(hxfp);
			H1_hxpfp_pid_acc->Fill(hxpfp / dtr);
			H1_hyfp_pid_acc->Fill(hyfp);
			H1_hypfp_pid_acc->Fill(hypfp / dtr);
			//----- Hadron Arm Reconstructed Quantities -----
			H1_hytar_pid_acc->Fill(hytar);
			H1_hyptar_pid_acc->Fill(hyptar / dtr);
			H1_hxptar_pid_acc->Fill(hxptar / dtr);
			H1_hdelta_pid_acc->Fill(hdelta);
			//----- Target Reconstruction (Hall Coord. System) -----
			H1_htarx_pid_acc->Fill(htarx);
			H1_htary_pid_acc->Fill(htary);
			H1_htarz_pid_acc->Fill(htarz);
			//----- HMS Collimator -----
			H1_hXColl_pid_acc->Fill(hXColl);
			H1_hYColl_pid_acc->Fill(hYColl);
			//----- Hadron Arm Focal Plane -----
			H1_exfp_pid_acc->Fill(exfp);
			H1_expfp_pid_acc->Fill(expfp / dtr);
			H1_eyfp_pid_acc->Fill(eyfp);
			H1_eypfp_pid_acc->Fill(eypfp / dtr);
			//----- Hadron Arm Reconstructed Quantities -----
			H1_eytar_pid_acc->Fill(eytar);
			H1_eyptar_pid_acc->Fill(eyptar / dtr);
			H1_exptar_pid_acc->Fill(exptar / dtr);
			H1_edelta_pid_acc->Fill(edelta);
			//----- Target Reconstruction (Hall Coord. System) -----
			H1_etarx_pid_acc->Fill(etarx);
			H1_etary_pid_acc->Fill(etary);
			H1_etarz_pid_acc->Fill(etarz);
			//----- SHMS Collimator -----
			H1_eXColl_pid_acc->Fill(eXColl);
			H1_eYColl_pid_acc->Fill(eYColl);
			//----- -----
			H1_ztar_diff_pid_acc->Fill(ztar_diff);
		}


		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z && kin)
		{
			//----- Hadron Arm Focal Plane -----
			H1_hxfp_pid_acc_kin->Fill(hxfp);
			H1_hxpfp_pid_acc_kin->Fill(hxpfp / dtr);
			H1_hyfp_pid_acc_kin->Fill(hyfp);
			H1_hypfp_pid_acc_kin->Fill(hypfp / dtr);
			//----- Hadron Arm Reconstructed Quantities -----
			H1_hytar_pid_acc_kin->Fill(hytar);
			H1_hyptar_pid_acc_kin->Fill(hyptar / dtr);
			H1_hxptar_pid_acc_kin->Fill(hxptar / dtr);
			H1_hdelta_pid_acc_kin->Fill(hdelta);
			//----- Target Reconstruction (Hall Coord. System) -----
			H1_htarx_pid_acc_kin->Fill(htarx);
			H1_htary_pid_acc_kin->Fill(htary);
			H1_htarz_pid_acc_kin->Fill(htarz);
			//----- HMS Collimator -----
			H1_hXColl_pid_acc_kin->Fill(hXColl);
			H1_hYColl_pid_acc_kin->Fill(hYColl);
			//----- Hadron Arm Focal Plane -----
			H1_exfp_pid_acc_kin->Fill(exfp);
			H1_expfp_pid_acc_kin->Fill(expfp / dtr);
			H1_eyfp_pid_acc_kin->Fill(eyfp);
			H1_eypfp_pid_acc_kin->Fill(eypfp / dtr);
			//----- Hadron Arm Reconstructed Quantities -----
			H1_eytar_pid_acc_kin->Fill(eytar);
			H1_eyptar_pid_acc_kin->Fill(eyptar / dtr);
			H1_exptar_pid_acc_kin->Fill(exptar / dtr);
			H1_edelta_pid_acc_kin->Fill(edelta);
			//----- Target Reconstruction (Hall Coord. System) -----
			H1_etarx_pid_acc_kin->Fill(etarx);
			H1_etary_pid_acc_kin->Fill(etary);
			H1_etarz_pid_acc_kin->Fill(etarz);
			//----- SHMS Collimator -----
			H1_eXColl_pid_acc_kin->Fill(eXColl);
			H1_eYColl_pid_acc_kin->Fill(eYColl);
			//----- -----
			H1_ztar_diff_pid_acc_kin->Fill(ztar_diff);
		}


		//--------------------------------------
		//Primary Electron Kinematics
		//--------------------------------------
		H1_Q2->Fill(Q2);
		H1_W->Fill(W);
		H1_nu->Fill(nu);
		H1_ph_q->Fill(ph_q / dtr);
		H1_q->Fill(q);
		H1_qx->Fill(q_x);
		H1_qy->Fill(q_y);
		H1_qz->Fill(q_z);
		H1_th_e->Fill(th_e / dtr);
		H1_th_q->Fill(th_q / dtr);
		H1_xbj->Fill(x_bj);

		if (pid_hms && pid_shms)
		{
			H1_Q2_pid->Fill(Q2);
			H1_W_pid->Fill(W);
			H1_nu_pid->Fill(nu);
			H1_ph_q_pid->Fill(ph_q / dtr);
			H1_q_pid->Fill(q);
			H1_qx_pid->Fill(q_x);
			H1_qy_pid->Fill(q_y);
			H1_qz_pid->Fill(q_z);
			H1_th_e_pid->Fill(th_e / dtr);
			H1_th_q_pid->Fill(th_q / dtr);
			H1_xbj_pid->Fill(x_bj);
		}


		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z)
		{
			H1_Q2_pid_acc->Fill(Q2);
			H1_W_pid_acc->Fill(W);
			H1_nu_pid_acc->Fill(nu);
			H1_ph_q_pid_acc->Fill(ph_q / dtr);
			H1_q_pid_acc->Fill(q);
			H1_qx_pid_acc->Fill(q_x);
			H1_qy_pid_acc->Fill(q_y);
			H1_qz_pid_acc->Fill(q_z);
			H1_th_e_pid_acc->Fill(th_e / dtr);
			H1_th_q_pid_acc->Fill(th_q / dtr);
			H1_xbj_pid_acc->Fill(x_bj);
		}
		

		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z && kin)
		{
			H1_Q2_pid_acc_kin->Fill(Q2);
			H1_W_pid_acc_kin->Fill(W);
			H1_nu_pid_acc_kin->Fill(nu);
			H1_ph_q_pid_acc_kin->Fill(ph_q / dtr);
			H1_q_pid_acc_kin->Fill(q);
			H1_qx_pid_acc_kin->Fill(q_x);
			H1_qy_pid_acc_kin->Fill(q_y);
			H1_qz_pid_acc_kin->Fill(q_z);
			H1_th_e_pid_acc_kin->Fill(th_e / dtr);
			H1_th_q_pid_acc_kin->Fill(th_q / dtr);
			H1_xbj_pid_acc_kin->Fill(x_bj);

			if (numLoops % 50 == 0)
			{
				cout << "Q2:" << Q2 << endl;
				cout << "W:" << W << endl;
				cout << "nu:" << nu << endl;
				cout << "phq:" << ph_q / dtr << endl;
				cout << "q3m:" << q << endl;
				cout << "qx:" << q_x << endl;
				cout << "qy:" << q_y << endl;
				cout << "qz:" << q_z << endl;
				cout << "theta_e:" << th_e / dtr << endl;
				cout << "th_q:" << th_q / dtr << endl;
				cout << "xbj:" << x_bj << endl;
			}
		}


		//--------------------------------------
		//Secondary Hadron Kinematics
		//--------------------------------------
		H1_prec_x->Fill(prec_x);
		H1_prec_y->Fill(prec_y);
		H1_prec_z->Fill(prec_z);
		//H1_emiss->Fill(emiss);
		H1_emiss_nuc->Fill(emiss_nuc);
		H1_ph_rq->Fill(ph_rq / dtr);
		H1_ph_pq->Fill(ph_pq / dtr);
		H1_pmiss->Fill(pmiss);
		H1_pmiss_x->Fill(pmiss_x);
		H1_pmiss_y->Fill(pmiss_y);
		H1_pmiss_z->Fill(pmiss_z);
		H1_Tr->Fill(Tr);
		H1_th_rq->Fill(cth_rq / dtr);
		H1_cth_rq->Fill(cos(cth_rq));
		H1_th_pq->Fill(th_pq / dtr);
		H1_Tx->Fill(Tx);
		H1_xangle->Fill(xangle / dtr);
		H1_mmiss->Fill(mmiss);
		//cout << mmiss << endl;////////////////////
		H1_th_p->Fill(th_p / dtr);
		//H1_omega->Fill(omega);

		if (pid_hms && pid_shms)
		{
			H1_prec_x_pid->Fill(prec_x);
			H1_prec_y_pid->Fill(prec_y);
			H1_prec_z_pid->Fill(prec_z);
			//H1_emiss_pid->Fill(emiss);
			H1_emiss_nuc_pid->Fill(emiss_nuc);
			H1_ph_rq_pid->Fill(ph_rq / dtr);
			H1_ph_pq_pid->Fill(ph_pq / dtr);
			H1_pmiss_pid->Fill(pmiss);
			H1_pmiss_x_pid->Fill(pmiss_x);
			H1_pmiss_y_pid->Fill(pmiss_y);
			H1_pmiss_z_pid->Fill(pmiss_z);
			H1_Tr_pid->Fill(Tr);
			H1_th_rq_pid->Fill(cth_rq / dtr);
			H1_cth_rq_pid->Fill(cos(cth_rq));
			H1_th_pq_pid->Fill(th_pq / dtr);
			H1_Tx_pid->Fill(Tx);
			H1_xangle_pid->Fill(xangle / dtr);
			H1_mmiss_pid->Fill(mmiss);
			H1_th_p_pid->Fill(th_p / dtr);
			//H1_omega_pid->Fill(omega);
		}


		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z)
		{
			H1_prec_x_pid_acc->Fill(prec_x);
			H1_prec_y_pid_acc->Fill(prec_y);
			H1_prec_z_pid_acc->Fill(prec_z);
			//H1_emiss_pid_acc->Fill(emiss);
			H1_emiss_nuc_pid_acc->Fill(emiss_nuc);
			H1_ph_rq_pid_acc->Fill(ph_rq / dtr);
			H1_ph_pq_pid_acc->Fill(ph_pq / dtr);
			H1_pmiss_pid_acc->Fill(pmiss);
			H1_pmiss_x_pid_acc->Fill(pmiss_x);
			H1_pmiss_y_pid_acc->Fill(pmiss_y);
			H1_pmiss_z_pid_acc->Fill(pmiss_z);
			H1_Tr_pid_acc->Fill(Tr);
			H1_th_rq_pid_acc->Fill(cth_rq / dtr);
			H1_cth_rq_pid_acc->Fill(cos(cth_rq));
			H1_th_pq_pid_acc->Fill(th_pq / dtr);
			H1_Tx_pid_acc->Fill(Tx);
			H1_xangle_pid_acc->Fill(xangle / dtr);
			H1_mmiss_pid_acc->Fill(mmiss);
			H1_th_p_pid_acc->Fill(th_p / dtr);
			//H1_omega_pid_acc->Fill(omega);
		}


		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z && kin)
		{
			H1_prec_x_pid_acc_kin->Fill(prec_x);
			H1_prec_y_pid_acc_kin->Fill(prec_y);
			H1_prec_z_pid_acc_kin->Fill(prec_z);
			//H1_emiss_pid_acc_kin->Fill(emiss);
			H1_emiss_nuc_pid_acc_kin->Fill(emiss_nuc);
			H1_ph_rq_pid_acc_kin->Fill(ph_rq / dtr);
			H1_ph_pq_pid_acc_kin->Fill(ph_pq / dtr);
			H1_pmiss_pid_acc_kin->Fill(pmiss);
			H1_pmiss_x_pid_acc_kin->Fill(pmiss_x);
			H1_pmiss_y_pid_acc_kin->Fill(pmiss_y);
			H1_pmiss_z_pid_acc_kin->Fill(pmiss_z);
			H1_Tr_pid_acc_kin->Fill(Tr);
			H1_th_rq_pid_acc_kin->Fill(cth_rq / dtr);
			H1_cth_rq_pid_acc_kin->Fill(cos(cth_rq));
			H1_th_pq_pid_acc_kin->Fill(th_pq / dtr);
			H1_Tx_pid_acc_kin->Fill(Tx);
			H1_xangle_pid_acc_kin->Fill(xangle / dtr);
			H1_mmiss_pid_acc_kin->Fill(mmiss);
			H1_th_p_pid_acc_kin->Fill(th_p / dtr);
			//H1_omega_pid_acc_kin->Fill(omega);

			if (numLoops % 50 == 0)
			{
				cout << "emiss_nuc:" << emiss_nuc << endl;
				cout << "ph_rq:" << ph_rq / dtr << endl;
				cout << "ph_pq:" << ph_pq / dtr << endl;
				cout << "pmiss:" << pmiss << endl;
				cout << "pmiss_x:" << pmiss_x << endl;
				cout << "pmiss_y:" << pmiss_y << endl;
				cout << "pmiss_z:" << pmiss_z << endl;
				cout << "th_rq:" << cth_rq / dtr << endl;
				cout << "th_pq:" << th_pq / dtr << endl;
				cout << "mmisss:" << mmiss << endl;
				cout << "th_p:" << th_p / dtr << endl;
			}
		}
		

		//--------------------------------------
		//2D Histograms
		//--------------------------------------
		H2_hxfp_hyfp->Fill(hyfp, hxfp);
		H2_exfp_eyfp->Fill(eyfp, exfp);
		H2_hxptar_exptar->Fill(exptar / dtr, hxptar / dtr);
		H2_hyptar_eyptar->Fill(eyptar / dtr, hyptar / dtr);
		H2_hdelta_edelta->Fill(edelta, hdelta);
		H2_hXColl_hYColl->Fill(hYColl, hXColl);
		H2_eXColl_eYColl->Fill(eYColl, eXColl);
		H2_xbj_Q2->Fill(Q2, x_bj);
		//H2_the_kf->Fill(hyfp, th_e / dtr);
		H2_Em_Pm->Fill(pmiss, emiss_nuc);
		H2_thrq_pmiss->Fill(pmiss, cth_rq /dtr);
		H2_cthrq_pmiss->Fill(pmiss, cos(cth_rq));
		//H2_thp_Pf->Fill(hyfp, th_p / dtr);
		H2_W_pmiss->Fill(pmiss, W);

		if (pid_hms && pid_shms)
		{
			H2_hxfp_hyfp_pid->Fill(hyfp, hxfp);
			H2_exfp_eyfp_pid->Fill(eyfp, exfp);
			H2_hxptar_exptar_pid->Fill(exptar / dtr, hxptar / dtr);
			H2_hyptar_eyptar_pid->Fill(eyptar / dtr, hyptar / dtr);
			H2_hdelta_edelta_pid->Fill(edelta, hdelta);
			H2_hXColl_hYColl_pid->Fill(hYColl, hXColl);
			H2_eXColl_eYColl_pid->Fill(eYColl, eXColl);
			H2_xbj_Q2_pid->Fill(Q2, x_bj);
			//H2_the_kf_pid->Fill(hyfp, th_e / dtr);
			H2_Em_Pm_pid->Fill(pmiss, emiss_nuc);
			H2_thrq_pmiss_pid->Fill(pmiss, cth_rq / dtr);
			H2_cthrq_pmiss_pid->Fill(pmiss, cos(cth_rq));
			//H2_thp_Pf_pid->Fill(hyfp, th_p / dtr);
			H2_W_pmiss_pid->Fill(pmiss, W);
		}

		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z)
		{
			H2_hxfp_hyfp_pid_acc->Fill(hyfp, hxfp);
			H2_exfp_eyfp_pid_acc->Fill(eyfp, exfp);
			H2_hxptar_exptar_pid_acc->Fill(exptar / dtr, hxptar / dtr);
			H2_hyptar_eyptar_pid_acc->Fill(eyptar / dtr, hyptar / dtr);
			H2_hdelta_edelta_pid_acc->Fill(edelta, hdelta);
			H2_hXColl_hYColl_pid_acc->Fill(hYColl, hXColl);
			H2_eXColl_eYColl_pid_acc->Fill(eYColl, eXColl);
			H2_xbj_Q2_pid_acc->Fill(Q2, x_bj);
			//H2_the_kf_pid_acc->Fill(hyfp, th_e / dtr);
			H2_Em_Pm_pid_acc->Fill(pmiss, emiss_nuc);
			H2_thrq_pmiss_pid_acc->Fill(pmiss, cth_rq / dtr);
			H2_cthrq_pmiss_pid_acc->Fill(pmiss, cos(cth_rq));
			//H2_thp_Pf_pid_acc->Fill(hyfp, th_p / dtr);
			H2_W_pmiss_pid_acc->Fill(pmiss, W);
		}
		
		if (pid_hms && pid_shms && accp_hms && accp_shms && accp_z && kin)
		{
			H2_hxfp_hyfp_pid_acc_kin->Fill(hyfp, hxfp);
			H2_exfp_eyfp_pid_acc_kin->Fill(eyfp, exfp);
			H2_hxptar_exptar_pid_acc_kin->Fill(exptar / dtr, hxptar / dtr);
			H2_hyptar_eyptar_pid_acc_kin->Fill(eyptar / dtr, hyptar / dtr);
			H2_hdelta_edelta_pid_acc_kin->Fill(edelta, hdelta);
			H2_hXColl_hYColl_pid_acc_kin->Fill(hYColl, hXColl);
			H2_eXColl_eYColl_pid_acc_kin->Fill(eYColl, eXColl);
			H2_xbj_Q2_pid_acc_kin->Fill(Q2, x_bj);
			//H2_the_kf_pid_acc_kin->Fill(hyfp, th_e / dtr);
			H2_Em_Pm_pid_acc_kin->Fill(pmiss, emiss_nuc);
			H2_thrq_pmiss_pid_acc_kin->Fill(pmiss, cth_rq / dtr);
			H2_cthrq_pmiss_pid_acc_kin->Fill(pmiss, cos(cth_rq));
			//H2_thp_Pf_pid_acc_kin->Fill(hyfp, th_p / dtr);
			H2_W_pmiss_pid_acc_kin->Fill(pmiss, W);
		}




		if(numLoops > max_loops){ cout << "Exceeded max_loops" << endl;	break;} //leave the loop after completing the loop max_loops times
		numLoops++;
	}//End Loop



	
	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	//-----------------------------------------------------------------------Write Out Histograms-----------------------------------------------------------------------
	//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
	TFile* outROOT;// = new TFile("RunResults.root", "RECREATE");																//Define output
	
	if (user_runtype == "MF")
	{
		if (user_target == "LD2")
		{
			outROOT = new TFile("MF_D2_Results.root", "RECREATE");
		}
		else if (user_target == "Be9")
		{
			outROOT = new TFile("MF_Be9_Results.root", "RECREATE");
		}
		else if (user_target == "B10")
		{
			outROOT = new TFile("MF_B10_Results.root", "RECREATE");
		}
		else if (user_target == "B11")
		{
			outROOT = new TFile("MF_B11_Results.root", "RECREATE");
		}
		else if (user_target == "C12")
		{
			outROOT = new TFile("MF_C12_Results.root", "RECREATE");
		}
		else if (user_target == "Ca40")
		{
			outROOT = new TFile("MF_Ca40_Results.root", "RECREATE");
		}
		else if (user_target == "Ca48")
		{
			outROOT = new TFile("MF_Ca48_Results.root", "RECREATE");
		}
		else if (user_target == "Fe54")
		{
			outROOT = new TFile("MF_Fe54_Results.root", "RECREATE");
		}
	}
	else if (user_runtype == "SRC")
	{
		if (user_target == "LD2")
		{
			outROOT = new TFile("SRC_D2_Results.root", "RECREATE");
		}
		else if (user_target == "Be9")
		{
			outROOT = new TFile("SRC_Be9_Results.root", "RECREATE");
		}
		else if (user_target == "B10")
		{
			outROOT = new TFile("SRC_B10_Results.root", "RECREATE");
		}
		else if (user_target == "B11")
		{
			outROOT = new TFile("SRC_B11_Results.root", "RECREATE");
		}
		else if (user_target == "C12")
		{
			outROOT = new TFile("SRC_C12_Results.root", "RECREATE");
		}
		else if (user_target == "Ca40")
		{
			outROOT = new TFile("SRC_Ca40_Results.root", "RECREATE");
		}
		else if (user_target == "Ca48")
		{
			outROOT = new TFile("SRC_Ca48_Results.root", "RECREATE");
		}
		else if (user_target == "Fe54")
		{
			outROOT = new TFile("SRC_Fe54_Results.root", "RECREATE");
		}
		else
		{
			cout << "Couldn't find file." << endl;
		}
	}
	
	outROOT->mkdir("PID");										//Make directories to store histograms based on Kinematic
	outROOT->cd("PID");											//Write Kinematics histos to kin_plots directory
	PIDList->Write();

	outROOT->mkdir("Prim_Kin");									//Make directories to store histograms based on Kinematic
	outROOT->cd("Prim_Kin");									//Write Kinematics histos to kin_plots directory
	PKList->Write();

	outROOT->mkdir("Sec_Kin");									//Make directories to store histograms based on Kinematic
	outROOT->cd("Sec_Kin");										//Write Kinematics histos to kin_plots directory
	SKList->Write();

	outROOT->mkdir("HMS_Accp");									//Make directories to store histograms based on Kinematic
	outROOT->cd("HMS_Accp");									//Write Kinematics histos to kin_plots directory
	HAccList->Write();

	outROOT->mkdir("SHMS_Accp");								//Make directories to store histograms based on Kinematic
	outROOT->cd("SHMS_Accp");									//Write Kinematics histos to kin_plots directory
	SAccList->Write();
	
	outROOT->mkdir("2D_Hist");									//Make directories to store histograms based on Kinematic
	outROOT->cd("2D_Hist");										//Write Kinematics histos to kin_plots directory
	HList2->Write();
	
	
	//--------------------------------------
	//PID Histograms Bins
	//--------------------------------------
	//----- HMS -----
	double y2f = 0;
	double y1f = 0;

	be_and_af(H1_hCerNpeSum, H1_hCerNpeSum_pid, H1_hCerNpeSum_pid_acc, H1_hCerNpeSum_pid_acc_kin, 0, 0, 0, y1f, 0.5, 0, 0.5, y2f, "hCerNpeSum.png", "hCerNpeSum", true, outROOT);
	be_and_af(H1_hCalEtotNorm, H1_hCalEtotNorm_pid, H1_hCalEtotNorm_pid_acc, H1_hCalEtotNorm_pid_acc_kin, 0, 0, 0, y1f, 0.6, 0, 0.6, y2f, "hCalEtotNorm.png", "hCalEtotNorm", true, outROOT);
	be_and_af(H1_hCalEtotTrkNorm, H1_hCalEtotTrkNorm_pid, H1_hCalEtotTrkNorm_pid_acc, H1_hCalEtotTrkNorm_pid_acc_kin, 0, 0, 0, y1f, 0.6, 0, 0.6, y2f, "hCalEtotTrkNorm.png", "hCalEtotTrkNorm", true, outROOT);
	be_and_af(H1_hHodBetaNtrk, H1_hHodBetaNtrk_pid, H1_hHodBetaNtrk_pid_acc, H1_hHodBetaNtrk_pid_acc_kin, 0.5, 0, 0.5, y1f, 1.5, 0, 1.5, y2f, "hHodBetaNtrk.png", "hHodBetaNtrk", true, outROOT);
	//----- SHMS -----
	be_and_af(H1_pNGCerNpeSum, H1_pNGCerNpeSum_pid, H1_pNGCerNpeSum_pid_acc, H1_pNGCerNpeSum_pid_acc_kin, 1, 0, 1, y1f, 100, 0, 100, y2f, "pNGCerNpeSum.png", "pNGCerNpeSum", true, outROOT);
	//be_and_af(H1_pHGCerNpeSum, H1_pHGCerNpeSum_pid, H1_pHGCerNpeSum_pid_acc, H1_pHGCerNpeSum_pid_acc_kin, 0, 0, 0, y1f, 1.5, 0, 1.5, y2f, "pHGCerNpeSum.png", "pHGCerNpeSum", true, outROOT);
	be_and_af(H1_pCalEtotNorm, H1_pCalEtotNorm_pid, H1_pCalEtotNorm_pid_acc, H1_pCalEtotNorm_pid_acc_kin, 0.85, 0, 0.85, y1f, 1.25, 0, 1.25, y2f, "pCalEtotNorm.png", "pCalEtotNorm", true, outROOT);
	be_and_af(H1_pCalEtotTrkNorm, H1_pCalEtotTrkNorm_pid, H1_pCalEtotTrkNorm_pid_acc, H1_pCalEtotTrkNorm_pid_acc_kin, 0.85, 0, 0.85, y1f, 1.25, 0, 1.25, y2f, "pCalEtotTrkNorm.png", "pCalEtotTrkNorm", true, outROOT);
	be_and_af(H1_pHodBetaNtrk, H1_pHodBetaNtrk_pid, H1_pHodBetaNtrk_pid_acc, H1_pHodBetaNtrk_pid_acc_kin, 0.5, 0, 0.5, y1f, 1.5, 0, 1.5, y2f, "pHodBetaNtrk.png", "pHodBetaNtrk", true, outROOT);
	//----- -----
	//be_and_af(H1_ep_ctime, H1_ep_ctime_pid, H1_ep_ctime_pid_acc, H1_ep_ctime_pid_acc_kin, -2.5, 0, -2.5, y1f, 2.5, 0, 2.5, y2f, "ep_ctime.png", "ep_ctime", true, outROOT);



	Double_t hdelta_min = Cuts::hdelta(0), hdelta_max = Cuts::hdelta(1);
	Double_t hxptar_min = Cuts::hxptar(0), hxptar_max = Cuts::hxptar(1);
	Double_t hyptar_min = Cuts::hyptar(0), hyptar_max = Cuts::hyptar(1);
	Double_t edelta_min = Cuts::edelta(0), edelta_max = Cuts::edelta(1);
	Double_t exptar_min = Cuts::exptar(0), exptar_max = Cuts::exptar(1);
	Double_t eyptar_min = Cuts::eyptar(0), eyptar_max = Cuts::eyptar(1);
	Double_t ztarDiff_min = Cuts::ztarDiff(0), ztarDiff_max = Cuts::ztarDiff(1);
	//----- Hadron Arm Focal Plane -----
	be_and_af(H1_hxfp, H1_hxfp_pid, H1_hxfp_pid_acc, H1_hxfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hxfp.png", "hxfp", false, outROOT);
	be_and_af(H1_hxpfp, H1_hxpfp_pid, H1_hxpfp_pid_acc, H1_hxpfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hxpfp.png", "hxpfp", false, outROOT);
	be_and_af(H1_hyfp, H1_hyfp_pid, H1_hyfp_pid_acc, H1_hyfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hyfp.png", "hyfp", false, outROOT);
	be_and_af(H1_hypfp, H1_hypfp_pid, H1_hypfp_pid_acc, H1_hypfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hypfp.png", "hypfp", false, outROOT);
	//----- Hadron Arm Reconstructed Quantities -----
	be_and_af(H1_hytar, H1_hytar_pid, H1_hytar_pid_acc, H1_hytar_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hytar.png", "hytar", false, outROOT);
	be_and_af(H1_hyptar, H1_hyptar_pid, H1_hyptar_pid_acc, H1_hyptar_pid_acc_kin, hyptar_min, 0, hyptar_min, y1f, hyptar_max, 0, hyptar_max, y2f, "hyptar.png", "hyptar", false, outROOT);
	be_and_af(H1_hxptar, H1_hxptar_pid, H1_hxptar_pid_acc, H1_hxptar_pid_acc_kin, hxptar_min, 0, hxptar_min, y1f, hxptar_max, 0, hxptar_max, y2f, "hxptar.png", "hxptar", false, outROOT);
	be_and_af(H1_hdelta, H1_hdelta_pid, H1_hdelta_pid_acc, H1_hdelta_pid_acc_kin, hdelta_min, 0, hdelta_min, y1f, hdelta_max, 0, hdelta_max, y2f, "hdelta.png", "hdelta", true, outROOT);
	//----- Target Reconstruction (Hall Coord. System) -----
	be_and_af(H1_htarx, H1_htarx_pid, H1_htarx_pid_acc, H1_htarx_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htarx.png", "htarx", false, outROOT);
	be_and_af(H1_htary, H1_htary_pid, H1_htary_pid_acc, H1_htary_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htary.png", "htary", false, outROOT);
	be_and_af(H1_htarz, H1_htarz_pid, H1_htarz_pid_acc, H1_htarz_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "htarz.png", "htarz", false, outROOT);
	//----- HMS Collimator -----
	be_and_af(H1_hXColl, H1_hXColl_pid, H1_hXColl_pid_acc, H1_hXColl_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hXColl.png", "hXColl", false, outROOT);
	be_and_af(H1_hYColl, H1_hYColl_pid, H1_hYColl_pid_acc, H1_hYColl_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "hYColl.png", "hYColl", false, outROOT);

	//----- Hadron Arm Focal Plane -----
	be_and_af(H1_exfp, H1_exfp_pid, H1_exfp_pid_acc, H1_exfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "exfp.png", "exfp", false, outROOT);
	be_and_af(H1_expfp, H1_expfp_pid, H1_expfp_pid_acc, H1_expfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "expfp.png", "expfp", false, outROOT);
	be_and_af(H1_eyfp, H1_eyfp_pid, H1_eyfp_pid_acc, H1_eyfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eyfp.png", "eyfp", false, outROOT);
	be_and_af(H1_eypfp, H1_eypfp_pid, H1_eypfp_pid_acc, H1_eypfp_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eypfp.png", "eypfp", false, outROOT);
	//----- Hadron Arm Reconstructed Quantities -----
	be_and_af(H1_eytar, H1_eytar_pid, H1_eytar_pid_acc, H1_eytar_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eytar.png", "eytar", false, outROOT);
	be_and_af(H1_eyptar, H1_eyptar_pid, H1_eyptar_pid_acc, H1_eyptar_pid_acc_kin, eyptar_min, 0, eyptar_min, y1f, eyptar_max, 0, eyptar_max, y2f, "eyptar.png", "eyptar", false, outROOT);
	be_and_af(H1_exptar, H1_exptar_pid, H1_exptar_pid_acc, H1_exptar_pid_acc_kin, exptar_min, 0, exptar_min, y1f, exptar_max, 0, exptar_max, y2f, "exptar.png", "exptar", false, outROOT);
	be_and_af(H1_edelta, H1_edelta_pid, H1_edelta_pid_acc, H1_edelta_pid_acc_kin, edelta_min, 0, edelta_min, y1f, edelta_max, 0, edelta_max, y2f, "edelta.png", "edelta", true, outROOT);
	//----- Target Reconstruction (Hall Coord. System) -----
	be_and_af(H1_etarx, H1_etarx_pid, H1_etarx_pid_acc, H1_etarx_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etarx.png", "etarx", false, outROOT);
	be_and_af(H1_etary, H1_etary_pid, H1_etary_pid_acc, H1_etary_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etary.png", "etary", false, outROOT);
	be_and_af(H1_etarz, H1_etarz_pid, H1_etarz_pid_acc, H1_etarz_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "etarz.png", "etarz", false, outROOT);
	//----- SHMS Collimator -----
	be_and_af(H1_eXColl, H1_eXColl_pid, H1_eXColl_pid_acc, H1_eXColl_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eXColl.png", "eXColl", false, outROOT);
	be_and_af(H1_eYColl, H1_eYColl_pid, H1_eYColl_pid_acc, H1_eYColl_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "eYColl.png", "eYColl", false, outROOT);

	//----- -----
	be_and_af(H1_ztar_diff, H1_ztar_diff_pid, H1_ztar_diff_pid_acc, H1_ztar_diff_pid_acc_kin, ztarDiff_min, 0, ztarDiff_min, y1f, ztarDiff_max, 0, ztarDiff_max, y2f, "ztar_diff.png", "ztar_diff", true, outROOT);



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

	
	be_and_af(H1_Q2, H1_Q2_pid, H1_Q2_pid_acc, H1_Q2_pid_acc_kin, Q2_min, 0, Q2_min, y1f, Q2_max, 0, Q2_max, y2f, "Q2.png", "Q2", true, outROOT);
	be_and_af(H1_xbj, H1_xbj_pid, H1_xbj_pid_acc, H1_xbj_pid_acc_kin, xbj_min, 0, xbj_min, y1f, xbj_max, 0, xbj_max, y2f, "xbj.png", "xbj", true, outROOT);
	be_and_af(H1_W, H1_W_pid, H1_W_pid_acc, H1_W_pid_acc_kin, W_min, 0, W_min, y1f, W_max, 0, W_max, y2f, "W.png", "W", true, outROOT);
	be_and_af(H1_nu, H1_nu_pid, H1_nu_pid_acc, H1_nu_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "nu.png", "nu", false, outROOT);
	be_and_af(H1_ph_q, H1_ph_q_pid, H1_ph_q_pid_acc, H1_ph_q_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_q.png", "ph_q", false, outROOT);
	be_and_af(H1_q, H1_q_pid, H1_q_pid_acc, H1_q_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "q.png", "q", false, outROOT);
	be_and_af(H1_qx, H1_qx_pid, H1_qx_pid_acc, H1_qx_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "qx.png", "qx", false, outROOT);
	be_and_af(H1_qy, H1_qy_pid, H1_qy_pid_acc, H1_qy_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "qy.png", "qy", false, outROOT);
	be_and_af(H1_qz, H1_qz_pid, H1_qz_pid_acc, H1_qz_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "qz.png", "qz", false, outROOT);
	be_and_af(H1_th_e, H1_th_e_pid, H1_th_e_pid_acc, H1_th_e_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "the.png", "the", false, outROOT);
	be_and_af(H1_th_q, H1_th_q_pid, H1_th_q_pid_acc, H1_th_q_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "thq.png", "thq", false, outROOT);

	

	//--------------------------------------
	//Secondary Hadron Kinematics
	//--------------------------------------
	be_and_af(H1_pmiss, H1_pmiss_pid, H1_pmiss_pid_acc, H1_pmiss_pid_acc_kin, Pm_min, 0, Pm_min, y1f, Pm_max, 0, Pm_max, y2f, "pmiss.png", "pmiss", true, outROOT);
	be_and_af(H1_th_rq, H1_th_rq_pid, H1_th_rq_pid_acc, H1_th_rq_pid_acc_kin, Src_thrq_min, 0, Src_thrq_min, y1f, Src_thrq_max, 0, Src_thrq_max, y2f, "th_rq.png", "th_rq", true, outROOT);
	be_and_af(H1_cth_rq, H1_cth_rq_pid, H1_cth_rq_pid_acc, H1_cth_rq_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "cth_rq.png", "cth_rq", false, outROOT);//true
	//be_and_af(H1_emiss, H1_emiss_pid, H1_emiss_pid_acc, H1_emiss_pid_acc_kin, Em_min, 0, Em_min, y1f, Em_max, 0, Em_max, y2f, "emiss.png", "emiss", true, outROOT);
	be_and_af(H1_emiss_nuc, H1_emiss_nuc_pid, H1_emiss_nuc_pid_acc, H1_emiss_nuc_pid_acc_kin, Em_min, 0, Em_min, y1f, Em_max, 0, Em_max, y2f, "emiss_nuc.png", "emiss_nuc", true, outROOT);
	be_and_af(H1_prec_x, H1_prec_x_pid, H1_prec_x_pid_acc, H1_prec_x_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Prec_x.png", "Prec_x", false, outROOT);
	be_and_af(H1_prec_y, H1_prec_y_pid, H1_prec_y_pid_acc, H1_prec_y_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Prec_y.png", "Prec_y", false, outROOT);
	be_and_af(H1_prec_z, H1_prec_z_pid, H1_prec_z_pid_acc, H1_prec_z_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Prec_z.png", "Prec_z", false, outROOT);
	be_and_af(H1_ph_rq, H1_ph_rq_pid, H1_ph_rq_pid_acc, H1_ph_rq_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_rq.png", "ph_rq", false, outROOT);
	be_and_af(H1_ph_pq, H1_ph_pq_pid, H1_ph_pq_pid_acc, H1_ph_pq_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "ph_pq.png", "ph_pq", false, outROOT);
	be_and_af(H1_pmiss_x, H1_pmiss_x_pid, H1_pmiss_x_pid_acc, H1_pmiss_x_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_x.png", "pmiss_x", false, outROOT);
	be_and_af(H1_pmiss_y, H1_pmiss_y_pid, H1_pmiss_y_pid_acc, H1_pmiss_y_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_y.png", "pmiss_y", false, outROOT);
	be_and_af(H1_pmiss_z, H1_pmiss_z_pid, H1_pmiss_z_pid_acc, H1_pmiss_z_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "pmiss_z.png", "pmiss_z", false, outROOT);
	be_and_af(H1_Tr, H1_Tr_pid, H1_Tr_pid_acc, H1_Tr_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Tr.png", "Tr", false, outROOT);
	be_and_af(H1_th_pq, H1_th_pq_pid, H1_th_pq_pid_acc, H1_th_pq_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "th_pq.png", "th_pq", false, outROOT);
	be_and_af(H1_xangle, H1_xangle_pid, H1_xangle_pid_acc, H1_xangle_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "xangle.png", "xangle", false, outROOT);
	be_and_af(H1_Tx, H1_Tx_pid, H1_Tx_pid_acc, H1_Tx_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "Tx.png", "Tx", false, outROOT);
	be_and_af(H1_mmiss, H1_mmiss_pid, H1_mmiss_pid_acc, H1_mmiss_pid_acc_kin, MM_min, 0, MM_min, y1f, MM_max, 0, MM_max, y2f, "mmiss.png", "mmiss", true, outROOT);
	be_and_af(H1_th_p, H1_th_p_pid, H1_th_p_pid_acc, H1_th_p_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "th_p.png", "th_p", false, outROOT);
	//be_and_af(H1_omega, H1_omega_pid, H1_omega_pid_acc, H1_omega_pid_acc_kin, 0, 0, 0, 0, 0, 0, 0, 0, "omega.png", "omega", false, outROOT);

	be_and_af2d(H2_hxfp_hyfp, H2_hxfp_hyfp_pid, H2_hxfp_hyfp_pid_acc, H2_hxfp_hyfp_pid_acc_kin, "H2_hxfp_hyfp.png", "H2_hxfp_hyfp", outROOT);
	be_and_af2d(H2_exfp_eyfp, H2_exfp_eyfp_pid, H2_exfp_eyfp_pid_acc, H2_exfp_eyfp_pid_acc_kin, "H2_exfp_eyfp.png", "H2_exfp_eyfp", outROOT);
	be_and_af2d(H2_hxptar_exptar, H2_hxptar_exptar_pid, H2_hxptar_exptar_pid_acc, H2_hxptar_exptar_pid_acc_kin, "H2_hxptar_exptar.png", "H2_hxptar_exptar", outROOT);
	be_and_af2d(H2_hyptar_eyptar, H2_hyptar_eyptar_pid, H2_hyptar_eyptar_pid_acc, H2_hyptar_eyptar_pid_acc_kin, "H2_hyptar_eyptar.png", "H2_hyptar_eyptar", outROOT);
	be_and_af2d(H2_hdelta_edelta, H2_hdelta_edelta_pid, H2_hdelta_edelta_pid_acc, H2_hdelta_edelta_pid_acc_kin, "H2_hdelta_edelta.png", "H2_hdelta_edelta", outROOT);
	be_and_af2d(H2_hXColl_hYColl, H2_hXColl_hYColl_pid, H2_hXColl_hYColl_pid_acc, H2_hXColl_hYColl_pid_acc_kin, "H2_hXColl_hYColl.png", "H2_hXColl_hYColl", outROOT);
	be_and_af2d(H2_eXColl_eYColl, H2_eXColl_eYColl_pid, H2_eXColl_eYColl_pid_acc, H2_eXColl_eYColl_pid_acc_kin, "H2_eXColl_eYColl.png", "H2_eXColl_eYColl", outROOT);
	be_and_af2d(H2_xbj_Q2, H2_xbj_Q2_pid, H2_xbj_Q2_pid_acc, H2_xbj_Q2_pid_acc_kin, "H2_xbj_Q2.png", "H2_xbj_Q2", outROOT);
	//be_and_af2d(H2_the_kf, H2_the_kf_pid, H2_the_kf_pid_acc, H2_the_kf_pid_acc_kin, "H2_the_kf.png", "H2_the_kf", outROOT);
	be_and_af2d(H2_Em_Pm, H2_Em_Pm_pid, H2_Em_Pm_pid_acc, H2_Em_Pm_pid_acc_kin, "H2_Em_Pm.png", "H2_Em_Pm", outROOT);
	be_and_af2d(H2_thrq_pmiss, H2_thrq_pmiss_pid, H2_thrq_pmiss_pid_acc, H2_thrq_pmiss_pid_acc_kin, "H2_thrq_pmiss.png", "H2_thrq_pmiss", outROOT);
	be_and_af2d(H2_cthrq_pmiss, H2_cthrq_pmiss_pid, H2_cthrq_pmiss_pid_acc, H2_cthrq_pmiss_pid_acc_kin, "H2_cthrq_pmiss.png", "H2_cthrq_pmiss", outROOT);
	//be_and_af2d(H2_thp_Pf, H2_thp_Pf_pid, H2_thp_Pf_pid_acc, H2_thp_Pf_pid_acc_kin, "H2_thp_Pf.png", "H2_thp_Pf", outROOT);
	be_and_af2d(H2_W_pmiss, H2_W_pmiss_pid, H2_W_pmiss_pid_acc, H2_W_pmiss_pid_acc_kin, "H2_W_pmiss.png", "H2_W_pmiss", outROOT);



	outROOT->Close();											//Close File
	//good->Print();											//Print reduced cloned tree
	//good.Write();												//Write reduced cloned tree
}
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
//-----------------------------------------------------------------------------End Main-----------------------------------------------------------------------------
//*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=



void plotstuff(TH1F* hist1, const char* name, const char* title, TFile* outRoot)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();
	
	hist1->Draw();

	outRoot->cd(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}



void be_and_af(TH1F* hist1, TH1F* hist2, TH1F* hist3, TH1F* hist4, double x1i, double y1i, double x1f, double y1f, double x2i, double y2i, double x2f, double y2f, const char* name, const char* title, bool line, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->Divide(2,2);

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

	c->cd(3); hist3->Draw("HIST SAME");

	if (line)
	{
		y1f = hist3->GetMaximum();
		y2f = hist3->GetMaximum();
		TLine* line4 = new TLine(x1i, y1i, x1f, y1f); line4->SetLineColor(2); line4->SetLineWidth(2); line4->Draw("SAME");
		TLine* line5 = new TLine(x2i, y2i, x2f, y2f); line5->SetLineColor(2); line5->SetLineWidth(2); line5->Draw("SAME");
	}

	c->cd(4); hist4->Draw("HIST SAME");

	if (line)
	{
		y1f = hist4->GetMaximum();
		y2f = hist4->GetMaximum();
		TLine* line6 = new TLine(x1i, y1i, x1f, y1f); line6->SetLineColor(2); line6->SetLineWidth(2); line6->Draw("SAME");
		TLine* line7 = new TLine(x2i, y2i, x2f, y2f); line7->SetLineColor(2); line7->SetLineWidth(2); line7->Draw("SAME");
	}

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}



void be_and_af2d(TH2F* hist1, TH2F* hist2, TH2F* hist3, TH2F* hist4, const char* name, const char* title, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->Divide(2, 2);

	c->cd(1); hist1->Draw("Col");
	c->cd(2); hist2->Draw("Col SAME");
	c->cd(3); hist3->Draw("Col SAME");
	c->cd(4); hist4->Draw("Col SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}