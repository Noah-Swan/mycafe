//eventual residual fits, width of calibrated one
//shms_dc_anomwire_1u2.png y-axis range
//shms_hodo_beta_vs_xfp y-axis label
//shms_hodo_beta_vs_xfp y-axis label
//everything xfocal plane

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include "TROOT.h"
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
using namespace std;

//-----------------------
// Constants
//------------------------
Double_t pi = 3.141592654;
Double_t dtr = pi/180.;
Double_t MP = 0.938272; //Proton Mass GeV
Double_t MD = 1.87561; //GeV
Double_t MN = 0.939566; //Neutron Mass GeV
Double_t me = 0.000510998; //GeV

void Draw2d(TH2F*, const char*, const char*, TFile*);
void HotWire(TH2F*, const char*, const char*, double, double, double, double, const char*, const char*, TFile*);
void OverLap2(TH1F*, TH1F*, double, double, double, double, const char*, const char*, bool, TFile*);
void Residuals(TH1F*, TH1F*, double, double, double, double, const char*, const char*, bool, TFile*);
void OverLap3(TH1F*, TH1F*, TH1F*, double, double, double, double, const char*, const char*, const char*, bool, TFile*);
void OverLap5(TH1F*, TH1F*, TH1F*, TH1F*, TH1F*, double, double, double, double, const char*, const char*, const char*, bool, TFile*);
void OverLay3(TGraph*, TGraph*, TGraph*, const char*, const char*, const char*, bool, TFile*);
void OverLay5(TGraph*, TGraph*, TGraph*, TGraph*, TGraph*, const char*, const char*, const char*, bool, TFile*);
void Compare2(TH2F*, TH2F*, double, double, double, double, const char*, const char*, bool, TFile*);



void compare_histos_light(TString, TString, TString, TString, TString, TString, TString, TString, TString, TString,
	TString, TString, TString,
	TString, TString, TString, TString, TString,
	bool, TFile*);

/*
void compare_4histos_light(TString, TString, TString, TString, TString, TString, TString, TString,
	TString, TString, TString,
	TString, TString, TString, TString,
	bool, TFile*);
*/

void compare_histos_heavy(TString, TString, TString, TString, TString, TString,
	TString, TString, TString,
	TString, TString, TString,
	bool, TFile*);

/*
void compare_2histos_heavy(TString, TString, TString, TString,
	TString, TString, TString,
	TString, TString,
	bool, TFile*);
*/




//*****************************************************************************************************************************************************************************
//Begin Main
//*****************************************************************************************************************************************************************************
void mkhistos() 
{
	/*
	gStyle->SetOptStat(1);//have stat window
	gStyle->SetOptFit(0011);//fit window
	gStyle->SetOptTitle(1);//have title
	gStyle->SetLabelSize(0.03, "xyz");
	gStyle->SetTitleSize(0.04, "xyz");
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadTopMargin(0.06);
	gStyle->SetPadLeftMargin(0.14);
	gStyle->SetPadRightMargin(0.2);
	gStyle->SetImageScaling(3.);
	*/

	TFile* outROOT = new TFile("Calibration_Results.root", "RECREATE");
	
	TH1F* hist1 = 0;	TH1F* hist2 = 0;	TH1F* hist3 = 0;	TH1F* hist4 = 0;	TH1F* hist5 = 0;
	TH2F* temp21 = 0; TH2F* hist21 = 0; TH2F* hist22 = 0; TH2F* hist23 = 0; TH2F* hist24 = 0; TH2F* hist25 = 0;
	TGraph* gtemp = 0; TGraph* g1 = 0; TGraph* g2 = 0; TGraph* g3 = 0; TGraph* g4 = 0; TGraph* g5 = 0;

	TH1F** h1 = new TH1F*[5]; h1[0] = hist1; h1[1] = hist2; h1[2] = hist3; h1[3] = hist4; h1[4] = hist5;
	TH2F** h2 = new TH2F*[5]; h2[0] = hist21; h2[1] = hist22; h2[2] = hist23; h2[3] = hist24; h2[4] = hist25;

	/*
	TH1F* histm0 = 0; TH1F* histm1 = 0; TH1F* histm2 = 0; TH1F* histm3 = 0;
	TH1F* histm4 = 0; TH1F* histm5 = 0; TH1F* histm6 = 0; TH1F* histm7 = 0;
	TH1F* histm8 = 0; TH1F *histm9 = 0; TH1F* histm10 = 0; TH1F* histm11 = 0;
	TH2F* hist2m0 = 0; TH2F* hist2m1 = 0; TH2F* hist2m2 = 0; TH2F* hist2m3 = 0;
	TH2F* hist2m4 = 0; TH2F* hist2m5 = 0; TH2F* hist2m6 = 0; TH2F* hist2m7 = 0;
	TH2F* hist2m8 = 0; TH2F* hist2m9 = 0; TH2F* hist2m10 = 0; TH2F* hist2m11 = 0;

	TH1F** h1m = new TH1F*[12];
	TH2F** h2m = new TH2F* [12];

	h1m[0] =histm0; h1m[1] =histm1; h1m[2] =histm2; h1m[3] =histm3;
	h1m[4] =histm4; h1m[5] =histm5; h1m[6] =histm6; h1m[7] =histm7;
	h1m[8] =histm8; h1m[9] =histm9; h1m[10] =histm10; h1m[11] = histm11;	
	h2m[0] = hist2m0; h2m[1] = hist2m1; h2m[2] = hist2m2; h2m[3] = hist2m3;
	h2m[4] = hist2m4; h2m[5] = hist2m5; h2m[6] = hist2m6; h2m[7] = hist2m7;
	h2m[8] = hist2m8; h2m[9] = hist2m9; h2m[10] = hist2m10; h2m[11] = hist2m11;
	*/

	//TCanvas* res = new TCanvas("res", "res", 1920, 1080); c->Divide(6, 2);
	//TCanvas* dist = new TCanvas("dist", "dist", 1920, 1080); c->Divide(6, 2);
	//TCanvas* wire1 = new TCanvas("wire1", "wire1", 1920, 1080); c->Divide(6, 2);
	//TCanvas* wire2 = new TCanvas("wire2", "wire2", 1920, 1080); c->Divide(6, 2);

	double yf = 0;

	//----- HMS Histograms: Enabled/Disabled -----
	bool hms_hod_1d = true;
	bool hms_hod_2d = true;
	bool hms_dc_res = false;
	bool hms_dc_dist = false;
	bool hms_dc_dt_v_wn = false;
	bool hms_dc_hotwire = false;
	bool hms_cal_1d = true;
	bool hms_cal_2d = true;
	bool hms_light_mf_1d = true;
	bool hms_light_mf_2d = true;
	bool hms_heavy_mf_1d = true;
	bool hms_heavy_mf_2d = true;
	bool hms_light_src_1d = true;
	bool hms_light_src_2d = true;
	bool hms_heavy_src_1d = true;
	bool hms_heavy_src_2d = true;
	bool unfocused = true;

	//----- SHMS Histograms: Enabled/Disabled -----
	bool shms_hod_1d = true;
	bool shms_hod_2d = true;
	bool shms_dc_res = false;
	bool shms_dc_dist = false;
	bool shms_dc_wn_v_dt = true;
	bool shms_dc_hotwire = true;
	bool shms_cal_1d = true;
	bool shms_cal_2d = true;
	bool shms_light_mf_1d = true;
	bool shms_light_mf_2d = true;
	bool shms_heavy_mf_1d = true;
	bool shms_heavy_mf_2d = true;
	bool shms_light_src_1d = true;
	bool shms_light_src_2d = true;
	bool shms_heavy_src_1d = true;
	bool shms_heavy_src_2d = true;

	//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
	//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= Parameter Files =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
	//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
	
	//--------------------------------------------------------------------------
	//------------------------------ HMS Hodoscope -----------------------------
	//--------------------------------------------------------------------------
	
	//----- Files -----
	//Calibrated
	TFile* hms_replay_4407_hodo_c = new TFile("../../../../OFFLINE/CALIBRATIONS/hodoscopes/hms_replay_production_all_4407_100000_hodCalib.root", "READ"); //Open calibrated unfocused run 4407, electrons incident on hodoscopes
	TTree* thms_replay_4407_hodo_c = (TTree*)hms_replay_4407_hodo_c->Get("T");																		//Open relevant tree

	//Uncalibrated
	TFile* hms_replay_4407_hodo_uc = new TFile("../../hms_replay_production_all_4407_1000000.root", "READ");												//Open uncalibrated unfocused run 4407, electrons incident on hodoscopes
	TTree* thms_replay_4407_hodo_uc = (TTree*)hms_replay_4407_hodo_uc->Get("T");																	//Open relevant tree

	//----- 1d -----
	if (hms_hod_1d)
	{
		//Beta
		hist1 = new TH1F("hist1", "HMS Hodo Beta; Beta (p/E)", 100, 0.1, 1.5);																		//Create Histogram for calibrated HMS hodo beta
		hist2 = new TH1F("hist2", "HMS Hodo Beta; Beta (p/E)", 100, 0.1, 1.5);																		//Create Histogram for uncalibrated HMS hodo beta
		thms_replay_4407_hodo_c->Draw("H.hod.beta>>hist1");
		thms_replay_4407_hodo_uc->Draw("H.hod.beta>>hist2");
		OverLap2(hist1, hist2, 1, 0, 1, yf, "hms_hodo_beta.png", "HMS Hodo Beta", true, outROOT);														//Draw calibrated & uncalibrated HMS hodo beta with a line
		delete hist1; delete hist2;																													//Clear hist1 & hist2 for future use
	}

	//----- 2d -----
	if (hms_hod_2d)
	{
		//Beta vs X_fp
		hist21 = new TH2F("hist21", "HMS Calibrated Hodo; X Focal Plane [cm]; Beta (p/E)", 100, -50, 50, 100, 0.1, 1.5);
		hist22 = new TH2F("hist22", "HMS Uncalibrated Hodo; X Focal Plane [cm]; Beta (p/E)", 100, -50, 50, 100, 0.1, 1.5);
		thms_replay_4407_hodo_c->Draw("H.hod.beta:H.dc.x_fp>>hist21");
		thms_replay_4407_hodo_uc->Draw("H.hod.beta:H.dc.x_fp>>hist22");
		Compare2(hist21, hist22, -50, 1, 50, 1, "hms_hodo_beta_vs_xfp.png", "HMS Hodo Beta vs X_fp", true, outROOT);
		delete hist21; delete hist22;

		//Beta vs Y_fp
		hist21 = new TH2F("hist21", "HMS Calibrated Hodo; Y Focal Plane [cm]; Beta (p/E)", 100, -50, 50, 100, 0.1, 1.5);
		hist22 = new TH2F("hist22", "HMS Uncalibrated Hodo; Y Focal Plane [cm]; Beta (p/E)", 100, -50, 50, 100, 0.1, 1.5);
		thms_replay_4407_hodo_c->Draw("H.hod.beta:H.dc.y_fp>>hist21");
		thms_replay_4407_hodo_uc->Draw("H.hod.beta:H.dc.y_fp>>hist22");
		Compare2(hist21, hist22, -50, 1, 50, 1, "hms_hodo_beta_vs_yfp.png", "HMS Hodo Beta vs Y_fp", true, outROOT);
		delete hist21; delete hist22;
	}
	


	//---------------------------------------------------------------------------
	//------------------------------ SHMS Hodoscope -----------------------------
	//---------------------------------------------------------------------------

	//----- Files -----
	//Calibrated
	TFile* shms_replay_17243_hodo_c = new TFile("../../../../OFFLINE/CALIBRATIONS/hodoscopes/shms_replay_production_all_17243_100000_hodCalib.root", "READ");
	TTree* tshms_replay_17243_hodo_c = (TTree*)shms_replay_17243_hodo_c->Get("T");

	//Uncalibrated
	TFile* shms_replay_17243_hodo_uc = new TFile("../../shms_replay_production_all_17243_1000000.root", "READ");
	TTree* tshms_replay_17243_hodo_uc = (TTree*)shms_replay_17243_hodo_uc->Get("T");

	//----- 1d -----
	if (shms_hod_1d)
	{
		//Beta
		hist1 = new TH1F("hist1", "SHMS Hodo Beta; Beta (p/E)", 100, 0.1, 1.5);
		hist2 = new TH1F("hist2", "SHMS Hodo Beta; Beta (p/E)", 100, 0.1, 1.5);
		tshms_replay_17243_hodo_c->Draw("P.hod.beta>>hist1");
		tshms_replay_17243_hodo_uc->Draw("P.hod.beta>>hist2");
		OverLap2(hist1, hist2, 1, 0, 1, yf, "shms_hodo_beta.png", "SHMS Hodo Beta", true, outROOT);
		delete hist1; delete hist2;
	}
	
	//----- 2d -----
	if (shms_hod_2d)
	{
		//Beta vs X_fp
		hist21 = new TH2F("hist21", "SHMS Calibrated Hodo; X Focal Plane [cm]; Beta (p/E)", 100, -50, 50, 100, 0.1, 1.5);
		hist22 = new TH2F("hist22", "SHMS Uncalibrated Hodo; X Focal Plane [cm]; Beta (p/E)", 100, -50, 50, 100, 0.1, 1.5);
		tshms_replay_17243_hodo_c->Draw("P.hod.beta:P.dc.x_fp>>hist21");
		tshms_replay_17243_hodo_uc->Draw("P.hod.beta:P.dc.x_fp>>hist22");
		Compare2(hist21, hist22, -50, 1, 50, 1, "shms_hodo_beta_vs_xfp.png", "HMS Hodo Beta vs X_fp", true, outROOT);
		delete hist21;	delete hist22;

		//Beta vs Y_fp
		hist21 = new TH2F("hist21", "SHMS Calibrated Hodo; Y Focal Plane [cm]; Beta (p/E)", 100, -50, 50, 100, 0.1, 1.5);
		hist22 = new TH2F("hist22", "SHMS Uncalibrated Hodo; Y Focal Plane [cm]; Beta (p/E)", 100, -50, 50, 100, 0.1, 1.5);
		tshms_replay_17243_hodo_c->Draw("P.hod.beta:P.dc.y_fp>>hist21");
		tshms_replay_17243_hodo_uc->Draw("P.hod.beta:P.dc.y_fp>>hist22");
		Compare2(hist21, hist22, -50, 1, 50, 1, "shms_hodo_beta_vs_yfp.png", "SHMS Hodo Beta vs Y_fp", true, outROOT);
		delete hist21;	delete hist22;
	}



	//------------------------------------------------------------------------------
	//------------------------------ HMS Drift Chamber -----------------------------
	//------------------------------------------------------------------------------

	//----- Files -----
	//Calibrated
	TFile* hms_replay_4407_dc_c = new TFile("../../../../OFFLINE/CALIBRATIONS/drift_chambers/hms_replay_production_all_4407_100000_dcCalib.root", "READ");
	TTree* thms_replay_4407_dc_c = (TTree*)hms_replay_4407_dc_c->Get("T");

	//Uncalibrated
	TFile* hms_replay_4407_dc_uc = new TFile("../../../../OFFLINE/CALIBRATIONS/drift_chambers/hms_replay_production_all_4407_-1_dcUnCalib.root", "READ");
	TTree* thms_replay_4407_dc_uc = (TTree*)hms_replay_4407_dc_uc->Get("T");

	//----- 1d -----
	if (hms_dc_res)
	{
		const char* hms_dc_res_name1[12] = { "HMS DC Residual (1u1); Residual [cm]", "HMS DC Residual (1u2); Residual [cm]", "HMS DC Residual (2u1); Residual [cm]", "HMS DC Residual (2u2); Residual [cm]", "HMS DC Residual (1x1); Residual [cm]", "HMS DC Residual (1x2); Residual [cm]", "HMS DC Residual (2x1); Residual [cm]", "HMS DC Residual (2x2); Residual [cm]", "HMS DC Residual (1v1); Residual [cm]", "HMS DC Residual (1v2); Residual [cm]", "HMS DC Residual (2v1); Residual [cm]", "HMS DC Residual (2v2); Residual [cm]" };
		const char* hms_dc_res_name2[12] = { "H.dc.residualExclPlane[0]>>hist1", "H.dc.residualExclPlane[1]>>hist1", "H.dc.residualExclPlane[2]>>hist1", "H.dc.residualExclPlane[3]>>hist1", "H.dc.residualExclPlane[4]>>hist1", "H.dc.residualExclPlane[5]>>hist1", "H.dc.residualExclPlane[6]>>hist1", "H.dc.residualExclPlane[7]>>hist1", "H.dc.residualExclPlane[8]>>hist1", "H.dc.residualExclPlane[9]>>hist1", "H.dc.residualExclPlane[10]>>hist1", "H.dc.residualExclPlane[11]>>hist1" };
		const char* hms_dc_res_name3[12] = { "H.dc.residualExclPlane[0]>>hist2", "H.dc.residualExclPlane[1]>>hist2", "H.dc.residualExclPlane[2]>>hist2", "H.dc.residualExclPlane[3]>>hist2", "H.dc.residualExclPlane[4]>>hist2", "H.dc.residualExclPlane[5]>>hist2", "H.dc.residualExclPlane[6]>>hist2", "H.dc.residualExclPlane[7]>>hist2", "H.dc.residualExclPlane[8]>>hist2", "H.dc.residualExclPlane[9]>>hist2", "H.dc.residualExclPlane[10]>>hist2", "H.dc.residualExclPlane[11]>>hist2" };
		const char* hms_dc_res_name4[12] = { "hms_dc_Residual0.png", "hms_dc_Residual1.png", "hms_dc_Residual2.png", "hms_dc_Residual3.png", "hms_dc_Residual4.png", "hms_dc_Residual5.png", "hms_dc_Residual6.png", "hms_dc_Residual7.png", "hms_dc_Residual8.png", "hms_dc_Residual9.png", "hms_dc_Residual10.png", "hms_dc_Residual11.png" };
		const char* hms_dc_res_name5[12] = { "HMS DC Residual (1u1)", "HMS DC Residual (1u2)", "HMS DC Residual (2u1)", "HMS DC Residual (2u2)", "HMS DC Residual (1x1)", "HMS DC Residual (1x2)", "HMS DC Residual (2x1)", "HMS DC Residual (2x2)", "HMS DC Residual (1v1)", "HMS DC Residual (1v2)", "HMS DC Residual (2v1)", "HMS DC Residual (2v2)" };

		for (int i = 0; i < 12; i++){
			h1[0] = new TH1F("hist1", hms_dc_res_name1[i], 100, -0.3, 0.3);
			h1[1] = new TH1F("hist2", hms_dc_res_name1[i], 100, -0.3, 0.3);
			thms_replay_4407_dc_c->Draw(hms_dc_res_name2[i]);
			thms_replay_4407_dc_uc->Draw(hms_dc_res_name3[i]);
			Residuals(h1[0], h1[1], 0, 0, 0, yf, hms_dc_res_name4[i], hms_dc_res_name5[i], true, outROOT);
			delete h1[0]; delete h1[1];
		}
	}

	if (hms_dc_dist)
	{
		//Drift Distance
		const char* hms_dc_dist_name1[12] = { "HMS DC Drift Distance (1u1); Drift Distance [cm]", "HMS DC Drift Distance (1u2); Drift Distance [cm]", "HMS DC Drift Distance (2u1); Drift Distance [cm]", "HMS DC Drift Distance (2u2); Drift Distance [cm]", "HMS DC Drift Distance (1x1); Drift Distance [cm]", "HMS DC Drift Distance (1x2); Drift Distance [cm]", "HMS DC Drift Distance (2x1); Drift Distance [cm]", "HMS DC Drift Distance (2x2); Drift Distance [cm]", "HMS DC Drift Distance (1v1); Drift Distance [cm]", "HMS DC Drift Distance (1v2); Drift Distance [cm]", "HMS DC Drift Distance (2v1); Drift Distance [cm]", "HMS DC Drift Distance (2v2); Drift Distance [cm]" };
		const char* hms_dc_dist_name2[12] = { "H.dc.1u1.dist>>hist1", "H.dc.1u2.dist>>hist1", "H.dc.2u1.dist>>hist1", "H.dc.2u2.dist>>hist1", "H.dc.1x1.dist>>hist1", "H.dc.1x2.dist>>hist1", "H.dc.2x1.dist>>hist1", "H.dc.2x2.dist>>hist1", "H.dc.1v1.dist>>hist1", "H.dc.1v2.dist>>hist1", "H.dc.2v1.dist>>hist1", "H.dc.2v2.dist>>hist1" };
		const char* hms_dc_dist_name3[12] = { "H.dc.1u1.dist>>hist2", "H.dc.1u2.dist>>hist2", "H.dc.2u1.dist>>hist2", "H.dc.2u2.dist>>hist2", "H.dc.1x1.dist>>hist2", "H.dc.1x2.dist>>hist2", "H.dc.2x1.dist>>hist2", "H.dc.2x2.dist>>hist2", "H.dc.1v1.dist>>hist2", "H.dc.1v2.dist>>hist2", "H.dc.2v1.dist>>hist2", "H.dc.2v2.dist>>hist2" };
		const char* hms_dc_dist_name4[12] = { "hms_dc_dist1u1.png", "hms_dc_dist1u2.png", "hms_dc_dist2u1.png", "hms_dc_dist2u2.png", "hms_dc_dist1x1.png", "hms_dc_dist1x2.png", "hms_dc_dist2x1.png", "hms_dc_dist2x2.png", "hms_dc_dist1v1.png", "hms_dc_dist1v2.png", "hms_dc_dist2v1.png", "hms_dc_dist2v2.png" };
		const char* hms_dc_dist_name5[12] = { "HMS DC Drift Distance (1u1)", "HMS DC Drift Distance (1u2)", "HMS DC Drift Distance (2u1)", "HMS DC Drift Distance (2u2)", "HMS DC Drift Distance (1x1)", "HMS DC Drift Distance (1x2)", "HMS DC Drift Distance (2x1)", "HMS DC Drift Distance (2x2)", "HMS DC Drift Distance (1v1)", "HMS DC Drift Distance (1v2)", "HMS DC Drift Distance (2v1)", "HMS DC Drift Distance (2v2)" };
		const char* hms_dc_dist_name6[12] = { "H.dc.1u1.nhit == 1", "H.dc.1u2.nhit == 1", "H.dc.2u1.nhit == 1", "H.dc.2u2.nhit == 1", "H.dc.1x1.nhit == 1", "H.dc.1x2.nhit == 1", "H.dc.2x1.nhit == 1", "H.dc.2x2.nhit == 1", "H.dc.1v1.nhit == 1", "H.dc.1v2.nhit == 1", "H.dc.2v1.nhit == 1", "H.dc.2v2.nhit == 1" };

		for (int i = 0; i < 12; i++){
			h1[0] = new TH1F("hist1", hms_dc_dist_name1[i], 100, -0.1, 0.55);
			h1[1] = new TH1F("hist2", hms_dc_dist_name1[i], 100, -0.1, 0.55);
			thms_replay_4407_dc_c->Draw(hms_dc_dist_name2[i], hms_dc_dist_name6[i], "");
			thms_replay_4407_dc_uc->Draw(hms_dc_dist_name3[i], hms_dc_dist_name6[i], "");
			OverLap2(h1[0], h1[1], 0, 0, 0, 0, hms_dc_dist_name4[i], hms_dc_dist_name5[i], false, outROOT);
			delete h1[0]; delete h1[1];
		}
	}

	//----- 2d -----
	if (unfocused)
	{
		h2[0] = new TH2F("hist21", "HMS DC Y vs X Focal Plane Run 4407; X Focal Plane [cm]; Y Focal Plane [cm]", 100, -50, 50, 100, -50, 50);
		thms_replay_4407_dc_c->Draw("H.dc.y_fp:H.dc.x_fp>>hist21");
		Draw2d(h2[0], "hms_dc_4407_unfocused.png", "HMS DC Y vs X Focal Plane Run 4407", outROOT);
		delete h2[0];
	}

	if (hms_dc_hotwire)
	{
		h2[0] = new TH2F("hist21", "HMS Calibrated DC (1v2) Anomalous Wire; Wire Number; Drift Time [ns]", 120, 0, 120, 300, -50, 250);
		thms_replay_4407_dc_c->Draw("H.dc.1v2.time:H.dc.1v2.wirenum>>hist21");
		HotWire(h2[0], "hms_dc_anomwire_1v2.png", "HMS DC (1v2) Anomalous Wire", 64, 65, 65, 66, "Y Projection Anomalous Wire [x=65,66]", "Y Projection Normal Wire [x=64,65]", outROOT);
		delete h2[0];
	}
	
	if (hms_dc_dt_v_wn)
	{
		//Drift Time vs Wire Number
		const char* hms_dc_dt_v_wn_name1[12] = { "HMS Calibrated DC (1u1); Wire Number; Drift Time [ns]", "HMS Calibrated DC (1u2); Wire Number; Drift Time [ns]","HMS Calibrated DC (2u1); Wire Number; Drift Time [ns]","HMS Calibrated DC (2u2); Wire Number; Drift Time [ns]","HMS Calibrated DC (1x1); Wire Number; Drift Time [ns]","HMS Calibrated DC (1x2); Wire Number; Drift Time [ns]","HMS Calibrated DC (2x1); Wire Number; Drift Time [ns]","HMS Calibrated DC (2x2); Wire Number; Drift Time [ns]","HMS Calibrated DC (1v1); Wire Number; Drift Time [ns]","HMS Calibrated DC (1v2); Wire Number; Drift Time [ns]","HMS Calibrated DC (2v1); Wire Number; Drift Time [ns]","HMS Calibrated DC (2v2); Wire Number; Drift Time [ns]"};
		const char* hms_dc_dt_v_wn_name2[12] = { "H.dc.1u1.time:H.dc.1u1.wirenum>>hist21", "H.dc.1u2.time:H.dc.1u2.wirenum>>hist21","H.dc.2u1.time:H.dc.2u1.wirenum>>hist21","H.dc.2u2.time:H.dc.2u2.wirenum>>hist21","H.dc.1x1.time:H.dc.1x1.wirenum>>hist21","H.dc.1x2.time:H.dc.1x2.wirenum>>hist21","H.dc.2x1.time:H.dc.2x1.wirenum>>hist21","H.dc.2x2.time:H.dc.2x2.wirenum>>hist21","H.dc.1v1.time:H.dc.1v1.wirenum>>hist21","H.dc.1v2.time:H.dc.1v2.wirenum>>hist21","H.dc.2v1.time:H.dc.2v1.wirenum>>hist21","H.dc.2v2.time:H.dc.2v2.wirenum>>hist21"};
		const char* hms_dc_dt_v_wn_name3[12] = { "HMS Uncalibrated DC (1u1); Wire Number; Drift Time [ns]", "HMS Uncalibrated DC (1u2); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (2u1); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (2u2); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (1x1); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (1x2); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (2x1); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (2x2); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (1v1); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (1v2); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (2v1); Wire Number; Drift Time [ns]","HMS Uncalibrated DC (2v2); Wire Number; Drift Time [ns]"};
		const char* hms_dc_dt_v_wn_name4[12] = { "H.dc.1u1.time:H.dc.1u1.wirenum>>hist22", "H.dc.1u2.time:H.dc.1u2.wirenum>>hist22","H.dc.2u1.time:H.dc.2u1.wirenum>>hist22","H.dc.2u2.time:H.dc.2u2.wirenum>>hist22","H.dc.1x1.time:H.dc.1x1.wirenum>>hist22","H.dc.1x2.time:H.dc.1x2.wirenum>>hist22","H.dc.2x1.time:H.dc.2x1.wirenum>>hist22","H.dc.2x2.time:H.dc.2x2.wirenum>>hist22","H.dc.1v1.time:H.dc.1v1.wirenum>>hist22","H.dc.1v2.time:H.dc.1v2.wirenum>>hist22","H.dc.2v1.time:H.dc.2v1.wirenum>>hist22","H.dc.2v2.time:H.dc.2v2.wirenum>>hist22"};
		const char* hms_dc_dt_v_wn_name5[12] = { "hms_dc_wn_v_dt_1u1.png", "hms_dc_wn_v_dt_1u2.png","hms_dc_wn_v_dt_2u1.png","hms_dc_wn_v_dt_2u2.png","hms_dc_wn_v_dt_1x1.png","hms_dc_wn_v_dt_1x2.png","hms_dc_wn_v_dt_2x1.png","hms_dc_wn_v_dt_2x2.png","hms_dc_wn_v_dt_1v1.png","hms_dc_wn_v_dt_1v2.png","hms_dc_wn_v_dt_2v1.png","hms_dc_wn_v_dt_2v2.png"};
		const char* hms_dc_dt_v_wn_name6[12] = { "HMS DC Wire Number (1u1) vs Drift Time", "HMS DC Wire Number (1u2) vs Drift Time","HMS DC Wire Number (2u1) vs Drift Time","HMS DC Wire Number (2u2) vs Drift Time","HMS DC Wire Number (1x1) vs Drift Time","HMS DC Wire Number (1x2) vs Drift Time","HMS DC Wire Number (2x1) vs Drift Time","HMS DC Wire Number (2x2) vs Drift Time","HMS DC Wire Number (1v1) vs Drift Time","HMS DC Wire Number (1v2) vs Drift Time","HMS DC Wire Number (2v1) vs Drift Time","HMS DC Wire Number (2v2) vs Drift Time"};

		for (int i = 0; i < 12; i++){
			h2[0] = new TH2F("hist21", hms_dc_dt_v_wn_name1[i], 120, 0, 120, 300, -50, 250);
			h2[1] = new TH2F("hist22", hms_dc_dt_v_wn_name3[i], 120, 0, 120, 300, -50, 250);
			thms_replay_4407_dc_c->Draw(hms_dc_dt_v_wn_name2[i]);
			thms_replay_4407_dc_uc->Draw(hms_dc_dt_v_wn_name4[i]);
			Compare2(h2[0], h2[1], 0, 0, 0, 0, hms_dc_dt_v_wn_name5[i], hms_dc_dt_v_wn_name6[i], false, outROOT);
			delete h2[0]; delete h2[1];
		}
	}
	


	//-------------------------------------------------------------------------------
	//------------------------------ SHMS Drift Chamber -----------------------------
	//-------------------------------------------------------------------------------

	//----- Files -----
	//Calibrated
	TFile* shms_replay_17243_dc_c = new TFile("../../../../OFFLINE/CALIBRATIONS/drift_chambers/shms_replay_production_all_17243_100000_dcCalib.root", "READ");
	TTree* tshms_replay_17243_dc_c = (TTree*)shms_replay_17243_dc_c->Get("T");

	//Uncalibrated
	TFile* shms_replay_17243_dc_uc = new TFile("../../../../OFFLINE/CALIBRATIONS/drift_chambers/shms_replay_production_all_17243_-1_dcUnCalib.root", "READ");
	TTree* tshms_replay_17243_dc_uc = (TTree*)shms_replay_17243_dc_uc->Get("T");

	//----- 1d -----
	if (shms_dc_res)
	{
		//Residuals
		const char* shms_dc_res_name1[12] = { "SHMS DC Residual (1u1); Residual [cm]", "SHMS DC Residual (1u2); Residual [cm]", "SHMS DC Residual (2u1); Residual [cm]", "SHMS DC Residual (2u2); Residual [cm]", "SHMS DC Residual (1x1); Residual [cm]", "SHMS DC Residual (1x2); Residual [cm]", "SHMS DC Residual (2x1); Residual [cm]", "SHMS DC Residual (2x2); Residual [cm]", "SHMS DC Residual (1v1); Residual [cm]", "SHMS DC Residual (1v2); Residual [cm]", "SHMS DC Residual (2v1); Residual [cm]", "SHMS DC Residual (2v2); Residual [cm]" };
		const char* shms_dc_res_name2[12] = { "P.dc.residualExclPlane[0]>>hist1", "P.dc.residualExclPlane[1]>>hist1", "P.dc.residualExclPlane[2]>>hist1", "P.dc.residualExclPlane[3]>>hist1", "P.dc.residualExclPlane[4]>>hist1", "P.dc.residualExclPlane[5]>>hist1", "P.dc.residualExclPlane[6]>>hist1", "P.dc.residualExclPlane[7]>>hist1", "P.dc.residualExclPlane[8]>>hist1", "P.dc.residualExclPlane[9]>>hist1", "P.dc.residualExclPlane[10]>>hist1", "P.dc.residualExclPlane[11]>>hist1" };
		const char* shms_dc_res_name3[12] = { "P.dc.residualExclPlane[0]>>hist2", "P.dc.residualExclPlane[1]>>hist2", "P.dc.residualExclPlane[2]>>hist2", "P.dc.residualExclPlane[3]>>hist2", "P.dc.residualExclPlane[4]>>hist2", "P.dc.residualExclPlane[5]>>hist2", "P.dc.residualExclPlane[6]>>hist2", "P.dc.residualExclPlane[7]>>hist2", "P.dc.residualExclPlane[8]>>hist2", "P.dc.residualExclPlane[9]>>hist2", "P.dc.residualExclPlane[10]>>hist2", "P.dc.residualExclPlane[11]>>hist2" };
		const char* shms_dc_res_name4[12] = { "shms_dc_Residual0.png", "shms_dc_Residual1.png", "shms_dc_Residual2.png", "shms_dc_Residual3.png", "shms_dc_Residual4.png", "shms_dc_Residual5.png", "shms_dc_Residual6.png", "shms_dc_Residual7.png", "shms_dc_Residual8.png", "shms_dc_Residual9.png", "shms_dc_Residual10.png", "shms_dc_Residual11.png" };
		const char* shms_dc_res_name5[12] = { "SHMS DC Residual (1u1)", "SHMS DC Residual (1u2)", "SHMS DC Residual (2u1)", "SHMS DC Residual (2u2)", "SHMS DC Residual (1x1)", "SHMS DC Residual (1x2)", "SHMS DC Residual (2x1)", "SHMS DC Residual (2x2)", "SHMS DC Residual (1v1)", "SHMS DC Residual (1v2)", "SHMS DC Residual (2v1)", "SHMS DC Residual (2v2)" };

		for (int i = 0; i < 12; i++){
			h1[0] = new TH1F("hist1", shms_dc_res_name1[i], 100, -0.3, 0.3);
			h1[1] = new TH1F("hist2", shms_dc_res_name1[i], 100, -0.3, 0.3);
			tshms_replay_17243_dc_c->Draw(shms_dc_res_name2[i]);
			tshms_replay_17243_dc_uc->Draw(shms_dc_res_name3[i]);
			Residuals(h1[0], h1[1], 0, 0, 0, yf, shms_dc_res_name4[i], shms_dc_res_name5[i], true, outROOT);
			delete h1[1]; delete h1[0];
		}
	}
	
	if (shms_dc_dist)
	{
		//Drift Distance
		const char* shms_dc_dist_name1[12] = { "SHMS DC Drift Distance (1u1); Drift Distance [cm]", "SHMS DC Drift Distance (1u2); Drift Distance [cm]", "SHMS DC Drift Distance (2u1); Drift Distance [cm]", "SHMS DC Drift Distance (2u2); Drift Distance [cm]", "SHMS DC Drift Distance (1x1); Drift Distance [cm]", "SHMS DC Drift Distance (1x2); Drift Distance [cm]", "SHMS DC Drift Distance (2x1); Drift Distance [cm]", "SHMS DC Drift Distance (2x2); Drift Distance [cm]", "SHMS DC Drift Distance (1v1); Drift Distance [cm]", "SHMS DC Drift Distance (1v2); Drift Distance [cm]", "SHMS DC Drift Distance (2v1); Drift Distance [cm]", "SHMS DC Drift Distance (2v2); Drift Distance [cm]" };
		const char* shms_dc_dist_name2[12] = { "P.dc.1u1.dist>>hist1", "P.dc.1u2.dist>>hist1", "P.dc.2u1.dist>>hist1", "P.dc.2u2.dist>>hist1", "P.dc.1x1.dist>>hist1", "P.dc.1x2.dist>>hist1", "P.dc.2x1.dist>>hist1", "P.dc.2x2.dist>>hist1", "P.dc.1v1.dist>>hist1", "P.dc.1v2.dist>>hist1", "P.dc.2v1.dist>>hist1", "P.dc.2v2.dist>>hist1" };
		const char* shms_dc_dist_name3[12] = { "P.dc.1u1.dist>>hist2", "P.dc.1u2.dist>>hist2", "P.dc.2u1.dist>>hist2", "P.dc.2u2.dist>>hist2", "P.dc.1x1.dist>>hist2", "P.dc.1x2.dist>>hist2", "P.dc.2x1.dist>>hist2", "P.dc.2x2.dist>>hist2", "P.dc.1v1.dist>>hist2", "P.dc.1v2.dist>>hist2", "P.dc.2v1.dist>>hist2", "P.dc.2v2.dist>>hist2" };
		const char* shms_dc_dist_name4[12] = { "shms_dc_dist1u1.png", "shms_dc_dist1u2.png", "shms_dc_dist2u1.png", "shms_dc_dist2u2.png", "shms_dc_dist1x1.png", "shms_dc_dist1x2.png", "shms_dc_dist2x1.png", "shms_dc_dist2x2.png", "shms_dc_dist1v1.png", "shms_dc_dist1v2.png", "shms_dc_dist2v1.png", "shms_dc_dist2v2.png" };
		const char* shms_dc_dist_name5[12] = { "SHMS DC Drift Distance (1u1)", "SHMS DC Drift Distance (1u2)", "SHMS DC Drift Distance (2u1)", "SHMS DC Drift Distance (2u2)", "SHMS DC Drift Distance (1x1)", "SHMS DC Drift Distance (1x2)", "SHMS DC Drift Distance (2x1)", "SHMS DC Drift Distance (2x2)", "SHMS DC Drift Distance (1v1)", "SHMS DC Drift Distance (1v2)", "SHMS DC Drift Distance (2v1)", "SHMS DC Drift Distance (2v2)" };
		const char* shms_dc_dist_name6[12] = { "P.dc.1u1.nhit == 1", "P.dc.1u2.nhit == 1", "P.dc.2u1.nhit == 1", "P.dc.2u2.nhit == 1", "P.dc.1x1.nhit == 1", "P.dc.1x2.nhit == 1", "P.dc.2x1.nhit == 1", "P.dc.2x2.nhit == 1", "P.dc.1v1.nhit == 1", "P.dc.1v2.nhit == 1", "P.dc.2v1.nhit == 1", "P.dc.2v2.nhit == 1" };

		for (int i = 0; i < 12; i++){
			h1[0] = new TH1F("hist1", shms_dc_dist_name1[i], 100, -0.1, 0.55);
			h1[1] = new TH1F("hist2", shms_dc_dist_name1[i], 100, -0.1, 0.55);
			tshms_replay_17243_dc_c->Draw(shms_dc_dist_name2[i], shms_dc_dist_name6[i], "");
			tshms_replay_17243_dc_uc->Draw(shms_dc_dist_name3[i], shms_dc_dist_name6[i], "");
			OverLap2(h1[0], h1[1], 0, 0, 0, 0, shms_dc_dist_name4[i], shms_dc_dist_name5[i], false, outROOT);
			delete h1[0]; delete h1[1];
		}
	}
	
	//----- 2d -----
	if (unfocused)
	{
		h2[0] = new TH2F("hist21", "SHMS DC Y vs X Focal Plane Run 17243; X Focal Plane [cm]; Y Focal Plane [cm]", 100, -50, 50, 100, -50, 50);
		tshms_replay_17243_dc_c->Draw("P.dc.y_fp:P.dc.x_fp>>hist21");
		Draw2d(h2[0], "shms_dc_17243_unfocused.png", "SHMS DC Y vs X Focal Plane Run 17243", outROOT);
		delete h2[0];
	}
	
	if (shms_dc_hotwire)
	{
		h2[0] = new TH2F("hist21", "SHMS Calibrated DC (1u2) Anomalous Wire; Wire Number; Drift Time [ns]", 120, 0, 120, 300, -50, 250);
		tshms_replay_17243_dc_c->Draw("P.dc.1u2.time:P.dc.1u2.wirenum>>hist21");
		HotWire(h2[0], "shms_dc_anomwire_1u2.png", "SHMS DC (1u2) Anomalous Wire", 25, 26, 24, 25, "Y Projection Anomalous Wire [x=24,25]", "Y Projection Normal Wire [x=25,26]", outROOT);
		delete h2[0];
	}
	
	if (shms_dc_wn_v_dt)//1u2 hot wire
	{
		//Drift Time vs Wire Number
		const char* shms_dc_dt_v_wn_name1[12] = { "SHMS Calibrated DC (1u1); Wire Number; Drift Time [ns]", "SHMS Calibrated DC (1u2); Wire Number; Drift Time [ns]","SHMS Calibrated DC (2u1); Wire Number; Drift Time [ns]","SHMS Calibrated DC (2u2); Wire Number; Drift Time [ns]","SHMS Calibrated DC (1x1); Wire Number; Drift Time [ns]","SHMS Calibrated DC (1x2); Wire Number; Drift Time [ns]","SHMS Calibrated DC (2x1); Wire Number; Drift Time [ns]","SHMS Calibrated DC (2x2); Wire Number; Drift Time [ns]","SHMS Calibrated DC (1v1); Wire Number; Drift Time [ns]","SHMS Calibrated DC (1v2); Wire Number; Drift Time [ns]","SHMS Calibrated DC (2v1); Wire Number; Drift Time [ns]","SHMS Calibrated DC (2v2); Wire Number; Drift Time [ns]" };
		const char* shms_dc_dt_v_wn_name2[12] = { "P.dc.1u1.time:P.dc.1u1.wirenum>>hist21", "P.dc.1u2.time:P.dc.1u2.wirenum>>hist21","P.dc.2u1.time:P.dc.2u1.wirenum>>hist21","P.dc.2u2.time:P.dc.2u2.wirenum>>hist21","P.dc.1x1.time:P.dc.1x1.wirenum>>hist21","P.dc.1x2.time:P.dc.1x2.wirenum>>hist21","P.dc.2x1.time:P.dc.2x1.wirenum>>hist21","P.dc.2x2.time:P.dc.2x2.wirenum>>hist21","P.dc.1v1.time:P.dc.1v1.wirenum>>hist21","P.dc.1v2.time:P.dc.1v2.wirenum>>hist21","P.dc.2v1.time:P.dc.2v1.wirenum>>hist21","P.dc.2v2.time:P.dc.2v2.wirenum>>hist21" };
		const char* shms_dc_dt_v_wn_name3[12] = { "SHMS Uncalibrated DC (1u1); Wire Number; Drift Time [ns]", "SHMS Uncalibrated DC (1u2); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (2u1); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (2u2); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (1x1); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (1x2); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (2x1); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (2x2); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (1v1); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (1v2); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (2v1); Wire Number; Drift Time [ns]","SHMS Uncalibrated DC (2v2); Wire Number; Drift Time [ns]" };
		const char* shms_dc_dt_v_wn_name4[12] = { "P.dc.1u1.time:P.dc.1u1.wirenum>>hist22", "P.dc.1u2.time:P.dc.1u2.wirenum>>hist22","P.dc.2u1.time:P.dc.2u1.wirenum>>hist22","P.dc.2u2.time:P.dc.2u2.wirenum>>hist22","P.dc.1x1.time:P.dc.1x1.wirenum>>hist22","P.dc.1x2.time:P.dc.1x2.wirenum>>hist22","P.dc.2x1.time:P.dc.2x1.wirenum>>hist22","P.dc.2x2.time:P.dc.2x2.wirenum>>hist22","P.dc.1v1.time:P.dc.1v1.wirenum>>hist22","P.dc.1v2.time:P.dc.1v2.wirenum>>hist22","P.dc.2v1.time:P.dc.2v1.wirenum>>hist22","P.dc.2v2.time:P.dc.2v2.wirenum>>hist22" };
		const char* shms_dc_dt_v_wn_name5[12] = { "shms_dc_wn_v_dt_1u1.png", "shms_dc_wn_v_dt_1u2.png","shms_dc_wn_v_dt_2u1.png","shms_dc_wn_v_dt_2u2.png","shms_dc_wn_v_dt_1x1.png","shms_dc_wn_v_dt_1x2.png","shms_dc_wn_v_dt_2x1.png","shms_dc_wn_v_dt_2x2.png","shms_dc_wn_v_dt_1v1.png","shms_dc_wn_v_dt_1v2.png","shms_dc_wn_v_dt_2v1.png","shms_dc_wn_v_dt_2v2.png" };
		const char* shms_dc_dt_v_wn_name6[12] = { "SHMS DC Wire Number (1u1) vs Drift Time", "SHMS DC Wire Number (1u2) vs Drift Time","SHMS DC Wire Number (2u1) vs Drift Time","SHMS DC Wire Number (2u2) vs Drift Time","SHMS DC Wire Number (1x1) vs Drift Time","SHMS DC Wire Number (1x2) vs Drift Time","SHMS DC Wire Number (2x1) vs Drift Time","SHMS DC Wire Number (2x2) vs Drift Time","SHMS DC Wire Number (1v1) vs Drift Time","SHMS DC Wire Number (1v2) vs Drift Time","SHMS DC Wire Number (2v1) vs Drift Time","SHMS DC Wire Number (2v2) vs Drift Time" };

		for (int i = 0; i < 12; i++){
			h2[0] = new TH2F("hist21", shms_dc_dt_v_wn_name1[i], 120, 0, 120, 300, -50, 250);
			h2[1] = new TH2F("hist22", shms_dc_dt_v_wn_name3[i], 120, 0, 120, 300, -50, 250);
			tshms_replay_17243_dc_c->Draw(shms_dc_dt_v_wn_name2[i]);
			tshms_replay_17243_dc_uc->Draw(shms_dc_dt_v_wn_name4[i]);
			Compare2(h2[0], h2[1], 0, 0, 0, 0, shms_dc_dt_v_wn_name5[i], shms_dc_dt_v_wn_name6[i], false, outROOT);
			delete h2[0]; delete h2[1];
		}
	}
	


	//----------------------------------------------------------------------------
	//------------------------------ HMS Calorimeter -----------------------------
	//----------------------------------------------------------------------------

	//----- Files -----
	//Calibrated
	TFile* hms_replay_14967_c = new TFile("../../calcalib/cafe_replay_calcalib_14967_c.root", "READ");
	TTree* thms_replay_14967_c = (TTree*)hms_replay_14967_c->Get("T");

	//Uncalibrated
	TFile* hms_replay_14967_uc = new TFile("../../calcalib/cafe_replay_calcalib_14967_uc.root", "READ");
	TTree* thms_replay_14967_uc = (TTree*)hms_replay_14967_uc->Get("T");
	
	//----- 1d -----
	if (hms_cal_1d)
	{
		//eTrkNorm
		hist1 = new TH1F("hist1", "HMS Cal eTrkNorm; eTrkNorm [E_{dep}/p_{track}]", 100, 0.1, 2);
		hist2 = new TH1F("hist2", "HMS Cal eTrkNorm; eTrkNorm [E_{dep}/p_{track}]", 100, 0.1, 2);
		thms_replay_14967_c->Draw("H.cal.etottracknorm>>hist1");
		thms_replay_14967_uc->Draw("H.cal.etottracknorm>>hist2");
		OverLap2(hist1, hist2, 1, 0, 1, yf, "hms_cal_eTrkNorm.png", "HMS Cal eTrkNorm", true, outROOT);
		delete hist1;	delete hist2;
	}

	//----- 2d -----
	if (unfocused)
	{
		h2[0] = new TH2F("hist21", "HMS DC Y vs X Focal Plane Run 14967; X Focal Plane [cm]; Y Focal Plane [cm]", 100, -50, 50, 100, -50, 50);
		thms_replay_14967_c->Draw("H.dc.y_fp:H.dc.x_fp>>hist21");
		Draw2d(h2[0], "hms_dc_14967_unfocused.png", "HMS DC Y vs X Focal Plane Run 14967", outROOT);
		delete h2[0];
	}

	if (hms_cal_2d)
	{
		//eTrkNorm vs X_fp
		hist21 = new TH2F("hist21", "HMS Calibrated Cal; X Focal Plane [cm]; eTrkNorm [E_{dep}/p_{track}]", 100, -50, 50, 100, 0.1, 2);
		hist22 = new TH2F("hist22", "HMS Uncalibrated Cal; X Focal Plane [cm]; eTrkNorm [E_{dep}/p_{track}]", 100, -50, 50, 100, 0.1, 2);
		thms_replay_14967_c->Draw("H.cal.etottracknorm:H.dc.x_fp>>hist21");
		thms_replay_14967_uc->Draw("H.cal.etottracknorm:H.dc.x_fp>>hist22");
		Compare2(hist21, hist22, -50, 1, 50, 1, "hms_cal_etotTrkNorm_vs_xfp.png", "HMS Cal etotTrkNorm vs X_fp", true, outROOT);
		delete hist21;	delete hist22;

		//eTrkNorm vs Y_fp
		hist21 = new TH2F("hist21", "HMS Calibrated Cal; Y Focal Plane [cm]; eTrkNorm [E_{dep}/p_{track}]", 100, -50, 50, 100, 0.1, 2);
		hist22 = new TH2F("hist22", "HMS Uncalibrated Cal; Y Focal Plane [cm]; eTrkNorm [E_{dep}/p_{track}]", 100, -50, 50, 100, 0.1, 2);
		thms_replay_14967_c->Draw("H.cal.etottracknorm:H.dc.y_fp>>hist21");
		thms_replay_14967_uc->Draw("H.cal.etottracknorm:H.dc.y_fp>>hist22");
		Compare2(hist21, hist22, -50, 1, 50, 1, "hms_cal_etotTrkNorm_vs_yfp.png", "HMS Cal etotTrkNorm vs Y_fp", true, outROOT);
		delete hist21;	delete hist22;

		//eTrkNorm vs Momentum Acceptanace
		hist21 = new TH2F("hist21", "HMS Calibrated Cal; Momentum Acceptance [%]; eTrkNorm [E_{dep}/p_{track}]", 100, -30, 30, 100, 0.1, 2);
		hist22 = new TH2F("hist22", "HMS Uncalibrated Cal; Momentum Acceptance [%]; eTrkNorm [E_{dep}/p_{track}]", 100, -30, 30, 100, 0.1, 2);
		thms_replay_14967_c->Draw("H.cal.etottracknorm:H.gtr.dp>>hist21");
		thms_replay_14967_uc->Draw("H.cal.etottracknorm:H.gtr.dp>>hist22");
		Compare2(hist21, hist22, -30, 1, 30, 1, "hms_cal_etotTrkNorm_vs_dp.png", "HMS Cal etotTrkNorm vs Momentum Acceptance", true, outROOT);
		delete hist21;	delete hist22;
	}
	


	//-----------------------------------------------------------------------------
	//------------------------------ SHMS Calorimeter -----------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	//Calibrated
	TFile* shms_replay_16962_c = new TFile("../../calcalib/cafe_replay_calcalib_16962_c.root", "READ");
	TTree* tshms_replay_16962_c = (TTree*)shms_replay_16962_c->Get("T");
	
	//Uncalibrated
	TFile* shms_replay_16962_uc = new TFile("../../calcalib/cafe_replay_calcalib_16962_uc.root", "READ");
	TTree* tshms_replay_16962_uc = (TTree*)shms_replay_16962_uc->Get("T");

	//----- 1d -----
	if (shms_cal_1d)
	{
		//EtotTrkNorm
		hist1 = new TH1F("hist1", "SHMS Cal eTrkNorm; eTrkNorm [E_{dep}/p_{track}]", 100, 0.1, 2);
		hist2 = new TH1F("hist2", "SHMS Cal eTrkNorm; eTrkNorm [E_{dep}/p_{track}]", 100, 0.1, 2);
		tshms_replay_16962_c->Draw("P.cal.etottracknorm>>hist1");
		tshms_replay_16962_uc->Draw("P.cal.etottracknorm>>hist2");
		OverLap2(hist1, hist2, 1, 0, 1, yf, "shms_cal_eTrkNorm.png", "SHMS Cal eTrkNorm", true, outROOT);
		delete hist1;	delete hist2;
	}

	//----- 2d -----
	if (unfocused)
	{
		h2[0] = new TH2F("hist21", "SHMS DC Y vs X Focal Plane Run 16962; X Focal Plane [cm]; Y Focal Plane [cm]", 100, -50, 50, 100, -50, 50);
		tshms_replay_16962_c->Draw("P.dc.y_fp:P.dc.x_fp>>hist21");
		Draw2d(h2[0], "shms_dc_16962_focused.png", "SHMS DC Y vs X Focal Plane Run 16962", outROOT);
		delete h2[0];
	}
	
	if (shms_cal_2d)
	{
		//SHMS Calorimeter etotTrkNorm vs X_fp
		hist21 = new TH2F("hist21", "Calibrated; X Focal Plane [cm]; eTrkNorm [E_{dep}/p_{track}]", 100, -50, 50, 100, 0.1, 2);
		hist22 = new TH2F("hist22", "Uncalibrated; X Focal Plane [cm]; eTrkNorm [E_{dep}/p_{track}]", 100, -50, 50, 100, 0.1, 2);
		tshms_replay_16962_c->Draw("P.cal.etottracknorm:P.dc.x_fp>>hist21");
		tshms_replay_16962_uc->Draw("P.cal.etottracknorm:P.dc.x_fp>>hist22");
		Compare2(hist21, hist22, -50, 1, 50, 1, "shms_cal_etotTrkNorm_vs_xfp.png", "SHMS Cal etotTrkNorm vs X_fp", true, outROOT);
		delete hist21;	delete hist22;

		//SHMS Calorimeter etotTrkNorm vs Y_fp
		hist21 = new TH2F("hist21", "Calibrated; Y Focal Plane [cm]; eTrkNorm [E_{dep}/p_{track}]", 100, -50, 50, 100, 0.1, 2);
		hist22 = new TH2F("hist22", "Uncalibrated; Y Focal Plane [cm]; eTrkNorm [E_{dep}/p_{track}]", 100, -50, 50, 100, 0.1, 2);
		tshms_replay_16962_c->Draw("P.cal.etottracknorm:P.dc.y_fp>>hist21");
		tshms_replay_16962_uc->Draw("P.cal.etottracknorm:P.dc.y_fp>>hist22");
		Compare2(hist21, hist22, -50, 1, 50, 1, "shms_cal_etotTrkNorm_vs_yfp.png", "SHMS Cal etotTrkNorm vs Y_fp", true, outROOT);
		delete hist21;	delete hist22;

		//SHMS Calorimeter etotTrkNorm vs Momentum Acceptanace
		hist21 = new TH2F("hist21", "Calibrated; Momentum Acceptance [%]; eTrkNorm [E_{dep}/p_{track}]", 100, -30, 30, 100, 0.1, 2);
		hist22 = new TH2F("hist22", "Uncalibrated; Momentum Acceptance [%]; eTrkNorm [E_{dep}/p_{track}]", 100, -30, 30, 100, 0.1, 2);
		tshms_replay_16962_c->Draw("P.cal.etottracknorm:P.gtr.dp>>hist21");
		tshms_replay_16962_uc->Draw("P.cal.etottracknorm:P.gtr.dp>>hist22");
		Compare2(hist21, hist22, -30, 1, 30, 1, "shms_cal_etotTrkNorm_vs_dp.png", "SHMS Cal etotTrkNorm vs Momentum Acceptance", true, outROOT);
		delete hist21;	delete hist22;
	}
	


	//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
	//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= Run Files =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
	//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
	
	//-----------------------------------------------------------------------------
	//-------------------------------- HMS Light MF -------------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	TFile* hms_LM_D2_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/hms_Calib_histos16975.root", "READ");
	TFile* hms_LM_Be9_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/hms_Calib_histos16983.root", "READ");
	TFile* hms_LM_B10_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/hms_Calib_histos16984.root", "READ");
	TFile* hms_LM_B11_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/hms_Calib_histos16986.root", "READ");
	TFile* hms_LM_C12_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/hms_Calib_histos16977.root", "READ");

	//----- 1d -----
	if (hms_light_mf_1d)
	{
		//Beta
		hist1 = (TH1F*)hms_LM_D2_calib->Get("hHod_Beta");
		hist2 = (TH1F*)hms_LM_Be9_calib->Get("hHod_Beta");
		hist3 = (TH1F*)hms_LM_B10_calib->Get("hHod_Beta");
		hist4 = (TH1F*)hms_LM_B11_calib->Get("hHod_Beta");
		hist5 = (TH1F*)hms_LM_C12_calib->Get("hHod_Beta");
		OverLap5(hist1, hist2, hist3, hist4, hist5, 0.8888, 0, 0.8888, yf, "hms_LM_beta.png", "HMS MF Hodo Beta", "Beta [p/E]", true, outROOT);
		delete hist1; delete hist2; delete hist3; delete hist4; delete hist5;
	}
	
	//----- 2d -----
	if (hms_light_mf_2d)
	{
		//Residuals
		g1 = (TGraph*)hms_LM_D2_calib->Get("Graph;2");
		g2 = (TGraph*)hms_LM_Be9_calib->Get("Graph;2");
		g3 = (TGraph*)hms_LM_B10_calib->Get("Graph;2");
		g4 = (TGraph*)hms_LM_B11_calib->Get("Graph;2");
		g5 = (TGraph*)hms_LM_C12_calib->Get("Graph;2");
		OverLay5(g1, g2, g3, g4, g5, "hms_LM_residuals.png", "HMS MF DC Residuals", "HMS DC Planes", true, outROOT);
		delete g1; delete g2; delete g3; delete g4; delete g5;
	}
	


	//-----------------------------------------------------------------------------
	//------------------------------- SHMS Light MF -------------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	TFile* shms_LM_D2_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/shms_Calib_histos16975.root", "READ");
	TFile* shms_LM_Be9_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/shms_Calib_histos16983.root", "READ");
	TFile* shms_LM_B10_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/shms_Calib_histos16984.root", "READ");
	TFile* shms_LM_B11_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/shms_Calib_histos16986.root", "READ");
	TFile* shms_LM_C12_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_MF/shms_Calib_histos16977.root", "READ");

	//----- 1d -----
	if (shms_light_mf_1d)
	{
		//Beta
		hist1 = (TH1F*)shms_LM_D2_calib->Get("pHod_Beta");
		hist2 = (TH1F*)shms_LM_Be9_calib->Get("pHod_Beta");
		hist3 = (TH1F*)shms_LM_B10_calib->Get("pHod_Beta");
		hist4 = (TH1F*)shms_LM_B11_calib->Get("pHod_Beta");
		hist5 = (TH1F*)shms_LM_C12_calib->Get("pHod_Beta");
		OverLap5(hist1, hist2, hist3, hist4, hist5, 1, 0, 1, yf, "shms_LM_beta.png", "SHMS MF Hodo Beta", "Beta [p/E]", false, outROOT);
		delete hist1; delete hist2; delete hist3; delete hist4; delete hist5;
		
		//etotTrkNorm
		hist1 = (TH1F*)shms_LM_D2_calib->Get("pCal_eTrkNorm");
		hist2 = (TH1F*)shms_LM_Be9_calib->Get("pCal_eTrkNorm");
		hist3 = (TH1F*)shms_LM_B10_calib->Get("pCal_eTrkNorm");
		hist4 = (TH1F*)shms_LM_B11_calib->Get("pCal_eTrkNorm");
		hist5 = (TH1F*)shms_LM_C12_calib->Get("pCal_eTrkNorm");
		OverLap5(hist1, hist2, hist3, hist4, hist5, 1, 0, 1, yf, "shms_LM_eTrkNorm.png", "SHMS MF Cal eTrkNorm", "eTrkNorm [E_{dep}/p_{track}]", false, outROOT);
		delete hist1; delete hist2; delete hist3; delete hist4; delete hist5;
	}
	
	//----- 2d -----
	if (shms_light_mf_2d)
	{
		//Residuals
		g1 = (TGraph*)shms_LM_D2_calib->Get("Graph;2");
		g2 = (TGraph*)shms_LM_Be9_calib->Get("Graph;2");
		g3 = (TGraph*)shms_LM_B10_calib->Get("Graph;2");
		g4 = (TGraph*)shms_LM_B11_calib->Get("Graph;2");
		g5 = (TGraph*)shms_LM_C12_calib->Get("Graph;2");
		OverLay5(g1, g2, g3, g4, g5, "shms_LM_residuals.png", "SHMS MF DC Residuals", "SHMS DC Planes", false, outROOT);
		delete g1; delete g2; delete g3; delete g4; delete g5;
	}
	
	

	//-----------------------------------------------------------------------------
	//-------------------------------- HMS Heavy MF -------------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	TFile* hms_HM_Ca40_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_MF/hms_Calib_histos16980.root", "READ");
	TFile* hms_HM_Ca48_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_MF/hms_Calib_histos17096.root", "READ");
	TFile* hms_HM_Fe54_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_MF/hms_Calib_histos16982.root", "READ");

	//----- 1d -----
	if (hms_heavy_mf_1d)
	{
		//Beta
		hist1 = (TH1F*)hms_HM_Ca40_calib->Get("hHod_Beta");
		hist2 = (TH1F*)hms_HM_Ca48_calib->Get("hHod_Beta");
		hist3 = (TH1F*)hms_HM_Fe54_calib->Get("hHod_Beta");
		OverLap3(hist1, hist2, hist3, 0.8888, 0, 0.8888, yf, "hms_HM_beta.png", "HMS MF Hodo Beta", "Beta [p/E]", true, outROOT);
		delete hist1; delete hist2; delete hist3;
	}
	
	//----- 2d -----
	if (hms_heavy_mf_2d)
	{
		//Residuals
		g1 = (TGraph*)hms_HM_Ca40_calib->Get("Graph;2");
		g2 = (TGraph*)hms_HM_Ca48_calib->Get("Graph;2");
		g3 = (TGraph*)hms_HM_Fe54_calib->Get("Graph;2");
		OverLay3(g1, g2, g3, "hms_HM_residuals.png", "HMS MF DC Residuals", "HMS DC Planes", true, outROOT);
		delete g1; delete g2; delete g3;
	}


	//-----------------------------------------------------------------------------
	//------------------------------- SHMS Heavy MF -------------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	TFile* shms_HM_Ca40_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_MF/shms_Calib_histos16980.root", "READ");
	TFile* shms_HM_Ca48_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_MF/shms_Calib_histos17096.root", "READ");
	TFile* shms_HM_Fe54_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_MF/shms_Calib_histos16982.root", "READ");

	//----- 1d -----
	if (shms_heavy_mf_1d)
	{
		//Beta
		hist1 = (TH1F*)shms_HM_Ca40_calib->Get("pHod_Beta");
		hist2 = (TH1F*)shms_HM_Ca48_calib->Get("pHod_Beta");
		hist3 = (TH1F*)shms_HM_Fe54_calib->Get("pHod_Beta");
		OverLap3(hist1, hist2, hist3, 1, 0, 1, yf, "shms_HM_beta.png", "SHMS MF Hodo Beta", "Beta [p/E]", true, outROOT);
		delete hist1; delete hist2; delete hist3;
		
		//eTrkNorm
		hist1 = (TH1F*)shms_HM_Ca40_calib->Get("pCal_eTrkNorm");
		hist2 = (TH1F*)shms_HM_Ca48_calib->Get("pCal_eTrkNorm");
		hist3 = (TH1F*)shms_HM_Fe54_calib->Get("pCal_eTrkNorm");
		OverLap3(hist1, hist2, hist3, 1, 0, 1, yf, "shms_HM_eTrkNorm.png", "SHMS MF Cal eTrkNorm", "eTrkNorm [E_{dep}/[p_{track}]", true, outROOT);
		delete hist1;  delete hist2; delete hist3;
	}

	//----- 2d -----
	if (shms_heavy_mf_2d)
	{
		//Residuals
		g1 = (TGraph*)shms_HM_Ca40_calib->Get("Graph;2");
		g2 = (TGraph*)shms_HM_Ca48_calib->Get("Graph;2");
		g3 = (TGraph*)shms_HM_Fe54_calib->Get("Graph;2");
		OverLay3(g1, g2, g3, "shms_HM_residuals.png", "SHMS MF DC Residuals", "SHMS DC Planes", false, outROOT);
		delete g1; delete g2; delete g3;
	}


	//-----------------------------------------------------------------------------
	//-------------------------------- HMS Light SRC ------------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	TFile* hms_LS_D2_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/hms_Calib_histos17134.root", "READ");
	TFile* hms_LS_Be9_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/hms_Calib_histos17106.root", "READ");
	TFile* hms_LS_B10_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/hms_Calib_histos17114.root", "READ");
	TFile* hms_LS_B11_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/hms_Calib_histos17120.root", "READ");
	TFile* hms_LS_C12_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/hms_Calib_histos17085.root", "READ");

	//----- 1d -----
	if (hms_light_src_1d)
	{
		//Beta
		hist1 = (TH1F*)hms_LS_D2_calib->Get("hHod_Beta");
		hist2 = (TH1F*)hms_LS_Be9_calib->Get("hHod_Beta");
		hist3 = (TH1F*)hms_LS_B10_calib->Get("hHod_Beta");
		hist4 = (TH1F*)hms_LS_B11_calib->Get("hHod_Beta");
		hist5 = (TH1F*)hms_LS_C12_calib->Get("hHod_Beta");
		OverLap5(hist1, hist2, hist3, hist4, hist5, 0.8161, 0, 0.8161, yf, "hms_LS_beta.png", "HMS SRC Hodo Beta", "Beta [p/E]", true, outROOT);
		delete hist1; delete hist2; delete hist3; delete hist4; delete hist5;
	}
	
	//----- 2d -----
	if (hms_light_src_2d)
	{
		//Residuals
		g1 = (TGraph*)hms_LS_D2_calib->Get("Graph;2");
		g2 = (TGraph*)hms_LS_Be9_calib->Get("Graph;2");
		g3 = (TGraph*)hms_LS_B10_calib->Get("Graph;2");
		g4 = (TGraph*)hms_LS_B11_calib->Get("Graph;2");
		g5 = (TGraph*)hms_LS_C12_calib->Get("Graph;2");
		OverLay5(g1, g2, g3, g4, g5, "hms_LS_residuals.png", "HMS SRC DC Residuals", "HMS DC Planes", true, outROOT);
		delete g1; delete g2; delete g3; delete g4; delete g5;
	}



	//-----------------------------------------------------------------------------
	//------------------------------- SHMS Light SRC ------------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	TFile* shms_LS_D2_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/shms_Calib_histos17134.root", "READ");
	TFile* shms_LS_Be9_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/shms_Calib_histos17106.root", "READ");
	TFile* shms_LS_B10_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/shms_Calib_histos17114.root", "READ");
	TFile* shms_LS_B11_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/shms_Calib_histos17120.root", "READ");
	TFile* shms_LS_C12_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/light_SRC/shms_Calib_histos17085.root", "READ");

	//----- 1d -----
	if (shms_light_src_1d)
	{
		//Beta
		hist1 = (TH1F*)shms_LS_D2_calib->Get("pHod_Beta");
		hist2 = (TH1F*)shms_LS_Be9_calib->Get("pHod_Beta");
		hist3 = (TH1F*)shms_LS_B10_calib->Get("pHod_Beta");
		hist4 = (TH1F*)shms_LS_B11_calib->Get("pHod_Beta");
		hist5 = (TH1F*)shms_LS_C12_calib->Get("pHod_Beta");
		OverLap5(hist1, hist2, hist3, hist4, hist5, 1, 0, 1, yf, "shms_LS_beta.png", "SHMS SRC Hodo Beta", "Beta [p/E]", false, outROOT);
		delete hist1; delete hist2; delete hist3; delete hist4; delete hist5;
		
		//etotTrkNorm
		hist1 = (TH1F*)shms_LS_D2_calib->Get("pCal_eTrkNorm");	hist1->Scale(1. / hist1->Integral(), "width");
		hist2 = (TH1F*)shms_LS_Be9_calib->Get("pCal_eTrkNorm");	hist2->Scale(1. / hist2->Integral(), "width");
		hist3 = (TH1F*)shms_LS_B10_calib->Get("pCal_eTrkNorm");	hist3->Scale(1. / hist3->Integral(), "width");
		hist4 = (TH1F*)shms_LS_B11_calib->Get("pCal_eTrkNorm");	hist4->Scale(1. / hist4->Integral(), "width");
		hist5 = (TH1F*)shms_LS_C12_calib->Get("pCal_eTrkNorm");	hist5->Scale(1. / hist5->Integral(), "width");
		OverLap5(hist1, hist2, hist3, hist4, hist5, 1, 0, 1, yf, "shms_LS_eTrkNorm.png", "SHMS SRC Cal eTrkNorm", "eTrkNorm [E_{dep}/p_{track}]", false, outROOT);
		delete hist1; delete hist2; delete hist3; delete hist4; delete hist5;
	}
	
	//----- 2d -----
	if (shms_light_src_2d)
	{
		//Residuals
		g1 = (TGraph*)shms_LS_D2_calib->Get("Graph;2");
		g2 = (TGraph*)shms_LS_Be9_calib->Get("Graph;2");
		g3 = (TGraph*)shms_LS_B10_calib->Get("Graph;2");
		g4 = (TGraph*)shms_LS_B11_calib->Get("Graph;2");
		g5 = (TGraph*)shms_LS_C12_calib->Get("Graph;2");
		OverLay5(g1, g2, g3, g4, g5, "shms_LS_residuals.png", "SHMS SRC DC Residuals", "SHMS DC Planes", false, outROOT);
		delete g1; delete g2; delete g3; delete g4; delete g5;
	}
	
	
	
	//-----------------------------------------------------------------------------
	//-------------------------------- HMS Heavy SRC ------------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	TFile* hms_HS_Ca40_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_SRC/hms_Calib_histos17032.root", "READ");
	TFile* hms_HS_Ca48_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_SRC/hms_Calib_histos17056.root", "READ");
	TFile* hms_HS_Fe54_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_SRC/hms_Calib_histos17144.root", "READ");

	//----- 1d -----
	if (hms_heavy_src_1d)
	{
		//Beta
		hist1 = (TH1F*)hms_HS_Ca40_calib->Get("hHod_Beta");
		hist2 = (TH1F*)hms_HS_Ca48_calib->Get("hHod_Beta");
		hist3 = (TH1F*)hms_HS_Fe54_calib->Get("hHod_Beta");
		OverLap3(hist1, hist2, hist3, 0.8161, 0, 0.8161, yf, "hms_HS_beta.png", "HMS SRC Hodo Beta", "Beta [p/E]", true, outROOT);
		delete hist1; delete hist2; delete hist3;
	}
	
	//----- 2d -----
	if (hms_heavy_src_2d)
	{
		//Residuals
		g1 = (TGraph*)hms_HS_Ca40_calib->Get("Graph;2");
		g2 = (TGraph*)hms_HS_Ca48_calib->Get("Graph;2");
		g3 = (TGraph*)hms_HS_Fe54_calib->Get("Graph;2");
		OverLay3(g1, g2, g3, "hms_HS_residuals.png", "HMS SRC DC Residuals", "HMS DC Planes", true, outROOT);
		delete g1; delete g2; delete g3;
	}



	//-----------------------------------------------------------------------------
	//------------------------------- SHMS Heavy SRC ------------------------------
	//-----------------------------------------------------------------------------

	//----- Files -----
	TFile* shms_HS_Ca40_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_SRC/shms_Calib_histos17032.root", "READ");
	TFile* shms_HS_Ca48_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_SRC/shms_Calib_histos17056.root", "READ");
	TFile* shms_HS_Fe54_calib = new TFile("../../../../OFFLINE/CALIBRATIONS/cafe_checks/heavy_SRC/shms_Calib_histos17144.root", "READ");

	//----- ld -----
	if (shms_heavy_src_1d)
	{
		//Beta
		hist1 = (TH1F*)shms_HS_Ca40_calib->Get("pHod_Beta");
		hist2 = (TH1F*)shms_HS_Ca48_calib->Get("pHod_Beta");
		hist3 = (TH1F*)shms_HS_Fe54_calib->Get("pHod_Beta");
		OverLap3(hist1, hist2, hist3, 1, 0, 1, yf, "shms_HS_beta.png", "SHMS SRC Hodo Beta", "Beta [p/E]", true, outROOT);
		delete hist1; delete hist2; delete hist3;

		//eTrkNorm
		hist1 = (TH1F*)shms_HS_Ca40_calib->Get("pCal_eTrkNorm");
		hist2 = (TH1F*)shms_HS_Ca48_calib->Get("pCal_eTrkNorm");
		hist3 = (TH1F*)shms_HS_Fe54_calib->Get("pCal_eTrkNorm");
		OverLap3(hist1, hist2, hist3, 1, 0, 1, yf, "shms_HS_eTrkNorm.png", "SHMS SRC Cal eTrkNorm", "eTrkNorm [E_{dep}/p_{track}]", true, outROOT);
		delete hist1; delete hist2; delete hist3;
	}
	
	//----- 2d -----
	if (shms_heavy_src_2d)
	{
		//Residuals
		g1 = (TGraph*)shms_HS_Ca40_calib->Get("Graph;2");
		g2 = (TGraph*)shms_HS_Ca48_calib->Get("Graph;2");
		g3 = (TGraph*)shms_HS_Fe54_calib->Get("Graph;2");
		OverLay3(g1, g2, g3, "shms_HS_residuals.png", "SHMS SRC DC Residuals", "SHMS DC Planes", false, outROOT);
		delete g1; delete g2; delete g3;
	}


	outROOT->Close();	//Close File
}
//*****************************************************************************************************************************************************************************
//End Main
//*****************************************************************************************************************************************************************************




void Draw2d(TH2F* hist1, const char* name, const char* title, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");
	
	
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();
	
	hist1->Draw("Col");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void HotWire(TH2F* hist1, const char* name, const char* title, double x1i, double x1f, double x2i, double x2f, const char* title2, const char* title3, TFile* outROOT)
{
	
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.10);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");
	
	
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->Divide(2,2);

	c->cd(1);
	hist1->Draw("Col");

	c->cd(3);
	TH1D* hist2 = hist1->ProjectionY("Normal", x1i, x1f, "");
	hist2->SetTitle(title3); 
	
	TH1D* hist3 = hist1->ProjectionY("Anomalous", x2i, x2f, "");
	hist3->SetTitle(title2); 

	double yf = 0;
	if (hist2->GetMaximum() > yf) { yf = hist2->GetMaximum(); }
	if (hist3->GetMaximum() > yf) { yf = hist3->GetMaximum(); }
	
	hist2->GetYaxis()->SetRangeUser(0, yf);
	hist3->GetYaxis()->SetRangeUser(0, yf);
	yf = 0;

	hist2->Draw("SAME");

	c->cd(4);
	hist3->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}


void OverLap2(TH1F* hist1, TH1F* hist2, double xi, double yi, double xf, double yf, const char* name, const char* title, bool line, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->SetStats(0);
	hist1->SetLineWidth(2);
	hist1->SetLineColor(kRed);
	hist1->SetFillColorAlpha(kRed, 0.40);
	hist1->SetFillStyle(3004);
	hist1->Scale(1. / hist1->Integral(), "width");

	hist2->SetStats(0);
	hist2->SetLineWidth(2);
	hist2->SetLineColor(kBlue);
	hist2->SetFillColorAlpha(kBlue, 0.40);
	hist2->SetFillStyle(3004);
	hist2->Scale(1. / hist2->Integral(), "width");




	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");



	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = hist1->GetMaximum();
	if (yaxisrange < hist2->GetMaximum()) { yaxisrange = hist2->GetMaximum(); }
	hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.3 * yaxisrange));
	 
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	hist1->Draw("HIST");
	hist2->Draw("HIST SAME");

	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, "Calibrated", "lp"); legend->AddEntry(hist2, "Uncalibrated", "lp"); legend->Draw("SAME");

	if (line)
	{
		yf = hist1->GetMaximum();
		if (hist2->GetMaximum() > yf) { yf = hist2->GetMaximum(); }
		TLine* line = new TLine(xi, yi, xf, yf); line->SetLineColor(1); line->SetLineWidth(2); line->Draw("SAME");
	}
	
	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void Residuals(TH1F* hist1, TH1F* hist2, double xi, double yi, double xf, double yf, const char* name, const char* title, bool line, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->SetStats(0);
	hist1->SetLineWidth(2);
	hist1->SetLineColor(kRed);
	hist1->SetFillColorAlpha(kRed, 0.40);
	hist1->SetFillStyle(3004);
	hist1->Scale(1. / hist1->Integral(), "width");

	hist2->SetStats(0);
	hist2->SetLineWidth(2);
	hist2->SetLineColor(kBlue);
	hist2->SetFillColorAlpha(kBlue, 0.40);
	hist2->SetFillStyle(3005);
	hist2->Scale(1. / hist2->Integral(), "width");

	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");




	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = hist1->GetMaximum();
	if (yaxisrange < hist2->GetMaximum()) { yaxisrange = hist2->GetMaximum(); }
	hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.3 * yaxisrange));

	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	hist1->Draw("HIST");
	hist2->Draw("HIST SAME");

	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, "Calibrated", "lp"); legend->AddEntry(hist2, "Uncalibrated", "lp"); legend->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void OverLap3(TH1F* hist1, TH1F* hist2, TH1F* hist3, double xi, double yi, double xf, double yf, const char* name, const char* title, const char* xaxis, bool line, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);


	


	hist1->SetLineWidth(2);
	hist2->SetLineWidth(2);
	hist3->SetLineWidth(2);


	hist1->SetLineColor(kRed);
	hist1->SetFillColorAlpha(kRed, 0.40);
	hist1->SetFillStyle(3004);
	hist1->SetTitle(title);
	hist1->GetXaxis()->SetTitle(xaxis);
	hist1->Scale(1. / hist1->Integral(), "width");


	hist2->SetLineColor(kGreen);
	hist2->SetFillColorAlpha(kGreen, 0.40);
	hist1->SetFillStyle(3005);
	//hist2->SetTitle(title);
	//hist2->GetXaxis()->SetTitle(xaxis);
	hist2->Scale(1. / hist2->Integral(), "width");


	hist3->SetLineColor(kBlue);
	hist3->SetFillColorAlpha(kBlue, 0.40);
	hist1->SetFillStyle(3006);
	//hist3->SetTitle(title);
	//hist3->GetXaxis()->SetTitle(xaxis);
	hist3->Scale(1. / hist3->Integral(), "width");


	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");


	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();
	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = hist1->GetMaximum();
	if (yaxisrange < hist2->GetMaximum()) { yaxisrange = hist2->GetMaximum(); }
	if (yaxisrange < hist3->GetMaximum()) { yaxisrange = hist3->GetMaximum(); }
	hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.25 * yaxisrange));

	hist1->Draw("histE0");
	hist2->Draw("samehistE0");
	hist3->Draw("samehistE0");



	


	

	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, "Ca40", "lp"); legend->AddEntry(hist2, "Ca48", "lp"); legend->AddEntry(hist3, "Fe54", "lp"); legend->Draw("SAME");


	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void OverLap5(TH1F* hist1, TH1F* hist2, TH1F* hist3, TH1F* hist4, TH1F* hist5, double xi, double yi, double xf, double yf, const char* name, const char* title, const char* xaxis, bool mf, TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);

	hist1->SetStats(0);
	hist1->SetLineWidth(2);
	hist1->SetLineColor(kRed);
	hist1->SetFillColorAlpha(kRed, 0.40);
	hist1->SetFillStyle(3004);
	hist1->SetTitle(title);
	hist1->GetXaxis()->SetTitle(xaxis);
	hist1->Scale(1. / hist1->Integral(), "width");

	hist2->SetStats(0);
	hist2->SetLineWidth(2);
	hist2->SetLineColor(kMagenta);
	hist2->SetFillColorAlpha(kMagenta, 0.40);
	hist1->SetFillStyle(3005);
	hist2->SetTitle(title);
	hist2->GetXaxis()->SetTitle(xaxis);
	hist2->Scale(1. / hist2->Integral(), "width");

	hist3->SetStats(0);
	hist3->SetLineWidth(2);
	hist3->SetLineColor(kGreen);
	hist3->SetFillColorAlpha(kGreen, 0.40);
	hist1->SetFillStyle(3006);
	hist3->SetTitle(title);
	hist3->GetXaxis()->SetTitle(xaxis);
	hist3->Scale(1. / hist3->Integral(), "width");

	hist4->SetStats(0);
	hist4->SetLineWidth(2);
	hist4->SetLineColor(kCyan);
	hist4->SetFillColorAlpha(kCyan, 0.40);
	hist1->SetFillStyle(3007);
	hist4->SetTitle(title);
	hist4->GetXaxis()->SetTitle(xaxis);
	hist4->Scale(1. / hist4->Integral(), "width");

	hist5->SetStats(0);
	hist5->SetLineWidth(2);
	hist5->SetLineColor(kBlue);
	hist5->SetFillColorAlpha(kBlue, 0.40);
	hist1->SetFillStyle(3008);
	hist5->SetTitle(title);
	hist5->GetXaxis()->SetTitle(xaxis);
	hist5->Scale(1. / hist5->Integral(), "width");

	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = hist1->GetMaximum();
	if (yaxisrange < hist2->GetMaximum()) { yaxisrange = hist2->GetMaximum(); }
	if (yaxisrange < hist3->GetMaximum()) { yaxisrange = hist3->GetMaximum(); }
	if (yaxisrange < hist4->GetMaximum()) { yaxisrange = hist4->GetMaximum(); }
	if (yaxisrange < hist5->GetMaximum()) { yaxisrange = hist5->GetMaximum(); }
	hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.25 * yaxisrange));
	
	// set histogram titles/labels/font
	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");

	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();
	
	if (mf)	{hist1->Draw("HIST"); hist2->Draw("HIST SAME"); hist3->Draw("HIST SAME"); hist4->Draw("HIST SAME"); hist5->Draw("HIST SAME");}
	else	{hist3->Draw("HIST"); hist1->Draw("HIST SAME"); hist2->Draw("HIST SAME"); hist4->Draw("HIST SAME"); hist5->Draw("HIST SAME");}

	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, "D2", "lp"); legend->AddEntry(hist2, "Be9", "lp"); legend->AddEntry(hist3, "B10", "lp"); legend->AddEntry(hist4, "B11", "lp"); legend->AddEntry(hist5, "C12", "lp"); legend->Draw("SAME");

	TLine* line = new TLine(xi, yi, xf, yf); line->SetLineColor(1); line->SetLineWidth(2); line->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);
	
	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void Compare2(TH2F* hist1, TH2F* hist2, double xi, double yi, double xf, double yf, const char* name, const char* title, bool line, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->Divide(2, 1);

	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);
	
	hist1->GetXaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->SetLabelSize(0.04);
	hist1->GetYaxis()->CenterTitle();
	hist1->GetXaxis()->CenterTitle();
	hist1->SetLabelFont(font_type, "XY");
	hist1->SetTitleFont(font_type, "XY");
	hist1->SetTitleSize(0.05, "XY");
	hist1->SetTitleOffset(1., "XY");

	hist2->GetXaxis()->SetLabelSize(0.04);
	hist2->GetYaxis()->SetLabelSize(0.04);
	hist2->GetYaxis()->CenterTitle();
	hist2->GetXaxis()->CenterTitle();
	hist2->SetLabelFont(font_type, "XY");
	hist2->SetTitleFont(font_type, "XY");
	hist2->SetTitleSize(0.05, "XY");
	hist2->SetTitleOffset(1., "XY");

	if (line) {
		c->cd(2); hist1->Draw("Col"); TLine* line1 = new TLine(xi, yi, xf, yf); line1->SetLineColor(2); line1->SetLineWidth(2); line1->Draw("SAME");
		c->cd(1); hist2->Draw("Col"); TLine* line2 = new TLine(xi, yi, xf, yf); line2->SetLineColor(2); line2->SetLineWidth(2); line2->Draw("SAME");
	}
	else {
		c->cd(2); hist1->Draw("Col");
		c->cd(1); hist2->Draw("Col");
	}
	
	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void OverLay3(TGraph* hist1, TGraph* hist2, TGraph* hist3, const char* name, const char* title, const char* xaxis, bool hms, TFile* outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	if (hms) { hist1->GetYaxis()->SetRangeUser(250, 450); hist2->GetYaxis()->SetRangeUser(250, 450); hist3->GetYaxis()->SetRangeUser(250, 450); }
	else	 {hist1->GetYaxis()->SetRangeUser(200, 350); hist2->GetYaxis()->SetRangeUser(200, 350); hist3->GetYaxis()->SetRangeUser(200, 350);}
	
	//hist1->GetXaxis()->SetBinLabel(7, "1u1"); hist1->GetXaxis()->SetBinLabel(14, "1u2"); hist1->GetXaxis()->SetBinLabel(21, "2u1"); hist1->GetXaxis()->SetBinLabel(28, "2u2");
	//hist1->GetXaxis()->SetBinLabel(35, "1x1"); hist1->GetXaxis()->SetBinLabel(42, "1x2"); hist1->GetXaxis()->SetBinLabel(49, "2x1"); hist1->GetXaxis()->SetBinLabel(56, "2x2");
	//hist1->GetXaxis()->SetBinLabel(63, "1v1"); hist1->GetXaxis()->SetBinLabel(70, "1v2"); hist1->GetXaxis()->SetBinLabel(77, "2v1"); hist1->GetXaxis()->SetBinLabel(84, "2v2");

	hist1->SetMarkerColor(2); hist1->SetMarkerStyle(20); hist1->SetMarkerSize(2); hist1->SetLineColor(2); hist1->SetLineWidth(2); hist1->SetTitle(title); hist1->GetXaxis()->SetTitle(xaxis); hist1->Draw("");
	hist2->SetMarkerColor(8); hist2->SetMarkerStyle(21); hist2->SetMarkerSize(2); hist2->SetLineColor(8); hist2->SetLineWidth(2); hist2->SetTitle(title); hist2->GetXaxis()->SetTitle(xaxis); hist2->Draw("PL SAME");
	hist3->SetMarkerColor(4); hist3->SetMarkerStyle(22); hist3->SetMarkerSize(2); hist3->SetLineColor(4); hist3->SetLineWidth(2); hist3->SetTitle(title); hist3->GetXaxis()->SetTitle(xaxis); hist3->Draw("PL SAME");

	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box
	legend->AddEntry(hist1, "Ca40", "lp"); legend->AddEntry(hist2, "Ca48", "lp"); legend->AddEntry(hist3, "Fe54", "lp"); legend->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}




void OverLay5(TGraph* hist1, TGraph* hist2, TGraph* hist3, TGraph* hist4, TGraph* hist5, const char* name, const char* title, const char* xaxis, bool mf, TFile * outROOT)
{
	TCanvas* c = new TCanvas(name, title, 1920, 1080); c->cd();

	if (mf) {hist1->GetYaxis()->SetRangeUser(250, 450); hist2->GetYaxis()->SetRangeUser(250, 450); hist3->GetYaxis()->SetRangeUser(250, 450); hist4->GetYaxis()->SetRangeUser(250, 450); hist5->GetYaxis()->SetRangeUser(250, 450);}
	else	{hist1->GetYaxis()->SetRangeUser(200, 350); hist2->GetYaxis()->SetRangeUser(200, 350); hist3->GetYaxis()->SetRangeUser(200, 350); hist4->GetYaxis()->SetRangeUser(200, 350); hist5->GetYaxis()->SetRangeUser(200, 350);}

	//hist1->GetXaxis()->SetBinLabel(10, "1u1"); hist1->GetXaxis()->SetBinLabel(20, "1u2"); hist1->GetXaxis()->SetBinLabel(30, "2u1"); hist1->GetXaxis()->SetBinLabel(40, "2u2");
	//hist1->GetXaxis()->SetBinLabel(50, "1x1"); hist1->GetXaxis()->SetBinLabel(60, "1x2"); hist1->GetXaxis()->SetBinLabel(70, "2x1"); hist1->GetXaxis()->SetBinLabel(80, "2x2");
	//hist1->GetXaxis()->SetBinLabel(90, "1v1"); hist1->GetXaxis()->SetBinLabel(100, "1v2"); hist1->GetXaxis()->SetBinLabel(110, "2v1"); hist1->GetXaxis()->SetBinLabel(120, "2v2");

	hist1->SetMarkerColor(2); hist1->SetMarkerStyle(20); hist1->SetMarkerSize(2); hist1->SetLineColor(2); hist1->SetLineWidth(2); hist1->SetTitle(title); hist1->GetXaxis()->SetTitle(xaxis); hist1->Draw("");
	hist2->SetMarkerColor(6); hist2->SetMarkerStyle(21); hist2->SetMarkerSize(2); hist2->SetLineColor(6); hist2->SetLineWidth(2); hist2->SetTitle(title); hist2->GetXaxis()->SetTitle(xaxis); hist2->Draw("PL SAME");
	hist3->SetMarkerColor(8); hist3->SetMarkerStyle(22); hist3->SetMarkerSize(2); hist3->SetLineColor(8); hist3->SetLineWidth(2); hist3->SetTitle(title); hist3->GetXaxis()->SetTitle(xaxis); hist3->Draw("PL SAME");
	hist4->SetMarkerColor(7); hist4->SetMarkerStyle(23); hist4->SetMarkerSize(2); hist4->SetLineColor(7); hist4->SetLineWidth(2); hist4->SetTitle(title); hist4->GetXaxis()->SetTitle(xaxis); hist4->Draw("PL SAME");
	hist5->SetMarkerColor(4); hist5->SetMarkerStyle(29); hist5->SetMarkerSize(2); hist5->SetLineColor(4); hist5->SetLineWidth(2); hist5->SetTitle(title); hist5->GetXaxis()->SetTitle(xaxis); hist5->Draw("PL SAME");
	
	TLegend* legend = new TLegend(0.16, 0.68, 0.24, 0.88);//x1, y1, x2, y2, lower left and upper right corner of legends box,
	legend->AddEntry(hist1, "D2", "lp"); legend->AddEntry(hist2, "Be9", "lp"); legend->AddEntry(hist3, "B10", "lp"); legend->AddEntry(hist4, "B11", "lp"); legend->AddEntry(hist5, "C12", "lp"); legend->Draw("SAME");

	outROOT->cd(); c->Write(); c->Print(name);
	
	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}

























void compare_histos_light(
	TString file1_path, TString hist1,
	TString file2_path, TString hist2,
	TString file3_path, TString hist3,
	TString file4_path, TString hist4,
	TString file5_path, TString hist5,
	TString xlabel, TString ylabel, TString title,
	TString hist1_leg, TString hist2_leg, TString hist3_leg, TString hist4_leg, TString hist5_leg,
	bool norm,
	TFile* outROOT)
{

	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);


	//Open  ROOT files;
	TFile* file1 = NULL;
	TFile* file2 = NULL;
	TFile* file3 = NULL;
	TFile* file4 = NULL;
	TFile* file5 = NULL;

	file1 = new TFile(file1_path.Data());
	file2 = new TFile(file2_path.Data());
	file3 = new TFile(file3_path.Data());
	file4 = new TFile(file4_path.Data());
	file5 = new TFile(file5_path.Data());

	// declare 1D histos
	TH1F* H_hist1 = 0;
	TH1F* H_hist2 = 0;
	TH1F* H_hist3 = 0;
	TH1F* H_hist4 = 0;
	TH1F* H_hist5 = 0;

	// get histogram objects
	file1->cd();
	file1->GetObject(hist1.Data(), H_hist1);
	file2->cd();
	file2->GetObject(hist2.Data(), H_hist2);
	file3->cd();
	file3->GetObject(hist3.Data(), H_hist3);
	file4->cd();
	file4->GetObject(hist4.Data(), H_hist4);
	file5->cd();
	file5->GetObject(hist5.Data(), H_hist5);

	double h1_I, h2_I, h3_I, h4_I, h5_I;
	double h1_Ierr, h2_Ierr, h3_Ierr, h4_Ierr, h5_Ierr;
	double nbins = H_hist1->GetNbinsX();  //Get total number of bins (excluding overflow)
	h1_I = H_hist1->IntegralAndError(1, nbins, h1_Ierr);
	h2_I = H_hist2->IntegralAndError(1, nbins, h2_Ierr);
	h3_I = H_hist3->IntegralAndError(1, nbins, h3_Ierr);
	h4_I = H_hist4->IntegralAndError(1, nbins, h4_Ierr);
	h5_I = H_hist5->IntegralAndError(1, nbins, h5_Ierr);

	H_hist1->SetLineWidth(2);
	H_hist2->SetLineWidth(2);
	H_hist3->SetLineWidth(2);
	H_hist4->SetLineWidth(2);
	H_hist5->SetLineWidth(2);

	// set histos aethetics
	H_hist1->SetLineColor(kRed);
	H_hist1->SetFillColorAlpha(kRed, 0.40);
	H_hist1->SetFillStyle(3004);
	H_hist1->Scale(1. / H_hist1->Integral(), "width");

	H_hist2->SetLineColor(kMagenta);
	H_hist2->SetFillColorAlpha(kMagenta, 0.40);
	H_hist2->SetFillStyle(3005);
	H_hist2->Scale(1. / H_hist2->Integral(), "width");

	H_hist3->SetLineColor(kGreen);
	H_hist3->SetFillColorAlpha(kGreen, 0.40);
	H_hist3->SetFillStyle(3006);
	H_hist3->Scale(1. / H_hist3->Integral(), "width");

	H_hist4->SetLineColor(kCyan);
	H_hist4->SetFillColorAlpha(kCyan, 0.40);
	H_hist4->SetFillStyle(3007);
	H_hist4->Scale(1. / H_hist4->Integral(), "width");

	H_hist5->SetLineColor(kBlue);
	H_hist5->SetFillColorAlpha(kBlue, 0.40);
	H_hist5->SetFillStyle(3008);
	H_hist5->Scale(1. / H_hist5->Integral(), "width");

	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = H_hist1->GetMaximum();
	if (yaxisrange < H_hist2->GetMaximum()) { yaxisrange = H_hist2->GetMaximum(); }
	if (yaxisrange < H_hist3->GetMaximum()) { yaxisrange = H_hist3->GetMaximum(); }
	if (yaxisrange < H_hist4->GetMaximum()) { yaxisrange = H_hist4->GetMaximum(); }
	if (yaxisrange < H_hist5->GetMaximum()) { yaxisrange = H_hist5->GetMaximum(); }
	H_hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.25 * yaxisrange));


	// set histogram titles/labels/font
	H_hist1->SetTitle(title);

	H_hist1->GetXaxis()->SetLabelSize(0.04);
	H_hist1->GetYaxis()->SetLabelSize(0.04);

	H_hist1->GetYaxis()->SetTitle(ylabel);
	H_hist1->GetXaxis()->SetTitle(xlabel);

	H_hist1->GetYaxis()->CenterTitle();
	H_hist1->GetXaxis()->CenterTitle();

	H_hist1->SetLabelFont(font_type, "XY");
	H_hist1->SetTitleFont(font_type, "XY");
	H_hist1->SetTitleSize(0.05, "XY");
	H_hist1->SetTitleOffset(1., "XY");


	TCanvas* c = new TCanvas("c", "c", 1920, 1080);

	H_hist1->Draw("histE0");
	H_hist2->Draw("sameshistE0");
	H_hist3->Draw("sameshistE0");
	H_hist4->Draw("sameshistE0");
	H_hist5->Draw("sameshistE0");

	// create legend ( displays hist legend label and integral counts)
	TLegend* leg = new TLegend(0.14, 0.89, 0.25, 0.78);
	leg->AddEntry(H_hist1, Form("%s | Integral: %.3f", hist1_leg.Data(), h1_I), "f");
	leg->AddEntry(H_hist2, Form("%s | Integral: %.3f", hist2_leg.Data(), h2_I));
	leg->AddEntry(H_hist3, Form("%s | Integral: %.3f", hist3_leg.Data(), h3_I));
	leg->AddEntry(H_hist4, Form("%s | Integral: %.3f", hist4_leg.Data(), h4_I));
	leg->AddEntry(H_hist5, Form("%s | Integral: %.3f", hist5_leg.Data(), h5_I));
	// draw legend
	leg->Draw();

	outROOT->cd(); c->Write(); c->Print(Form("%s.png", title.Data()));

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}

















void compare_histos_heavy(
	TString file1_path, TString hist1,
	TString file2_path, TString hist2,
	TString file3_path, TString hist3,
	TString xlabel, TString ylabel, TString title,
	TString hist1_leg, TString hist2_leg, TString hist3_leg,
	bool norm,
	TFile* outROOT)
{
	int font_type = 132;

	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetLabelSize(0.05);//
	gStyle->SetTitleFont(font_type, "");
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFont(font_type);
	gStyle->SetLegendTextSize(0.03);


	//Open  ROOT files;
	TFile* file1 = NULL;
	TFile* file2 = NULL;
	TFile* file3 = NULL;

	file1 = new TFile(file1_path.Data());
	file2 = new TFile(file2_path.Data());
	file3 = new TFile(file3_path.Data());

	// declare 1D histos
	TH1F* H_hist1 = 0;
	TH1F* H_hist2 = 0;
	TH1F* H_hist3 = 0;

	// get histogram objects
	file1->cd();
	file1->GetObject(hist1.Data(), H_hist1);
	file2->cd();
	file2->GetObject(hist2.Data(), H_hist2);
	file3->cd();
	file3->GetObject(hist3.Data(), H_hist3);

	double h1_I, h2_I, h3_I;
	double h1_Ierr, h2_Ierr, h3_Ierr;
	double nbins = H_hist1->GetNbinsX();  //Get total number of bins (excluding overflow)
	h1_I = H_hist1->IntegralAndError(1, nbins, h1_Ierr);
	h2_I = H_hist2->IntegralAndError(1, nbins, h2_Ierr);
	h3_I = H_hist3->IntegralAndError(1, nbins, h3_Ierr);

	H_hist1->SetLineWidth(2);
	H_hist2->SetLineWidth(2);
	H_hist3->SetLineWidth(2);

	// set histos aethetics
	H_hist1->SetLineColor(kRed);
	H_hist1->SetFillColorAlpha(kRed, 0.40);
	H_hist1->SetFillStyle(3004);
	H_hist1->Scale(1. / H_hist1->Integral(), "width");

	H_hist2->SetLineColor(kGreen);
	H_hist2->SetFillColorAlpha(kGreen, 0.40);
	H_hist2->SetFillStyle(3005);
	H_hist2->Scale(1. / H_hist2->Integral(), "width");

	H_hist3->SetLineColor(kBlue);
	H_hist3->SetFillColorAlpha(kBlue, 0.40);
	H_hist3->SetFillStyle(3006);
	H_hist3->Scale(1. / H_hist3->Integral(), "width");

	// set y-range
	Double_t yaxisrange = 0;
	yaxisrange = H_hist1->GetMaximum();
	if (yaxisrange < H_hist2->GetMaximum()) { yaxisrange = H_hist2->GetMaximum(); }
	if (yaxisrange < H_hist3->GetMaximum()) { yaxisrange = H_hist3->GetMaximum(); }
	H_hist1->GetYaxis()->SetRangeUser(0, yaxisrange + 0.25 * yaxisrange);

	// set histogram titles/labels/font
	H_hist1->SetTitle(title);

	H_hist1->GetXaxis()->SetLabelSize(0.04);
	H_hist1->GetYaxis()->SetLabelSize(0.04);

	H_hist1->GetYaxis()->SetTitle(ylabel);
	H_hist1->GetXaxis()->SetTitle(xlabel);

	H_hist1->GetYaxis()->CenterTitle();
	H_hist1->GetXaxis()->CenterTitle();

	H_hist1->SetLabelFont(font_type, "XY");
	H_hist1->SetTitleFont(font_type, "XY");
	H_hist1->SetTitleSize(0.05, "XY");
	H_hist1->SetTitleOffset(1., "XY");


	TCanvas* c = new TCanvas("c", "c", 1920, 1080);

	H_hist1->Draw("histE0");
	H_hist2->Draw("sameshistE0");
	H_hist3->Draw("sameshistE0");

	// create legend ( displays hist legend label and integral counts)
	TLegend* leg = new TLegend(0.14, 0.89, 0.25, 0.78);
	leg->AddEntry(H_hist1, Form("%s | Integral: %.3f", hist1_leg.Data(), h1_I), "f");
	leg->AddEntry(H_hist2, Form("%s | Integral: %.3f", hist2_leg.Data(), h2_I));
	leg->AddEntry(H_hist3, Form("%s | Integral: %.3f", hist3_leg.Data(), h3_I));
	leg->Draw();

	outROOT->cd(); c->Write(); c->Print(Form("%s.png", title.Data()));

	if (c) { c->Close(); gSystem->ProcessEvents(); delete c; c = 0; }//delete canvas
}