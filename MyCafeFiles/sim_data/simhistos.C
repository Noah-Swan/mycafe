
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
void skimhistos() 
{
	TFile* outROOT = new TFile("outHist.root", "RECREATE");
	
	const char* Target[8] = { "D2", "Be9", "B10", "B11", "C12", "Ca40", "Ca48", "Fe54" };
	const char* MF_File[8] = { "../MF/Mf_Skim_D2/MF_D2_Results.root", "../MF/Mf_Skim_Be9/MF_Be9_Results.root" , "../MF/Mf_Skim_B10/MF_B10_Results.root" , "../MF/Mf_Skim_B11/MF_B11_Results.root" , "../MF/Mf_Skim_C12/MF_C12_Results.root" , "../MF/Mf_Skim_Ca40/MF_Ca40_Results.root" , "../MF/Mf_Skim_Ca48/MF_Ca48_Results.root" , "../MF/Mf_Skim_Fe54/MF_Fe54_Results.root" };
	const char* SRC_File[8] = { "../SRC/Src_Skim_D2/SRC_D2_Results.root", "../SRC/Src_Skim_Be9/SRC_Be9_Results.root" , "../SRC/Src_Skim_B10/SRC_B10_Results.root" , "../SRC/Src_Skim_B11/SRC_B11_Results.root" , "../SRC/Src_Skim_C12/SRC_C12_Results.root" , "../SRC/Src_Skim_Ca40/SRC_Ca40_Results.root" , "../SRC/Src_Skim_Ca48/SRC_Ca48_Results.root" , "../SRC/Src_Skim_Fe54/SRC_Fe54_Results.root" };


	bool Prim_Kin = true;
	bool Sec_Kin = true;
	bool HMS_Accp = true;
	bool SHMS_Accp = true;



	const char* Prim_KinHist[11][3] = {	
		{"th_e", "#theta_{e} [deg]", "Electron Scattering Angle"},
		{"W", "W [GeV]", "Invariant Mass"}, 
		{"Q2", "Q^{2} [(GeV/c)^{2}]", "4-Momentum Transfer"}, 
		{"xbj", "x_{B}", "x-Bjorken"}, 
		{"nu", "#nu [GeV/c]", "Energy Transfer"}, 
		{"q", "|#vec{q}| [GeV/c]", "3-Momentum Transfer"}, 
		{"qx", "|#vec{q}_{x}| [GeV/c]", "3-Momentum Transfer (X-Component)"}, 
		{"qy", "|#vec{q}_{y}| [GeV/c]", "3-Momentum Transfer (Y-Component)"}, 
		{"qz", "|#vec{q}_{z}| [GeV/c]", "3-Momentum Transfer (Z-Component)"}, 
		{"th_q", "#theta_{q} [deg]", "In-Plane angle between #vec{q} and beam (+Z)"}, 
		{"ph_q", "#phi_{q} [deg]", "Out-of-Plane angle between #vec{q} and the beam (+Z)"} 
		//{"ph_q", "#phi_{q} [deg]", "Phi of 3-Momentum Vector"},//final e- momentum
		//{"ph_q", "#phi_{q} [deg]", "Phi of 3-Momentum Vector"},//initial e- momentum
	};
	



	const char* Sec_KinHist[18][3] = {
		{"emiss_nuc", "E_{miss, nuc} [GeV]", "Missing Energy (Nuclear Physics)"},
		{"pmiss", "P_{miss} [GeV]", "Missing Momentum"},
		{"prec_x", "P_{miss, x} [GeV]", "P_{miss, x} (Lab)"},
		{"prec_y", "P_{miss, y} [GeV]", "P_{miss, y} (Lab)"},
		{"prec_z", "P_{miss, z} [GeV]", "P_{miss, z} (Lab)"},
		{"pmiss_x", "P_{miss, xq} [GeV]", "P_{miss, xq} (wrt #vec{q})"},
		{"pmiss_y", "P_{miss, yq} [GeV]", "P_{miss, yq} (wrt #vec{q})"},
		{"pmiss_z", "P_{miss, zq} [GeV]", "P_{miss, zq} (wrt #vec{q})"},
		{"Tx", "T_{p} [GeV]", "Kinetic Energy (detected)"},
		{"Tr", "T_{r} [GeV]", "Kinetic Energy (recoil)"},
		{"mmiss", "Mrecoil [GeV/c^{2}]", "Invariant Mass of the Recoil System"},
		{"th_pq", "#theta_{pq} [deg]", "In-Plane (detected) Angle"},
		{"th_rq", "#theta_{rq} [deg]", "In-Plane (recoil) Angle"},
		{"cth_rq", "cos(#theta_{rq})", "Cos of In-Plane (recoil) Angle"},
		{"ph_pq", "#phi_{pq} [deg]", "Out-of-Plane (detected) Angle"},
		{"ph_rq", "#phi_{rq} [deg]", "Out-of-Plane (recoil) Angle"},
		{"xangle", "xangle [deg]", "Angle Between the Detected Particle and Scattered Electron"},
		{"th_p", "#theta_{p} [deg]", "Hadron Scattering Angle (detected)"}
		//{"th_p", "#theta_{p} [deg]", "Hadron Scattering Angle (detected)"} final hadron momentum
	};



	const char* HMS_AccpHist[13][3] = {
		{"hxfp", "X_{fp} [cm]", "HMS X_{fp}"},
		{"hxpfp", "X'_{fp} [deg]", "HMS X'_{fp}"},
		{"hyfp", "Y_{fp} [cm]", "HMS Y_{fp}"},
		{"hypfp", "Y'_{fp} [deg]", "HMS Y'_{fp}"},
		{"hytar", "Y_{tar} [cm]", "HMS Y_{tar}"},
		{"hyptar", "Y'_{tar} [deg]", "HMS Y'_{tar}"},
		{"hxptar", "X'_{tar} [deg]", "HMS X'_{tar}"},
		{"hdelta", "#delta [%]", "HMS Momentum Acceptance"},
		{"htarx", "x-Target [cm]", "HMS x-Target (Lab)"},
		{"htary", "y-Target [cm]", "HMS y-Target (Lab)"},
		{"htarz", "z-Target [cm]", "HMS z-Target (Lab)"},
		{"hXColl", "X-Collimator [cm]", "HMS X-Collimator"},
		{"hYColl", "Y-Collimator [cm]", "HMS Y-Collimator"}
	};




	const char* SHMS_AccpHist[13][3] = {
		{"exfp", "X_{fp} [cm]", "SHMS X_{fp}"},
		{"expfp", "X'_{fp} [deg]", "SHMS X'_{fp}"},
		{"eyfp", "Y_{fp} [cm]", "SHMS Y_{fp}"},
		{"eypfp", "Y'_{fp} [deg]", "SHMS Y'_{fp}"},
		{"eytar", "Y_{tar} [cm]", "SHMS Y_{tar}"},
		{"eyptar", "Y'_{tar} [deg]", "SHMS Y'_{tar}"},
		{"exptar", "X'_{tar} [deg]", "SHMS X'_{tar}"},
		{"edelta", "#delta [%]", "SHMS Momentum Acceptance"},
		{"etarx", "x-Target [cm]", "SHMS x-Target (Lab)"},
		{"etary", "y-Target [cm]", "SHMS y-Target (Lab)"},
		{"etarz", "z-Target [cm]", "SHMS z-Target (Lab)"},
		{"eXColl", "X-Collimator [cm]", "SHMS X-Collimator"},
		{"eYColl", "Y-Collimator [cm]", "SHMS Y-Collimator"}
	};




	if (Prim_Kin)
	{
		for (int i = 0; i < 11; i++)
		{
			compare_histos_light(
				MF_File[0], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				MF_File[1], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				MF_File[2], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				MF_File[3], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				MF_File[4], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				Prim_KinHist[i][1], "", Form("Light MF %s", Prim_KinHist[i][2]),
				Target[0], Target[1], Target[2], Target[3], Target[4], false, outROOT);

			compare_histos_light(
				SRC_File[0], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				SRC_File[1], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				SRC_File[2], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				SRC_File[3], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				SRC_File[4], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				Prim_KinHist[i][1], "", Form("Light SRC %s", Prim_KinHist[i][2]),
				Target[0], Target[1], Target[2], Target[3], Target[4], false, outROOT);

			compare_histos_heavy(
				MF_File[5], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				MF_File[6], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				MF_File[7], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				Prim_KinHist[i][1], "", Form("Heavy MF %s", Prim_KinHist[i][2]),
				Target[5], Target[6], Target[7], false, outROOT);

			compare_histos_heavy(
				SRC_File[5], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				SRC_File[6], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				SRC_File[7], Form("Prim_Kin/H1_%s_pid_acc_kin", Prim_KinHist[i][0]),
				Prim_KinHist[i][1], "", Form("Heavy SRC %s", Prim_KinHist[i][2]),
				Target[5], Target[6], Target[7], false, outROOT);
		}
	}




	if (Sec_Kin)
	{
		for (int i = 0; i < 18; i++)
		{
			compare_histos_light(
				MF_File[0], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				MF_File[1], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				MF_File[2], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				MF_File[3], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				MF_File[4], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				Sec_KinHist[i][1], "", Form("Light MF %s", Sec_KinHist[i][2]),
				Target[0], Target[1], Target[2], Target[3], Target[4], false, outROOT);

			compare_histos_light(
				SRC_File[0], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				SRC_File[1], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				SRC_File[2], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				SRC_File[3], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				SRC_File[4], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				Sec_KinHist[i][1], "", Form("Light SRC %s", Sec_KinHist[i][2]),
				Target[0], Target[1], Target[2], Target[3], Target[4], false, outROOT);

			compare_histos_heavy(
				MF_File[5], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				MF_File[6], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				MF_File[7], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				Sec_KinHist[i][1], "", Form("Heavy MF %s", Sec_KinHist[i][2]),
				Target[5], Target[6], Target[7], false, outROOT);

			compare_histos_heavy(
				SRC_File[5], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				SRC_File[6], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				SRC_File[7], Form("Sec_Kin/H1_%s_pid_acc_kin", Sec_KinHist[i][0]),
				Sec_KinHist[i][1], "", Form("Heavy SRC %s", Sec_KinHist[i][2]),
				Target[5], Target[6], Target[7], false, outROOT);
		}
	}


	

	if (HMS_Accp)
	{
		for (int i = 0; i < 13; i++)
		{
			compare_histos_light(
				MF_File[0], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				MF_File[1], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				MF_File[2], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				MF_File[3], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				MF_File[4], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				HMS_AccpHist[i][1], "", Form("Light MF %s", HMS_AccpHist[i][2]),
				Target[0], Target[1], Target[2], Target[3], Target[4], false, outROOT);

			compare_histos_light(
				SRC_File[0], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				SRC_File[1], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				SRC_File[2], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				SRC_File[3], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				SRC_File[4], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				HMS_AccpHist[i][1], "", Form("Light SRC %s", HMS_AccpHist[i][2]),
				Target [0], Target[1], Target[2], Target[3], Target[4], false, outROOT);

			compare_histos_heavy(
				MF_File[5], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				MF_File[6], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				MF_File[7], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				HMS_AccpHist[i][1], "", Form("Heavy MF %s", HMS_AccpHist[i][2]),
				Target[5], Target[6], Target[7], false, outROOT);

			compare_histos_heavy(
				SRC_File[5], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				SRC_File[6], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				SRC_File[7], Form("HMS_Accp/H1_%s_pid_acc_kin", HMS_AccpHist[i][0]),
				HMS_AccpHist[i][1], "", Form("Heavy SRC %s", HMS_AccpHist[i][2]),
				Target[5], Target[6], Target[7], false, outROOT);
		}
	}

	


	if (SHMS_Accp)
	{
		for (int i = 0; i < 13; i++)
		{
			compare_histos_light(
				MF_File[0], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				MF_File[1], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				MF_File[2], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				MF_File[3], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				MF_File[4], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SHMS_AccpHist[i][1], "", Form("Light MF %s", SHMS_AccpHist[i][2]),
				Target[0], Target[1], Target[2], Target[3], Target[4], false, outROOT);

			compare_histos_light(
				SRC_File[0], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SRC_File[1], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SRC_File[2], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SRC_File[3], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SRC_File[4], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SHMS_AccpHist[i][1], "", Form("Light SRC %s", SHMS_AccpHist[i][2]),
				Target[0], Target[1], Target[2], Target[3], Target[4], false, outROOT);

			compare_histos_heavy(
				MF_File[5], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				MF_File[6], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				MF_File[7], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SHMS_AccpHist[i][1], "", Form("Heavy MF %s", SHMS_AccpHist[i][2]),
				Target[5], Target[6], Target[7], false, outROOT);

			compare_histos_heavy(
				SRC_File[5], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SRC_File[6], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SRC_File[7], Form("SHMS_Accp/H1_%s_pid_acc_kin", SHMS_AccpHist[i][0]),
				SHMS_AccpHist[i][1], "", Form("Heavy SRC %s", SHMS_AccpHist[i][2]),
				Target[5], Target[6], Target[7], false, outROOT);
		}
	}

	





	outROOT->Close();	//Close File
}
//*****************************************************************************************************************************************************************************
//End Main
//*****************************************************************************************************************************************************************************










void compare_histos_light(
			TString file1_path,	TString hist1,
			TString file2_path,	TString hist2,
			TString file3_path,	TString hist3,
			TString file4_path,	TString hist4,
			TString file5_path,  TString hist5,
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
	H_hist1->GetYaxis()->SetRangeUser(0, (yaxisrange + 0.25*yaxisrange));


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