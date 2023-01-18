
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
//Contains public member functions that allow for easy definition of histograms for various useful quantities for the CaFe experiment for both raw data and simulation data.
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

#ifndef Histo_H
#define Histo_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"

#include "Bins.h"
using namespace std;


class Histo
{
	public:
	
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//------------------------------------------------------------------------PID Histograms Bins-----------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//------------------
	//Coincidence (Calc)
	//-------------------------
	static TH1F* mk_epctime(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("CTime.epCoinTime_ROC2", &var);
		TH1F* ptr = new TH1F(name, "ep Coincidence Time; ep Coincidence Time [ns]; Counts ", Bins::coin(bin), Bins::coin(bin+1), Bins::coin(bin+2));
		list->Add(ptr);
		return ptr;}
	//H_ep_ctime_total = new TH1F("H_ep_ctime_total", "ep Coincidence Time; ep Coincidence Time [ns]; Counts ", Bins::coin(bin), Bins::coin(bin+1), Bins::coin(bin+2));
	//H_ep_ctime_total_noCUT = new TH1F("H_ep_ctime_total_noCUT", "ep Coincidence Time; ep Coincidence Time [ns]; Counts ", Bins::coin(bin), Bins::coin(bin+1), Bins::coin(bin+2));

	//---------------------
	//HMS DETECTORS
	//--------------------
	static TH1F* mk_hCerNpeSum(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.cer.npeSum", &var);
		TH1F* ptr = new TH1F(name, "HMS Cherenkov NPE Sum; Cherenkov NPE Sum; Counts ", Bins::hcer(bin), Bins::hcer(bin+1), Bins::hcer(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hCalEtotNorm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.cal.etotnorm", &var);
		TH1F* ptr = new TH1F(name, "HMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}; Counts ", Bins::hcal(bin), Bins::hcal(bin+1), Bins::hcal(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hCalEtotTrkNorm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.cal.etottracknorm", &var);
		TH1F* ptr = new TH1F(name, "HMS Calorimeter Total Normalized Track Energy; E_{tot} / P_{trk}; Counts ", Bins::hcal(bin), Bins::hcal(bin+1), Bins::hcal(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hHodBetaNtrk(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.hod.betanotrack", &var);
		TH1F* ptr = new TH1F(name, "HMS Hodo #beta (no track); #beta (no track); Counts ", Bins::hbeta(bin), Bins::hbeta(bin+1), Bins::hbeta(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hHodBetaTrk(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.gtr.beta", &var);
		TH1F* ptr = new TH1F(name, "HMS Hodo #beta (golden track); #beta (golden track); Counts ", Bins::hbeta(bin), Bins::hbeta(bin+1), Bins::hbeta(bin+2));
		list->Add(ptr);
		return ptr;}

	//-------------------
	//SHMS DETECTORS
	//--------------------
	static TH1F* mk_pNGCerNpeSum(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.ngcer.npeSum", &var);
		TH1F* ptr = new TH1F(name, "SHMS Noble Gas Cherenkov NPE Sum; Cherenkov NPE Sum; Counts  ", Bins::pngcer(bin), Bins::pngcer(bin+1), Bins::pngcer(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_pHGCerNpeSum(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.hgcer.npeSum", &var);
		TH1F* ptr = new TH1F(name, "SHMS Heavy Gas Cherenkov NPE Sum; Cherenkov NPE Sum; Counts  ", Bins::phgcer(bin), Bins::phgcer(bin+1), Bins::phgcer(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_pCalEtotNorm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.cal.etotnorm", &var);
		TH1F* ptr = new TH1F(name, "SHMS Calorimeter Normalized Total Energy; E_{tot} / P_{cent}; Counts ", Bins::pcal(bin), Bins::pcal(bin+1), Bins::pcal(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_pCalEtotTrkNorm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.cal.etottracknorm", &var);
		TH1F* ptr = new TH1F(name, "SHMS Calorimeter Total Normalized Track Energy; E_{tot} / P_{trk}; Counts ", Bins::pcal(bin), Bins::pcal(bin+1), Bins::pcal(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_pHodBetaNtrk(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.hod.betanotrack", &var);
		TH1F* ptr = new TH1F(name, "SHMS Hodo #beta (no track); #beta (no track); Counts ", Bins::pbeta(bin), Bins::pbeta(bin+1), Bins::pbeta(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_pHodBetaTrk(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.gtr.beta", &var);
		TH1F* ptr = new TH1F(name, "SHMS Hodo #beta (golden track); #beta (golden track); Counts ", Bins::pbeta(bin), Bins::pbeta(bin+1), Bins::pbeta(bin+2));
		list->Add(ptr);
		return ptr;}


	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-------------------------------------------------------------Primary Kinematics (electron kinematics)-------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	
	//================================================================================
	//-------------------------------- Raw .root Files -------------------------------
	//================================================================================

	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	static TH1F* mk_the(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.scat_ang_rad", &var);
		TH1F* ptr = new TH1F(name, "Electron Scattering Angle; #theta_{e} [deg]", Bins::the(bin), Bins::the(bin+1), Bins::the(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_W(TTree* input, double &var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.W", &var);
		TH1F* ptr = new TH1F(name, "Invariant Mass; W [GeV]", Bins::W(bin), Bins::W(bin+1), Bins::W(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_W2(TTree* input, double &var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.W2", &var);
		TH1F* ptr = new TH1F(name, "Invariant Mass; W^{2} [(GeV)^{2}]", Bins::W2(bin), Bins::W2(bin+1), Bins::W2(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Q2(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.Q2", &var);
		TH1F* ptr = new TH1F(name, "4-Momentum Transfer; Q^{2} [(GeV/c)^{2}]", Bins::Q2(bin), Bins::Q2(bin+1), Bins::Q2(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_xbj(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.x_bj", &var);
		TH1F* ptr = new TH1F(name, "x-Bjorken; x_{B}", Bins::X(bin), Bins::X(bin+1), Bins::X(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_nu(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.nu", &var);
		TH1F* ptr = new TH1F(name, "Energy Transfer; #nu [GeV/c]", Bins::nu(bin), Bins::nu(bin+1), Bins::nu(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_q(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.q3m", &var);
		TH1F* ptr = new TH1F(name, "3-Momentum Transfer; |#vec{q}| [GeV/c]", Bins::q(bin), Bins::q(bin+1), Bins::q(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_qx(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.q_x", &var);
		TH1F* ptr = new TH1F(name, "3-Momentum Transfer (X-Component); |#vec{q}_{x}| [GeV/c]", Bins::qx(bin), Bins::qx(bin+1), Bins::qx(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_qy(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.q_y", &var);
		TH1F* ptr = new TH1F(name, "3-Momentum Transfer (Y-Component); |#vec{q}_{y}| [GeV/c]", Bins::qy(bin), Bins::qy(bin+1), Bins::qy(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_qz(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.q_z", &var);
		TH1F* ptr = new TH1F(name, "3-Momentum Transfer (Z-Component); |#vec{q}_{z}| [GeV/c]", Bins::qz(bin), Bins::qz(bin+1), Bins::qz(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thq(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.th_q", &var);
		TH1F* ptr = new TH1F(name, "In-Plane angle between #vec{q} and beam (+Z); #theta_{q} [deg]", Bins::thq(bin), Bins::thq(bin+1), Bins::thq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_phq(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("P.kin.primary.ph_q", &var);
		TH1F* ptr = new TH1F(name, "Out-of-Plane angle between #vec{q} and the beam (+Z); #phi_{q} [deg]", Bins::phq(bin), Bins::phq(bin+1), Bins::phq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_epsilon(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.kin.primary.epsilon", &var);
		TH1F* ptr = new TH1F(name, "Virtual Photon Polarization Factor; #epsilon", Bins::epsilon(bin), Bins::epsilon(bin + 1), Bins::epsilon(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_omega(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.kin.primary.omega", &var);
		TH1F* ptr = new TH1F(name, "Energy Transfer; #omega [GeV]", Bins::omega(bin), Bins::omega(bin + 1), Bins::omega(bin + 2));
		list->Add(ptr);
		return ptr;}
	
	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	static TH1F* mk_ki(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "Initial e^{-} Momentum; K_{i} [GeV/c]", Bins::ki(bin), Bins::ki(bin + 1), Bins::ki(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_kf(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "Final e^{-} Momentum; K_{f} [GeV/c]", Bins::kf(bin), Bins::kf(bin+1), Bins::kf(bin+2));
		list->Add(ptr);
		return ptr;}
	
	
	
	//================================================================================
	//-------------------------------- Simulation Files ------------------------------
	//================================================================================

	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	static TH1F* mk_the_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("theta_e", &var);
		TH1F* ptr = new TH1F(name, "Electron Scattering Angle; #theta_{e} [deg]", Bins::the(bin), Bins::the(bin + 1), Bins::the(bin + 2));
		list->Add(ptr);
		return ptr;	}
	static TH1F* mk_Q2_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("Q2", &var);
		TH1F* ptr = new TH1F(name, "4-Momentum Transfer; Q^{2} [(GeV/c)^{2}]", Bins::Q2(bin), Bins::Q2(bin + 1), Bins::Q2(bin + 2));
		list->Add(ptr);
		return ptr;	}
	static TH1F* mk_nu_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("nu", &var);
		TH1F* ptr = new TH1F(name, "Energy Transfer; #nu [GeV/c]", Bins::nu(bin), Bins::nu(bin + 1), Bins::nu(bin + 2));
		list->Add(ptr);
		return ptr;	}
	static TH1F* mk_q_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("q", &var);
		TH1F* ptr = new TH1F(name, "3-Momentum Transfer; |#vec{q}| [GeV/c]", Bins::q(bin), Bins::q(bin + 1), Bins::q(bin + 2));
		list->Add(ptr);
		return ptr;	}
	static TH1F* mk_W_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("W", &var);
		TH1F* ptr = new TH1F(name, "Invariant Mass; W [GeV]", Bins::W(bin), Bins::W(bin + 1), Bins::W(bin + 2));
		list->Add(ptr);
		return ptr;	}
	static TH1F* mk_kf_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_pf", &var);
		TH1F* ptr = new TH1F(name, "Final e^{-} Momentum; K_{f} [GeV/c]", Bins::kf(bin), Bins::kf(bin + 1), Bins::kf(bin + 2));
		list->Add(ptr);
		return ptr;	}

	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	static TH1F* mk_ki_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "Initial e^{-} Momentum; K_{i} [GeV/c]", Bins::kf(bin), Bins::kf(bin + 1), Bins::kf(bin + 2));
		list->Add(ptr);
		return ptr;	}
	static TH1F* mk_xbj_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "x-Bjorken; x_{B}", Bins::X(bin), Bins::X(bin + 1), Bins::X(bin + 2));
		list->Add(ptr);
		return ptr;	}
	static TH1F* mk_thq_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "In-Plane angle between #vec{q} and beam (+Z); #theta_{q} [deg]", Bins::thq(bin), Bins::thq(bin + 1), Bins::thq(bin + 2));
		list->Add(ptr);
		return ptr;	}



	
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-------------------------------------------------------------------Secondary Kinematics (Hadron)------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	
	//================================================================================
	//-------------------------------- Raw .root Files -------------------------------
	//================================================================================

	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	static TH1F* mk_Em(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.emiss", &var);
		TH1F* ptr = new TH1F(name, "Missing Energy; emiss [GeV]", Bins::Em(bin), Bins::Em(bin+1), Bins::Em(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Em_nuc(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.emiss_nuc", &var);
		TH1F* ptr = new TH1F(name, "Missing Energy (Nuclear Physics); emiss_nuc [GeV]", Bins::Em_nuc(bin), Bins::Em_nuc(bin+1), Bins::Em_nuc(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.pmiss", &var);
		TH1F* ptr = new TH1F(name, "Missing Momentum; P_{miss} [GeV]", Bins::Pm(bin), Bins::Pm(bin+1), Bins::Pm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pmx_lab(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.Prec_x", &var);
		TH1F* ptr = new TH1F(name, "P_{miss, x} (Lab); P_{miss, x} [GeV]", Bins::Pmx_lab(bin), Bins::Pmx_lab(bin+1), Bins::Pmx_lab(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pmy_lab(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.Prec_y", &var);
		TH1F* ptr = new TH1F(name, "P_{miss, y} (Lab); P_{miss, y} [GeV]", Bins::Pmy_lab(bin), Bins::Pmy_lab(bin+1), Bins::Pmy_lab(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pmz_lab(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.Prec_z", &var);
		TH1F* ptr = new TH1F(name, "P_{miss, z} (Lab); P_{miss, z} [GeV]", Bins::Pmz_lab(bin), Bins::Pmz_lab(bin+1), Bins::Pmz_lab(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pmx_q(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.pmiss_x", &var);
		TH1F* ptr = new TH1F(name, "P_{miss, xq} (wrt #vec{q}); P_{miss, xq} [GeV]", Bins::Pmx_q(bin), Bins::Pmx_q(bin+1), Bins::Pmx_q(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pmy_q(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.pmiss_y", &var);
		TH1F* ptr = new TH1F(name, "P_{miss, yq} (wrt #vec{q}); P_{miss, yq} [GeV]", Bins::Pmy_q(bin), Bins::Pmy_q(bin+1), Bins::Pmy_q(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pmz_q(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.pmiss_z", &var);
		TH1F* ptr = new TH1F(name, "P_{miss, zq} (wrt #vec{q}); P_{miss, zq} [GeV]", Bins::Pmz_q(bin), Bins::Pmz_q(bin+1), Bins::Pmz_q(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Tx(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.tx", &var);
		TH1F* ptr = new TH1F(name, "Kinetic Energy, T_{p} (detected); T_{p} [GeV]", Bins::Tx(bin), Bins::Tx(bin+1), Bins::Tx(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Tr(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.tb", &var);
		TH1F* ptr = new TH1F(name, "Kinetic Energy, T_{r} (recoil); T_{r} [GeV]", Bins::Tr(bin), Bins::Tr(bin+1), Bins::Tr(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_MM(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.MMp", &var);
		TH1F* ptr = new TH1F(name, "Missing Mass; MMp [GeV/c^{2}]", Bins::MM(bin), Bins::MM(bin+1), Bins::MM(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thpq(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.th_xq", &var);
		TH1F* ptr = new TH1F(name, "In-Plane (detected) Angle; #theta_{pq} [deg]", Bins::thxq(bin), Bins::thxq(bin+1), Bins::thxq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thrq(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.th_bq", &var);
		TH1F* ptr = new TH1F(name, "In-Plane (recoil) Angle; #theta_{rq} [deg]", Bins::thrq(bin), Bins::thrq(bin+1), Bins::thrq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_cthrq(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.kin.secondary.th_bq", &var);
		TH1F* ptr = new TH1F(name, "Cos of In-Plane (recoil) Angle; cos(#theta_{rq})", Bins::thrq(bin), -1.1, 1.1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_phpq(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.ph_xq", &var);
		TH1F* ptr = new TH1F(name, "Out-of-Plane (detected) Angle; #phi_{pq} [deg]", Bins::phxq(bin), Bins::phxq(bin+1), Bins::phxq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_phrq(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.ph_bq", &var);
		TH1F* ptr = new TH1F(name, "Out-of-Plane (recoil) Angle; #phi_{rq} [deg]", Bins::phrq(bin), Bins::phrq(bin+1), Bins::phrq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_xangle(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.xangle", &var);
		TH1F* ptr = new TH1F(name, "Angle Between The Detected Particle And Scattered Electron; xangle [Deg]", Bins::xangle(bin), Bins::xangle(bin+1), Bins::xangle(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Tx_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.tx_cm", &var);
		TH1F* ptr = new TH1F(name, "Kinetic Energy, T_{p, cm} (detected); T_{p, cm} [GeV]", Bins::Tx_cm(bin), Bins::Tx_cm(bin+1), Bins::Tx_cm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Tr_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.tb_cm", &var);
		TH1F* ptr = new TH1F(name, "Kinetic Energy T_{r, cm} (recoil); T_{r, cm} [GeV]", Bins::Tr_cm(bin), Bins::Tr_cm(bin+1), Bins::Tr_cm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thxq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.thx_cm", &var);
		TH1F* ptr = new TH1F(name, "In-Plane (detected) Angle; #theta_{pq, cm} [deg]", Bins::thxq_cm(bin), Bins::thxq_cm(bin+1), Bins::thxq_cm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thrq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.thb_cm", &var);
		TH1F* ptr = new TH1F(name, "In-Plane (recoil) Angle; #theta_{rq, cm} [deg]", Bins::thrq_cm(bin), Bins::thrq_cm(bin+1), Bins::thrq_cm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_phxq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.phx_cm", &var);
		TH1F* ptr = new TH1F(name, "Out-of-Plane (detected) Angle; #phi_{pq, cm} [deg]", Bins::phxq_cm(bin), Bins::phxq_cm(bin+1), Bins::phxq_cm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_phrq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.phb_cm", &var);
		TH1F* ptr = new TH1F(name, "Out-of-Plane (recoil) Angle; #phi_{rq, cm} [cm]", Bins::phrq_cm(bin), Bins::phrq_cm(bin+1), Bins::phrq_cm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Ttot_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.t_tot_cm", &var);
		TH1F* ptr = new TH1F(name, "Total CM Kinetic Energy; T_{tot, cm} [GeV]", Bins::Ttot_cm(bin), Bins::Ttot_cm(bin+1), Bins::Ttot_cm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_MandelS(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.MandelS", &var);
		TH1F* ptr = new TH1F(name, "s-Mandelstam; MandelS [(Gev)^{2}]", Bins::MandelS(bin), Bins::MandelS(bin+1), Bins::MandelS(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_MandelT(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.MandelT", &var);
		TH1F* ptr = new TH1F(name, "t-Mandelstam; MandelS [(Gev)^{2}]", Bins::MandelT(bin), Bins::MandelT(bin+1), Bins::MandelT(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_MandelU(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("H.kin.secondary.MandelU", &var);
		TH1F* ptr = new TH1F(name, "u-Mandelstam; MandelS [(Gev)^{2}]", Bins::MandelU(bin), Bins::MandelU(bin+1), Bins::MandelU(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Erecoil(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.kin.secondary.Erecoil", &var);
		TH1F* ptr = new TH1F(name, "Total Energy of Recoil System; Erecoil [GeV/c]", Bins::Erecoil(bin), Bins::Erecoil(bin + 1), Bins::Erecoil(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_MMk(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.kin.secondary.MMk", &var);
		TH1F* ptr = new TH1F(name, "Missing Mass (kaon); MMk [GeV/c^{2}]", Bins::MMk(bin), Bins::MMk(bin + 1), Bins::MMk(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_MMpi(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.kin.secondary.MMpi", &var);
		TH1F* ptr = new TH1F(name, "Missing Mass (pion); MMpi [GeV/c^{2}]", Bins::MMpi(bin), Bins::MMpi(bin + 1), Bins::MMpi(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Mrecoil(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.kin.secondary.Mrecoil", &var);
		TH1F* ptr = new TH1F(name, "Invariant Mass of the Recoil System; Mrecoil [GeV/c^{2}]", Bins::Mrecoil(bin), Bins::Mrecoil(bin + 1), Bins::Mrecoil(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_px_cm(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.kin.secondary.px_cm", &var);
		TH1F* ptr = new TH1F(name, "Magnitude of X momentum in CM System; px_cm [GeV]", Bins::px_cm(bin), Bins::px_cm(bin + 1), Bins::px_cm(bin + 2));
		list->Add(ptr);
		return ptr;}

	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	static TH1F* mk_Em_src(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "Missing Energy (SRC Nuclear Physics); Em_src [GeV]", Bins::Em_nuc(bin), Bins::Em_nuc(bin+1), Bins::Em_nuc(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_MM2(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "Missing Mass Squared; MM2 [(GeV/c^{2})^{2}]", Bins::MM2(bin), Bins::MM2(bin+1), Bins::MM2(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thp(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "Hadron Scattering Angle (detected); #theta_{p} [deg]", Bins::thx(bin), Bins::thx(bin+1), Bins::thx(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pf(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "Final Hadron Momentum (detected); P_{f} [GeV/c]", Bins::Pf(bin), Bins::Pf(bin+1), Bins::Pf(bin+2));
		list->Add(ptr);
		return ptr;}

	


	//================================================================================
	//-------------------------------- Simulation Files ------------------------------
	//================================================================================

	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------
	static TH1F* mk_Pm_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("Pm", &var);
		TH1F* ptr = new TH1F(name, "Missing Momentum; P_{miss} [GeV/c]", Bins::Pm(bin), Bins::Pm(bin+1), Bins::Pm(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thpq_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("theta_pq", &var);
		TH1F* ptr = new TH1F(name, "In-Plane (detected) Angle; #theta_{pq} [deg]", Bins::thxq(bin), Bins::thxq(bin+1), Bins::thxq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thrq_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("theta_rq", &var);
		TH1F* ptr = new TH1F(name, "In-Plane (recoil) Angle; #theta_{rq} [deg]", Bins::thrq(bin), Bins::thrq(bin+1), Bins::thrq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_phpq_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("phi_pq", &var);
		TH1F* ptr = new TH1F(name, "Out-of-Plane (detected) Angle; #phi_{pq} [deg]", Bins::phxq(bin), Bins::phxq(bin+1), Bins::phxq(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Pf_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("h_pf", &var);
		TH1F* ptr = new TH1F(name, "Final Hadron Momentum (detected); P_{f} [GeV/c]", Bins::px_cm(bin), Bins::px_cm(bin + 1), Bins::px_cm(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_thp_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("theta_p", &var);
		TH1F* ptr = new TH1F(name, "Hadron Scattering Angle (detected); #theta_{p} [deg]", Bins::px_cm(bin), Bins::px_cm(bin + 1), Bins::px_cm(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Em_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		input->SetBranchAddress("Em", &var);
		TH1F* ptr = new TH1F(name, "Missing Energy; emiss [GeV]", Bins::Em(bin), Bins::Em(bin + 1), Bins::Em(bin + 2));
		list->Add(ptr);
		return ptr;}
	
	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	static TH1F* mk_MM_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "Missing Mass; MM [GeV/c^{2}]", Bins::MM(bin), Bins::MM(bin+1), Bins::MM(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_MM2_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "Missing Mass Squared; MM2 [(GeV/c^{2})^{2}]", Bins::MM2(bin), Bins::MM2(bin+1), Bins::MM2(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_Ep_sim(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "Final Proton Energy; Ep [GeV]", Bins::Ep(bin), Bins::Ep(bin + 1), Bins::Ep(bin + 2));
		list->Add(ptr);
		return ptr;}



	//================================================================================
	//(Cosine, Sine) Histos of detected AND recoil angles (range is fixed at: -1, 1) (Calculated)
	//================================================================================

	//----- LAB FRAME -----
	static TH1F* mk_cth_xq(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "cos(#theta_{pq})", Bins::thxq(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_cth_rq(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "cos(#theta_{rq})", Bins::thrq(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_sth_xq(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "sin(#theta_{pq})", Bins::thxq(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_sth_rq(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "sin(#theta_{rq})", Bins::thrq(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_cphi_xq(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "cos(#phi_{pq})", Bins::phxq(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_cphi_rq(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "cos(#phi_{rq})", Bins::phrq(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_sphi_xq(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "sin(#phi_{pq})", Bins::phxq(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_sphi_rq(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "sin(#phi_{rq})", Bins::phrq(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	
	//----- CM FRAME -----
	static TH1F* mk_cth_xq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "cos(#theta_{pq,cm})", Bins::thxq_cm(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_cth_rq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "cos(#theta_{rq,cm})", Bins::thrq_cm(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_sth_xq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "sin(#theta_{pq,cm})", Bins::thxq_cm(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_sth_rq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "sin(#theta_{rq,cm})", Bins::thrq_cm(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_cphi_xq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "cos(#phi_{pq,cm})", Bins::phxq_cm(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_cphi_rq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "cos(#phi_{rq,cm})", Bins::phrq_cm(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_sphi_xq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "sin(#phi_{pq,cm})", Bins::phxq_cm(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_sphi_rq_cm(TTree* input, double& var, const char* name, TList* list, int bin){
		TH1F* ptr = new TH1F(name, "sin(#phi_{rq,cm})", Bins::phrq_cm(bin), -1, 1);
		list->Add(ptr);
		return ptr;}
	



	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------------Acceptance Histogram Bins---------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	
	//-------------------
	//HMS
	//-------------------
	//----- Hadron Arm Focal Plane -----
	static TH1F* mk_hxfp(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.dc.x_fp", &var);
		TH1F* ptr = new TH1F(name, "HMS X_{fp}; X_{fp} [cm]; Counts", Bins::hxfp(bin), Bins::hxfp(bin+1), Bins::hxfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hxpfp(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.dc.xp_fp", &var);
		TH1F* ptr = new TH1F(name, "HMS X'_{fp}; X'_{fp} [deg]; Counts", Bins::hxpfp(bin), Bins::hxpfp(bin+1), Bins::hxpfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hyfp(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.dc.y_fp", &var);
		TH1F* ptr = new TH1F(name, "HMS Y_{fp}; Y_{fp} [cm]; Counts", Bins::hyfp(bin), Bins::hyfp(bin+1), Bins::hyfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hypfp(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.dc.yp_fp", &var);
		TH1F* ptr = new TH1F(name, "HMS Y'_{fp}; Y'_{fp} [deg]; Counts", Bins::hypfp(bin), Bins::hypfp(bin+1), Bins::hypfp(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- Hadron Arm Reconstructed Quantities -----
	static TH1F* mk_hytar(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.gtr.y", &var);
		TH1F* ptr = new TH1F(name, "HMS Y_{tar}; Y_{tar} [cm]; Counts", Bins::hytar(bin), Bins::hytar(bin+1), Bins::hytar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hyptar(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.gtr.ph", &var);
		TH1F* ptr = new TH1F(name, "HMS Y'_{tar}; Y'_{tar} [deg]; Counts", Bins::hyptar(bin), Bins::hyptar(bin+1), Bins::hyptar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hxptar(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.gtr.th", &var);
		TH1F* ptr = new TH1F(name, "HMS X'_{tar}; X'_{tar} [deg]; Counts", Bins::hxptar(bin), Bins::hxptar(bin+1), Bins::hxptar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hdelta(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.gtr.dp", &var);
		TH1F* ptr = new TH1F(name, "HMS Momentum Acceptance, #delta; #delta [%]; Counts", Bins::hdelta(bin), Bins::hdelta(bin+1), Bins::hdelta(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- Target Reconstruction (Hall Coord. System) -----
	static TH1F* mk_htarx(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.react.x", &var);
		TH1F* ptr = new TH1F(name, "HMS x-Target (Lab); x-Target [cm]; Counts", Bins::tarx(bin), Bins::tarx(bin+1), Bins::tarx(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_htary(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.react.y", &var);
		TH1F* ptr = new TH1F(name, "HMS y_Target (Lab); y-Target [cm]; Counts", Bins::tary(bin), Bins::tary(bin+1), Bins::tary(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_htarz(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.react.z", &var);
		TH1F* ptr = new TH1F(name, "HMS z-Target (Lab); z-Target [cm]; Counts", Bins::tarz(bin), Bins::tarz(bin+1), Bins::tarz(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- HMS Collimator -----
	static TH1F* mk_hXColl(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.extcor.xsieve", &var);
		TH1F* ptr = new TH1F(name, "HMS X Collimator; X-Collimator [cm]; Counts ", Bins::hXColl(bin), Bins::hXColl(bin+1), Bins::hXColl(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hYColl(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("H.extcor.ysieve", &var);
		TH1F* ptr = new TH1F(name, "HMS Y Collimator; Y-Collimator [cm]; Counts ", Bins::hYColl(bin), Bins::hYColl(bin+1), Bins::hYColl(bin+2));
		list->Add(ptr);
		return ptr;}

	//-------------------
	//SHMS
	//-------------------
	//----- Electron Arm Focal Plane
	static TH1F* mk_exfp(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.dc.x_fp", &var);
		TH1F* ptr = new TH1F(name, "SHMS X_{fp}; X_{fp} [cm]; Counts", Bins::exfp(bin), Bins::exfp(bin+1), Bins::exfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_expfp(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.dc.xp_fp", &var);
		TH1F* ptr = new TH1F(name, "SHMS X'_{fp}; X'_{fp} [deg]; Counts", Bins::expfp(bin), Bins::expfp(bin+1), Bins::expfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_eyfp(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.dc.y_fp", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y_{fp}; Y_{fp} [cm]; Counts", Bins::eyfp(bin), Bins::eyfp(bin+1), Bins::eyfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_eypfp(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.dc.yp_fp", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y'_{fp}; Y'_{fp} [deg]; Counts", Bins::eypfp(bin), Bins::eypfp(bin+1), Bins::eypfp(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- Electron Arm Reconstructed Quantities -----
	static TH1F* mk_eytar(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.gtr.y", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y_{tar}; Y_{tar} [cm]; Counts", Bins::eytar(bin), Bins::eytar(bin+1), Bins::eytar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_eyptar(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.gtr.ph", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y'_{tar}; Y'_{tar} [deg]; Counts", Bins::eyptar(bin), Bins::eyptar(bin+1), Bins::eyptar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_exptar(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.gtr.th", &var);
		TH1F* ptr = new TH1F(name, "SHMS X'_{tar}; X'_{tar} [deg]; Counts", Bins::exptar(bin), Bins::exptar(bin+1), Bins::exptar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_edelta(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.gtr.dp", &var);
		TH1F* ptr = new TH1F(name, "SHMS Momentum Acceptance, #delta; #delta [%]; Counts", Bins::edelta(bin), Bins::edelta(bin+1), Bins::edelta(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- Target Reconstruction (Hall Coord. System)
	static TH1F* mk_etarx(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.react.x", &var);
		TH1F* ptr = new TH1F(name, "SHMS x-Target (Lab); x-Target [cm]; Counts", Bins::tarx(bin), Bins::tarx(bin+1), Bins::tarx(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_etary(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.react.y", &var);
		TH1F* ptr = new TH1F(name, "SHMS y-Target (Lab); y-Target [cm]; Counts", Bins::tary(bin), Bins::tary(bin+1), Bins::tary(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_etarz(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.react.z", &var);
		TH1F* ptr = new TH1F(name, "SHMS z-Target (Lab); z-Target [cm]; Counts", Bins::tarz(bin), Bins::tarz(bin+1), Bins::tarz(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- SHMS Collimator -----
	static TH1F* mk_eXColl(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.extcor.xsieve", &var);
		TH1F* ptr = new TH1F(name, "SHMS X Collimator; X-Collimator [cm]; Counts", Bins::eXColl(bin), Bins::eXColl(bin+1), Bins::eXColl(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_eYColl(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("P.extcor.ysieve", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y Collimator; Y-Collimator [cm]; Counts", Bins::eYColl(bin), Bins::eYColl(bin+1), Bins::eYColl(bin+2));
		list->Add(ptr);
		return ptr;}

	//-------------------
	//Calculated Quantities
	//-------------------
	//----- difference in reaction vertex z
	static TH1F* mk_ztar_diff(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "z-Vertex Difference; (HMS-SHMS) z-Vertex Difference [cm]; Counts", Bins::ztar_diff(bin), Bins::ztar_diff(bin+1), Bins::ztar_diff(bin+2));
		list->Add(ptr);
		return ptr;}
	






	//================================================================================
	//-------------------------------- Simulation Files ------------------------------
	//================================================================================

	//----------------------------------------
	//-------------- Given Leafs -------------
	//----------------------------------------

	//--------------------
	//-------- HMS -------
	//--------------------
	//----- Hadron Arm Focal Plane -----
	static TH1F* mk_hxfp_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_xfp", &var);
		TH1F* ptr = new TH1F(name, "HMS X_{fp}; X_{fp} [cm]; Counts / mC", Bins::hxfp(bin), Bins::hxfp(bin+1), Bins::hxfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hxpfp_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_xpfp", &var);
		TH1F* ptr = new TH1F(name, "HMS X'_{fp}; X'_{fp} [deg]; Counts / mC", Bins::hxfp(bin), Bins::hxfp(bin+1), Bins::hxfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hyfp_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_yfp", &var);
		TH1F* ptr = new TH1F(name, "HMS X'_{fp}; X'_{fp} [cm]; Counts / mC", Bins::hyfp(bin), Bins::hyfp(bin+1), Bins::hyfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hypfp_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_ypfp", &var);
		TH1F* ptr = new TH1F(name, "HMS Y'_{fp}; Y'_{fp} [deg]; Counts / mC", Bins::hypfp(bin), Bins::hypfp(bin+1), Bins::hypfp(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- Hadron Arm Reconstructed Quantities -----
	static TH1F* mk_hytar_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_ytar", &var);
		TH1F* ptr = new TH1F(name, "HMS Y_{tar}; Y_{tar} [cm]; Counts / mC", Bins::hytar(bin), Bins::hytar(bin+1), Bins::hytar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hyptar_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_yptar", &var);
		TH1F* ptr = new TH1F(name, "HMS Y'_{tar}; Y'_{tar} [deg]; Counts / mC", Bins::hyptar(bin), Bins::hyptar(bin+1), Bins::hyptar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hxptar_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_xptar", &var);
		TH1F* ptr = new TH1F(name, "HMS X'_{tar}; X'_{tar} [deg]; Counts / mC", Bins::hxptar(bin), Bins::hxptar(bin+1), Bins::hxptar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hdelta_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_delta", &var);
		TH1F* ptr = new TH1F(name, "HMS Momentum Acceptance, #delta; #delta [%]; Counts / mC", Bins::hdelta(bin), Bins::hdelta(bin+1), Bins::hdelta(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- Target Reconstruction (Hall Coord. System) -----
	static TH1F* mk_htarx_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("tar_x", &var);
		TH1F* ptr = new TH1F(name, "HMS x-Target (Lab); x-Target [cm]; Counts / mC", Bins::tarx(bin), Bins::tarx(bin+1), Bins::tarx(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_htary_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_yv", &var);
		TH1F* ptr = new TH1F(name, "HMS y-Target (Lab); y-Target [cm]; Counts / mC", Bins::tary(bin), Bins::tary(bin+1), Bins::tary(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_htarz_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("h_zv", &var);
		TH1F* ptr = new TH1F(name, "HMS z-Target (Lab); z-Target [cm]; Counts / mC", Bins::tarz(bin), Bins::tarz(bin+1), Bins::tarz(bin+2));
		list->Add(ptr);
		return ptr;}
	

	//--------------------
	//------- SHMS -------
	//--------------------
	//----- Electron Arm Focal Plane
	static TH1F* mk_exfp_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_xfp", &var);
		TH1F* ptr = new TH1F(name, "SHMS X_{fp}; X_{fp} [cm]; Counts / mC", Bins::exfp(bin), Bins::exfp(bin+1), Bins::exfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_expfp_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_xpfp", &var);
		TH1F* ptr = new TH1F(name, "SHMS X'_{fp}; X'_{fp} [deg]; Counts / mC", Bins::expfp(bin), Bins::expfp(bin+1), Bins::expfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_eyfp_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_yfp", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y_{fp}; Y_{fp} [cm]; Counts / mC", Bins::eyfp(bin), Bins::eyfp(bin+1), Bins::eyfp(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_eypfp_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_ypfp", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y'_{fp}; Y'_{fp} [deg]; Counts / mC", Bins::eypfp(bin), Bins::eypfp(bin+1), Bins::eypfp(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- Electron Arm Reconstructed Quantities -----
	static TH1F* mk_eytar_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_ytar", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y_{tar}; Y_{tar} [cm]; Counts / mC", Bins::eytar(bin), Bins::eytar(bin+1), Bins::eytar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_eyptar_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_yptar", &var);
		TH1F* ptr = new TH1F(name, "SHMS Y'_{tar}; Y'_{tar} [deg]; Counts / mC", Bins::eyptar(bin), Bins::eyptar(bin+1), Bins::eyptar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_exptar_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_xptar", &var);
		TH1F* ptr = new TH1F(name, "SHMS X'_{tar}; X'_{tar} [deg]; Counts / mC", Bins::exptar(bin), Bins::exptar(bin+1), Bins::exptar(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_edelta_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_delta", &var);
		TH1F* ptr = new TH1F(name, "SHMS Momentum Acceptance, #delta; #delta [%]; Counts / mC", Bins::edelta(bin), Bins::edelta(bin+1), Bins::edelta(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- Target Reconstruction (Hall Coord. System)
	static TH1F* mk_etarx_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("tar_x", &var);
		TH1F* ptr = new TH1F(name, "SHMS x-Target (Lab); x-Target [cm]; Counts / mC", Bins::tarx(bin), Bins::tarx(bin+1), Bins::tarx(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_etary_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_yv", &var);
		TH1F* ptr = new TH1F(name, "SHMS y-Target (Lab); y-Target [cm]; Counts / mC", Bins::tary(bin), Bins::tary(bin+1), Bins::tary(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_etarz_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		input->SetBranchAddress("e_zv", &var);
		TH1F* ptr = new TH1F(name, "SHMS z-Target (Lab); z-Target [cm]; Counts / mC", Bins::tarz(bin), Bins::tarz(bin+1), Bins::tarz(bin+2));
		list->Add(ptr);
		return ptr;}
	
	//----------------------------------------
	//--------- Calculated Quantities --------
	//----------------------------------------
	//----- HMS Collimator -----
	static TH1F* mk_hXColl_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "HMS X Collimator; X-Collimator [cm]; Counts / mC", Bins::hXColl(bin), Bins::hXColl(bin+1), Bins::hXColl(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_hYColl_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "HMS Y Collimator; Y-Collimator [cm]; Counts / mC", Bins::hYColl(bin), Bins::hYColl(bin+1), Bins::hYColl(bin+2));
		list->Add(ptr);
		return ptr;}
	//----- SHMS Collimator -----
	static TH1F* mk_eXColl_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "SHMS X Collimator; X-Collimator [cm]; Counts / mC", Bins::eXColl(bin), Bins::eXColl(bin+1), Bins::eXColl(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_eYColl_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "SHMS Y Collimator; Y-Collimator [cm]; Counts / mC", Bins::eYColl(bin), Bins::eYColl(bin+1), Bins::eYColl(bin+2));
		list->Add(ptr);
		return ptr;}
	static TH1F* mk_ztar_diff_sim(TTree* input, double& var, const char* name, TList* list, int bin) {
		TH1F* ptr = new TH1F(name, "z-Vertex Difference; (HMS-SHMS) z-Vertex Difference [cm]; Counts / mC", Bins::ztar_diff(bin), Bins::ztar_diff(bin+1), Bins::ztar_diff(bin+2));
		list->Add(ptr);
		return ptr;}








	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//-------------------------------------------------------------------- Random COIN Background ----------------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//For Random COIN Background for ep_ctime, W, Q2, X, nu, q, Em, Em_nuc, Pm, MM, thxq, & thrq just use the above fcns


	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------- Random Subtracted COIN Background ----------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//For Random Subtracted COIN Background for ep_ctime, W, Q2, X, nu, q, Em, Em_nuc, Pm, MM, thxq, & thrq just use the above fcns






	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------- 2D Acceptance Histograms ----------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static TH2F* mk_cthrq_pmiss(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "cos(#theta_{rq}) vs. P_{m}; P_{m} [GeV/c]; cos(#theta_{rq})", Bins::Pm(bin), Bins::Pm(bin + 1), Bins::Pm(bin + 2), Bins::thrq(bin), -1.1, 1.1);
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_hXColl_hYColl(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "HMS Collimator; HMS Y-Collimator [cm]; HMS X-Collimator [cm]", Bins::hYColl(bin), Bins::hYColl(bin + 1), Bins::hYColl(bin + 2), Bins::hXColl(bin), Bins::hXColl(bin + 1), Bins::hXColl(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_eXColl_eYColl(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "SHMS Collimator; SHMS Y-Collimator [cm]; SHMS X-Collimator [cm]", Bins::eYColl(bin), Bins::eYColl(bin + 1), Bins::eYColl(bin + 2), Bins::eXColl(bin), Bins::eXColl(bin + 1), Bins::eXColl(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_hxfp_hyfp(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "HMS X_{fp} vs. Y_{fp}; Y_{fp} [cm]; X_{fp} [cm]", Bins::hyfp(bin), Bins::hyfp(bin + 1), Bins::hyfp(bin + 2), Bins::hxfp(bin), Bins::hxfp(bin + 1), Bins::hxfp(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_exfp_eyfp(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "SHMS X_{fp} vs. Y_{fp}; Y_{fp} [cm]; X_{fp} [cm]", Bins::eyfp(bin), Bins::eyfp(bin + 1), Bins::eyfp(bin + 2), Bins::exfp(bin), Bins::exfp(bin + 1), Bins::exfp(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_hxptar_exptar(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "SHMS vs. HMS; X'_{htar} [deg]; X'_{etar} [deg]", Bins::exptar(bin), Bins::exptar(bin + 1), Bins::exptar(bin + 2), Bins::hxptar(bin), Bins::hxptar(bin + 1), Bins::hxptar(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_hyptar_eyptar(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "SHMS vs. HMS; Y'_{htar} [deg]; Y'_{etar} [deg]", Bins::eyptar(bin), Bins::eyptar(bin + 1), Bins::eyptar(bin + 2), Bins::hyptar(bin), Bins::hyptar(bin + 1), Bins::hyptar(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_hdelta_edelta(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "SHMS vs. HMS; #delta_{h} [%]; #delta_{e} [%]", Bins::edelta(bin), Bins::edelta(bin + 1), Bins::edelta(bin + 2), Bins::hdelta(bin), Bins::hdelta(bin + 1), Bins::hdelta(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_W_pmiss(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "P_{m} vs. W; P_{m} [GeV/c]; W [GeV]", Bins::Pm(bin), Bins::Pm(bin + 1), Bins::Pm(bin + 2), Bins::W(bin), Bins::W(bin + 1), Bins::W(bin + 2));
		list->Add(ptr);
		return ptr;
	}


	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	//--------------------------------------------------------------- 2D Kinematic Histograms ----------------------------------------------------------------
	//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
	static TH2F* mk_thrq_pmiss(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "P_{m} vs. #theta_{rq} (yield); P_{m} [GeV/c]; #theta_{rq} [deg]", Bins::Pm(bin), Bins::Pm(bin + 1), Bins::Pm(bin + 2), Bins::thrq(bin), Bins::thrq(bin + 1), Bins::thrq(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_xbj_Q2(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "x_{B} vs. Q^{2}; Q^{2} [(GeV/c)^{2}]; x_{B}", Bins::Q2(bin), Bins::Q2(bin + 1), Bins::Q2(bin + 2), Bins::X(bin), Bins::X(bin + 1), Bins::X(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_Em_Pm(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "E_{m} vs. P_{m}; P_{m} [GeV/c]; E_{m} [GeV]", Bins::Pm(bin), Bins::Pm(bin + 1), Bins::Pm(bin + 2), Bins::Em_nuc(bin), Bins::Em_nuc(bin + 1), Bins::Em_nuc(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_the_kf(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "#theta_{e} vs. K_{f}; K_{f} [GeV/c]; #theta_{e} [deg]", Bins::kf(bin), Bins::kf(bin + 1), Bins::kf(bin + 2), Bins::the(bin), Bins::the(bin + 1), Bins::the(bin + 2));
		list->Add(ptr);
		return ptr;}
	static TH2F* mk_thp_Pf(TTree* input, const char* name, TList* list, int bin) {
		TH2F* ptr = new TH2F(name, "#theta_{p} vs. P_{f}; P_{f} [GeV/c]; #theta_{p} [deg]", Bins::Pf(bin), Bins::Pf(bin + 1), Bins::Pf(bin + 2), Bins::thq(bin), Bins::thq(bin + 1), Bins::thq(bin + 2));
		list->Add(ptr);
		return ptr;}

	private:
};

#endif