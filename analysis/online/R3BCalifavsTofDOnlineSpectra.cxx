/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BCalifavsTofDOnlineSpectra.h"
#include "R3BCalifaClusterData.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTofdHitData.h"
#include "R3BCalifaMappedData.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BFrsData.h"
#include "R3BLosHitData.h"
#include "R3BWRData.h"
#include "R3BFootHitData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

R3BCalifavsTofDOnlineSpectra::R3BCalifavsTofDOnlineSpectra()
    : R3BCalifavsTofDOnlineSpectra("CALIFAvsTofDOnlineSpectra", 1)
{
}

R3BCalifavsTofDOnlineSpectra::R3BCalifavsTofDOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHitItemsCalifa(NULL)
    , fHitItemsTofd(NULL)
    , fNEvents(0)
    , fTpat(-1)
    , fMinProtonE(50000.) // 50MeV
    , fZselection(5.)
{
}

R3BCalifavsTofDOnlineSpectra::~R3BCalifavsTofDOnlineSpectra()
{
    R3BLOG(debug1, "");
    if (fHitItemsCalifa)
        delete fHitItemsCalifa;
    if (fHitItemsTofd)
        delete fHitItemsTofd;
    if (fWRItemsCalifa)
	delete fWRItemsCalifa;
    if (fWRItemsMaster)
	delete fWRItemsMaster;
    if (fMapItemsCalifa)
	delete fMapItemsCalifa;
    if (fCalItemsCalifa)
	delete fCalItemsCalifa;
    if (fHitItemsLos)
	delete fHitItemsLos;
    if (fHitItemsFrs)
	delete fHitItemsFrs;
    if (fHitItemsFoot)
	delete fHitItemsFoot;
}

void R3BCalifavsTofDOnlineSpectra::SetParContainers()
{
    // Parameter Container
    // Reading amsStripCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");
}

bool R3BCalifavsTofDOnlineSpectra::isFootDetect(TVector3 hit)
{
    double theta = hit.Theta() * TMath::RadToDeg();
    double phi = hit.Phi() * TMath::RadToDeg();
    
    bool cond_theta = ((theta>20.0) && (theta<80.0));
    
    double h = 8.0265; // cm, perp dist from target
    double beta = 1.10915; // 63.55 degrees
    double foot_half_height = 5.0; // cm, foot half height
    double d = h / TMath::Cos(hit.Theta()-beta);
    double phiR_lim = 180.0 - TMath::ATan(foot_half_height / d) * TMath::RadToDeg();
    double phiL_lim = TMath::ATan(foot_half_height / d) * TMath::RadToDeg();
    bool cond_phi = ((TMath::Abs(phi)<phiL_lim) || (TMath::Abs(phi)>phiR_lim));

    return (cond_theta && cond_phi);
}

InitStatus R3BCalifavsTofDOnlineSpectra::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();

    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = (R3BEventHeader*)mgr->GetObject("EventHeader.");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // get access to Hit data
    fHitItemsCalifa = (TClonesArray*)mgr->GetObject("CalifaClusterData");
    R3BLOG_IF(fatal, !fHitItemsCalifa, "CalifaClusterData not found");

    fHitItemsTofd = (TClonesArray*)mgr->GetObject("TofdHit");
    R3BLOG_IF(warn, !fHitItemsTofd, "TofdHit not found");

    //gROOT->ProcessLine("#include <vector>");

    // create Root File
    outfile = new TFile("protons.root","RECREATE");

    /*outfile->cd();
    outtree = new TTree("genTbuffer","General Tree");

    // general vars
    outtree->Branch("tpat",&tpatval,"tpatval/I");

    // califa map vars
    outtree->Branch("map_crystalId",&map_crystalId);
    outtree->Branch("map_energy",&map_energy);
    outtree->Branch("map_Ns",&map_Ns);
    outtree->Branch("map_Nf",&map_Nf);
    outtree->Branch("map_febexTime",&map_febexTime);
    outtree->Branch("map_wrts",&map_wrts);
    outtree->Branch("map_overflow",&map_overflow);
    outtree->Branch("map_pileup",&map_pileup);
    outtree->Branch("map_discard",&map_discard);

    // califa cal vars
    outtree->Branch("cal_crystalId",&cal_crystalId);
    outtree->Branch("cal_energy",&cal_energy);
    outtree->Branch("cal_Ns",&cal_Ns);
    outtree->Branch("cal_Nf",&cal_Nf);
    outtree->Branch("cal_time",&cal_time);
    
    outtree->Branch("theta",&thetaList);
    outtree->Branch("phi",&phiList);
    outtree->Branch("energy",&energyList);
    outtree->Branch("Nf",&NfList);
    outtree->Branch("Ns",&NsList);
    outtree->Branch("crystalHits",&crystalHitsList);
    outtree->Branch("clusterId",&clusterIdList);
    outtree->Branch("time",&timeList);

    // tofd vars
    outtree->Branch("detectorId",&detectorId);
    outtree->Branch("charge",&charge);
    outtree->Branch("tof",&tof);

    outtree->Branch("loshits",&loshits,"loshits/I");
    outtree->Branch("frshits",&frshits,"frshits/I");

    outtree->Branch("wr",&wr);
    outtree->Branch("wrm",&wrm,"wrm/L");

    // foot vars
    outtree->Branch("foot_detectorId",&foot_detectorId);
    outtree->Branch("foot_nbHits",&foot_nbHits);
    outtree->Branch("foot_mulStrip",&foot_mulStrip);
    outtree->Branch("foot_pos",&foot_pos);
    outtree->Branch("foot_poslab",&foot_poslab);
    outtree->Branch("foot_theta",&foot_theta);
    outtree->Branch("foot_phi",&foot_phi);
    outtree->Branch("foot_energy",&foot_energy);

    //std::cout << "CalifavsTofDOnlineSpectra Init middle" << std::endl;*/

    // Create histograms for detectors

    // CANVAS Theta vs Phi
    /*cCalifa_angles = new TCanvas("Califa_Theta_vs_Phi_andTofD", "Theta vs Phi", 10, 10, 500, 500);
    cCalifa_angles->Divide(2, 2);
    //cCalifa_angles = new TCanvas("Califa_Theta_vs_Phi_andTofD", "Theta vs Phi", 10, 10, 500, 500);
    //cCalifa_angles->Divide(2, 1);

    int bins = 2000;
    double minE = 0.;
    double maxE = 1000000.; // high: 1 GeV low: 20 MeV

    for (int i = 0; i < 3; i++)
    {
        //cCalifa_angles->cd(i + 1);
        
	TString buf1;
	if (i == 0) buf1 = "fh2_theta_phi";
	if (i == 1) buf1 = "fh2_theta_phi_withTofD";
	if (i == 2) buf1 = "fh2_theta_phi_withTofDFoot";
        fh2_Califa_theta_phi[i] = new TH2F(buf1, buf1, 50, 0, 90, 180, -180, 180);
        fh2_Califa_theta_phi[i]->GetXaxis()->SetTitle("Theta [degrees]");
        fh2_Califa_theta_phi[i]->GetYaxis()->SetTitle("Phi [degrees]");
        fh2_Califa_theta_phi[i]->GetYaxis()->SetTitleOffset(1.2);
        fh2_Califa_theta_phi[i]->GetXaxis()->CenterTitle(true);
        fh2_Califa_theta_phi[i]->GetYaxis()->CenterTitle(true);
        fh2_Califa_theta_phi[i]->Draw("COLZ");

	TString buf2;
	if (i == 0) buf2 = "fh2_NsNf";
	if (i == 1) buf2 = "fh2_NsNf_withTofD";
	if (i == 2) buf2 = "fh2_NsNf_withTofDFoot";
	fh2_Califa_NsNf[i] =
		new TH2F(buf2, buf2, bins, minE, maxE, bins, minE, maxE);
	fh2_Califa_NsNf[i]->GetXaxis()->SetTitle("Nf Energy [keV]");
	fh2_Califa_NsNf[i]->GetYaxis()->SetTitle("Ns Energy [keV]");
	fh2_Califa_NsNf[i]->GetYaxis()->SetTitleOffset(1.4);
	fh2_Califa_NsNf[i]->GetXaxis()->CenterTitle(true);
	fh2_Califa_NsNf[i]->GetYaxis()->CenterTitle(true);
	gPad->SetLogz();
	fh2_Califa_NsNf[i]->Draw("COLZ");

	TString buf3;
	if (i == 0) buf3 = "fh2_total_energy";
	if (i == 1) buf3 = "fh2_total_energy_withTofD";
	if (i == 2) buf3 = "fh2_total_energy_withTofDFoot";
	fh2_Califa_total_energy[i] = new TH1F(buf3, buf3, bins, minE, maxE);
	fh2_Califa_total_energy[i]->GetXaxis()->SetTitle("Energy [keV]");
	fh2_Califa_total_energy[i]->GetYaxis()->SetTitle("Counts");
	fh2_Califa_total_energy[i]->GetYaxis()->SetTitleOffset(1.4);
	fh2_Califa_total_energy[i]->GetXaxis()->CenterTitle(true);
	fh2_Califa_total_energy[i]->GetYaxis()->CenterTitle(true);
	fh2_Califa_total_energy[i]->SetFillColor(29);
	fh2_Califa_total_energy[i]->SetLineColor(1);
	fh2_Califa_total_energy[i]->SetLineWidth(2);
	gPad->SetLogy();
	fh2_Califa_total_energy[i]->Draw("");

	TString buf4;
	// find range for CrystalHits
	if (i == 0) buf4 = "fh2_CrystalHits";
	if (i == 1) buf4 = "fh2_CrystalHits_withTofD";
	if (i == 2) buf4 = "fh2_CrystalHits_withTofDFoot";
	fh2_Califa_CrystalHits[i] = new TH1F(buf4, buf4, 100, 0, 100);
	fh2_Califa_CrystalHits[i]->GetXaxis()->SetTitle("Crystal Hits");
	fh2_Califa_CrystalHits[i]->GetYaxis()->SetTitle("Counts");
	fh2_Califa_CrystalHits[i]->GetXaxis()->CenterTitle(true);
	fh2_Califa_CrystalHits[i]->GetYaxis()->CenterTitle(true);
	fh2_Califa_CrystalHits[i]->GetYaxis()->SetTitleOffset(1.2);
	fh2_Califa_CrystalHits[i]->GetXaxis()->SetTitleOffset(1.2);
	fh2_Califa_CrystalHits[i]->SetFillColor(8);
	fh2_Califa_CrystalHits[i]->SetLineColor(1);
	fh2_Califa_CrystalHits[i]->SetLineWidth(2);
	fh2_Califa_CrystalHits[i]->Draw("");

	TString buf5;
	// find range for ClusterId
	if (i == 0) buf5 = "fh2_ClusterId";
	if (i == 1) buf5 = "fh2_ClusterId_withTofD";
	if (i == 2) buf5 = "fh2_ClusterId_withTofDFoot";
	fh2_Califa_ClusterId[i] = new TH1F(buf5, buf5, 100, 0, 100);
	fh2_Califa_ClusterId[i]->GetXaxis()->SetTitle("Cluster Id");
	fh2_Califa_ClusterId[i]->GetYaxis()->SetTitle("Counts");
	fh2_Califa_ClusterId[i]->GetXaxis()->CenterTitle(true);
	fh2_Califa_ClusterId[i]->GetYaxis()->CenterTitle(true);
	fh2_Califa_ClusterId[i]->GetYaxis()->SetTitleOffset(1.2);
	fh2_Califa_ClusterId[i]->GetXaxis()->SetTitleOffset(1.2);
	fh2_Califa_ClusterId[i]->SetFillColor(8);
	fh2_Califa_ClusterId[i]->SetLineColor(1);
	fh2_Califa_ClusterId[i]->SetLineWidth(2);
	fh2_Califa_ClusterId[i]->Draw("");

	TString buf6;
	if (i == 0) buf6 = "fh2_theta_energy";
	if (i == 1) buf6 = "fh2_theta_energy_withTofD";
	if (i == 2) buf6 = "fh2_theta_energy_withTofDFoot";
	fh2_Califa_theta_energy[i] = new TH2F(buf6, buf6, 500, 0, 90, bins, minE, maxE);
	fh2_Califa_theta_energy[i]->GetXaxis()->SetTitle("Theta [degrees]");
	fh2_Califa_theta_energy[i]->GetYaxis()->SetTitle("Energy [keV]");
	fh2_Califa_theta_energy[i]->GetYaxis()->SetTitleOffset(1.4);
	fh2_Califa_theta_energy[i]->GetXaxis()->CenterTitle(true);
	fh2_Califa_theta_energy[i]->GetYaxis()->CenterTitle(true);
	fh2_Califa_theta_energy[i]->Draw("COLZ");

	TString buf7;
	if (i == 0) buf7 = "fh2_openangle_CalifaTofD";
	if (i == 1) buf7 = "fh2_openangle_withTofD";
	if (i == 2) buf7 = "fh2_openangle_withTofDFoot";
	fh2_openangle[i] = new TH1F(buf7, buf7, 160, 10, 170);
	fh2_openangle[i]->GetXaxis()->SetTitle("Opening angle [degrees]");
	fh2_openangle[i]->GetYaxis()->SetTitle("Counts");
        fh2_openangle[i]->GetXaxis()->CenterTitle(true);
	fh2_openangle[i]->GetYaxis()->CenterTitle(true);
	fh2_openangle[i]->GetYaxis()->SetTitleOffset(1.2);
	fh2_openangle[i]->GetXaxis()->SetTitleOffset(1.2);
	fh2_openangle[i]->SetFillColor(8);
	fh2_openangle[i]->SetLineColor(1);
	fh2_openangle[i]->SetLineWidth(2);
	fh2_openangle[i]->Draw("");

	TString buf8;
	if (i == 0) buf8 = "fh2_theta_corr";
	if (i == 1) buf8 = "fh2_theta_corr_withTofD";
	if (i == 2) buf8 = "fh2_theta_corr_withTofDFoot";
	fh2_coinTheta[i] = new TH2F(buf8, buf8, 500, 0, 100, 500, 0, 100);
	fh2_coinTheta[i]->GetXaxis()->SetTitle("Theta [degrees]");
	fh2_coinTheta[i]->GetYaxis()->SetTitle("Theta [degrees]");
	fh2_coinTheta[i]->GetYaxis()->SetTitleOffset(1.2);
	fh2_coinTheta[i]->GetXaxis()->CenterTitle(true);
	fh2_coinTheta[i]->GetYaxis()->CenterTitle(true);
	fh2_coinTheta[i]->Draw("COLZ");

	TString buf9;
	if (i == 0) buf9 = "fh2_phi_corr";
	if (i == 1) buf9 = "fh2_phi_corr_withTofD";
	if (i == 2) buf9 = "fh2_phi_corr_withTofDFoot";
	fh2_coinPhi[i] =
		new TH2F(buf9, buf9, 600, -190, 190, 600, -190, 190);
	fh2_coinPhi[i]->GetXaxis()->SetTitle("Phi [degrees]");
	fh2_coinPhi[i]->GetYaxis()->SetTitle("Phi [degrees]");
	fh2_coinPhi[i]->GetYaxis()->SetTitleOffset(1.2);
	fh2_coinPhi[i]->GetXaxis()->CenterTitle(true);
	fh2_coinPhi[i]->GetYaxis()->CenterTitle(true);
	fh2_coinPhi[i]->Draw("COLZ");

	TString buf10;
	if (i == 0) buf10 = "fh2_E_corr";
	if (i == 1) buf10 = "fh2_E_corr_withTofD";
	if (i == 2) buf10 = "fh2_E_corr_withTofDFoot";
	fh2_coinE[i] =
		new TH2F(buf10, buf10, bins, minE, maxE, bins, minE, maxE);
	fh2_coinE[i]->GetXaxis()->SetTitle("Energy (keV)");
	fh2_coinE[i]->GetYaxis()->SetTitle("Energy (keV)");
	fh2_coinE[i]->GetYaxis()->SetTitleOffset(1.2);
	fh2_coinE[i]->GetXaxis()->CenterTitle(true);
	fh2_coinE[i]->GetYaxis()->CenterTitle(true);
	fh2_coinE[i]->Draw("COLZ");

	TString buf11;
	if (i==0) buf11 = "fh2_leftE_openangle";
	if (i==1) buf11 = "fh2_leftE_openangle_withTofD";
	if (i==2) buf11 = "fh2_leftE_openangle_withTofDFoot";
	fh2_leftE_openangle[i] =
		new TH2F(buf11,buf11,160,10,170,bins,minE,maxE);
	fh2_leftE_openangle[i]->GetXaxis()->SetTitle("Open angle (degree)");
	fh2_leftE_openangle[i]->GetYaxis()->SetTitle("Energy (keV)");
	fh2_leftE_openangle[i]->GetYaxis()->SetTitleOffset(1.2);
	fh2_leftE_openangle[i]->GetXaxis()->CenterTitle(true);
	fh2_leftE_openangle[i]->GetYaxis()->CenterTitle(true);
	fh2_leftE_openangle[i]->Draw("COLZ");

	TString buf12;
	if (i==0) buf12 = "fh2_rightE_openangle";
	if (i==1) buf12 = "fh2_rightE_openangle_withTofD";
	if (i==2) buf12 = "fh2_rightE_openangle_withTofDFoot";
	fh2_rightE_openangle[i] =
		new TH2F(buf12,buf12,160,10,170,bins,minE,maxE);
	fh2_rightE_openangle[i]->GetXaxis()->SetTitle("Open angle (degree)");
	fh2_rightE_openangle[i]->GetYaxis()->SetTitle("Energy (keV)");
	fh2_rightE_openangle[i]->GetYaxis()->SetTitleOffset(1.2);
	fh2_rightE_openangle[i]->GetXaxis()->CenterTitle(true);
	fh2_rightE_openangle[i]->GetYaxis()->CenterTitle(true);
	fh2_rightE_openangle[i]->Draw("COLZ");
    }

    TString buf = "fh2_NsNf_withTofDFoot_openangle";
    fh2_Califa_NsNf[3] =
	    new TH2F(buf, buf, bins, minE, maxE, bins, minE, maxE);
    fh2_Califa_NsNf[3]->GetXaxis()->SetTitle("Nf Energy [keV]");
    fh2_Califa_NsNf[3]->GetYaxis()->SetTitle("Ns Energy [keV]");
    fh2_Califa_NsNf[3]->GetYaxis()->SetTitleOffset(1.4);
    fh2_Califa_NsNf[3]->GetXaxis()->CenterTitle(true);
    fh2_Califa_NsNf[3]->GetYaxis()->CenterTitle(true);
    gPad->SetLogz();
    fh2_Califa_NsNf[3]->Draw("COLZ");

    TString buf13 = "fh2_Q_tof";
    fh2_Q_tof = new TH2F(buf13,buf13,120,60,90,60,0,12);
    fh2_Q_tof->GetXaxis()->SetTitle("ToF");
    fh2_Q_tof->GetYaxis()->SetTitle("Q");
    fh2_Q_tof->GetYaxis()->SetTitleOffset(1.4);
    fh2_Q_tof->GetXaxis()->CenterTitle(true);
    fh2_Q_tof->GetYaxis()->CenterTitle(true);
    fh2_Q_tof->Draw("COLZ");

    // MAIN FOLDER-Califa
    TFolder* mainfolCalifa = new TFolder("CALIFAvsTofD", "CALIFA vs TofD info");

    if (fHitItemsCalifa && fHitItemsTofd)
    {
        //mainfolCalifa->Add(cCalifa_angles);
    }
    run->AddObject(mainfolCalifa);*/

    // Register command to reset histograms
    //run->GetHttpServer()->RegisterCommand("Reset_CalifavsTofD", Form("/Objects/%s/->Reset_Histo()", GetName()));

    std::cout << "CalifavsTofDOnlineSpectra end" << std::endl;
    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BCalifavsTofDOnlineSpectra::ReInit()
{
    std::cout << "CalifavsTofDOnlineSpectra ReInit" << std::endl;
    SetParContainers();
    // SetParameter();
    return kSUCCESS;
}

void R3BCalifavsTofDOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "");
    /*for (int i = 0; i < 2; i++)
    {
        fh2_Califa_theta_phi[i]->Reset();
	fh2_Califa_NsNf[i]->Reset();
	fh2_Califa_total_energy[i]->Reset();
	fh2_Califa_CrystalHits[i]->Reset();
	fh2_Califa_ClusterId[i]->Reset();
	fh2_Califa_theta_energy[i]->Reset();

	fh2_openangle[i]->Reset();
	fh2_coinTheta[i]->Reset();
	fh2_coinPhi[i]->Reset();
	fh2_coinE[i]->Reset();
    }

    fh2_Califa_NsNf[3]->Reset();
    fh2_Q_tof->Reset();*/
}

void R3BCalifavsTofDOnlineSpectra::Exec(Option_t* option)
{
    //if ((fTpat >= 0) && (header) && ((header->GetTpat() & fTpat) != fTpat))
        //return;

    tpatval = header->GetTpat();
    wrm = header->GetTimeStamp();
    //std::cout << "wrm: " << wrm << std::endl;
    //trigger = header->GetTrigger();

    map_crystalId.clear();
    map_energy.clear();
    map_Ns.clear();
    map_Nf.clear();
    map_febexTime.clear();
    map_wrts.clear();
    map_overflow.clear();
    map_pileup.clear();
    map_discard.clear();

    cal_crystalId.clear();
    cal_energy.clear();
    cal_Ns.clear();
    cal_Nf.clear();
    cal_time.clear();

    thetaList.clear();
    phiList.clear();
    energyList.clear();
    NfList.clear();
    NsList.clear();
    crystalHitsList.clear();
    clusterIdList.clear();
    timeList.clear();

    charge.clear();
    tof.clear();
    detectorId.clear();

    foot_detectorId.clear();
    foot_nbHits.clear();
    foot_mulStrip.clear();
    foot_pos.clear();
    foot_poslab.clear();
    foot_theta.clear();
    foot_phi.clear();
    foot_energy.clear();

    wr.clear();

    //int tofdHits = 0;
    //int califaHits = 0;
    int tofdHits = fHitItemsTofd->GetEntriesFast();
    int califaHits = fHitItemsCalifa->GetEntriesFast();
    loshits = fHitItemsLos->GetEntriesFast();
    frshits = fHitItemsFrs->GetEntriesFast();

    //std::cout << "los hits: " << losHits << std::endl;

    if ((tofdHits==0) && (califaHits==0))
    {
	    return;
    }

    if (fWRItemsCalifa && fWRItemsCalifa->GetEntriesFast()>0) {
    	for (Int_t ihit=0; ihit<fWRItemsCalifa->GetEntriesFast(); ihit++)
    	{
		R3BWRData* hit = (R3BWRData*) fWRItemsCalifa->At(ihit);
		if (!hit) continue;
		wr.push_back(hit->GetTimeStamp());
    	}
    }

    /*if (fHitItemsFoot && fHitItemsFoot->GetEntriesFast()>0) {
	for (Int_t ihit=0; ihit<fHitItemsFoot->GetEntriesFast(); ihit++)
	{
	    R3BFootHitData* hit = (R3BFootHitData*) fHitItemsFoot->At(ihit);
	    if (!hit) continue;

	    std::cout << "foot theta: " << hit->GetTheta() << std::endl;
	    foot_detectorId.push_back(hit->GetDetId());
	    foot_nbHits.push_back(hit->GetNbHit());
	    foot_mulStrip.push_back(hit->GetMulStrip());
	    foot_pos.push_back(hit->GetPos());
	    foot_poslab.push_back(hit->GetPosLab());
	    foot_theta.push_back(hit->GetTheta());
	    foot_phi.push_back(hit->GetPhi());
	    foot_energy.push_back(hit->GetEnergy());
	}
    }*/

    /*if (fWRItemsMaster && fWRItemsMaster->GetEntriesFast() > 0)
    {
	int WRhits = fWRItemsMaster->GetEntriesFast();
	for (Int_t ihit = 0; ihit < WRhits; ihit++)
	{
	    R3BWRData* hit = (R3BWRData*)fWRItemsMaster->At(ihit);
	    if (!hit)		
   		    continue;
	    wrm = hit->GetTimeStamp();	    
	}
    }*/

    for (Int_t ihit = 0; ihit < tofdHits; ihit++)
    {
        auto hit = (R3BTofdHitData*)fHitItemsTofd->At(ihit);
        if (!hit)
            continue;

	detectorId.push_back(hit->GetDetId());
	charge.push_back(hit->GetEloss());
	tof.push_back(hit->GetTof());
    }

    for (Int_t ihit = 0; ihit<fMapItemsCalifa->GetEntriesFast(); ihit++)
    {
        auto hit = (R3BCalifaMappedData*)fMapItemsCalifa->At(ihit);
	if (!hit)
	    continue;

	map_crystalId.push_back(hit->GetCrystalId());
	map_energy.push_back(hit->GetEnergy());
	map_Ns.push_back(hit->GetNs());
	map_Nf.push_back(hit->GetNf());
	map_febexTime.push_back(hit->GetFebexTime());
	map_wrts.push_back(hit->GetWrts());
	map_overflow.push_back(hit->GetOverFlow());
	map_pileup.push_back(hit->GetPileup());
	map_discard.push_back(hit->GetDiscard());
    }

    for (Int_t ihit =0; ihit<fCalItemsCalifa->GetEntriesFast(); ihit++)
    {
        auto hit = (R3BCalifaCrystalCalData*)fCalItemsCalifa->At(ihit);
	if (!hit)
	    continue;

	//std::cout << "ihit: " << ihit << std::endl;

	cal_crystalId.push_back(hit->GetCrystalId());
	cal_energy.push_back(hit->GetEnergy());
	cal_Ns.push_back(hit->GetNs());
	cal_Nf.push_back(hit->GetNf());
	cal_time.push_back(hit->GetTime());
    }

    Double_t califa_theta[califaHits];
    Double_t califa_phi[califaHits];
    Double_t califa_e[califaHits];
    Double_t califa_Ns[califaHits];
    Double_t califa_Nf[califaHits];

    for (Int_t ihit = 0; ihit < califaHits; ihit++)
    {
        auto hit = (R3BCalifaClusterData*)fHitItemsCalifa->At(ihit);
        if (hit->GetEnergy() < fMinProtonE)
            continue;

	thetaList.push_back(hit->GetTheta());
	phiList.push_back(hit->GetPhi());
	energyList.push_back(hit->GetEnergy());
	NsList.push_back(hit->GetNs());
	NfList.push_back(hit->GetNf());
	crystalHitsList.push_back(hit->GetNbOfCrystalHits());
	clusterIdList.push_back(hit->GetClusterId());
	timeList.push_back(hit->GetTime());
	
    }

        /*fh2_Califa_theta_phi[0]->Fill(theta, phi); // always
	fh2_Califa_NsNf[0]->Fill(hit->GetNs(), hit->GetNf());
	//fh2_Califa_total_energy[0]->Fill(hit->GetEnergy());
	fh2_Califa_CrystalHits[0]->Fill(hit->GetNbOfCrystalHits());
	fh2_Califa_ClusterId[0]->Fill(hit->GetClusterId());
	fh2_Califa_theta_energy[0]->Fill(theta, hit->GetEnergy());*/
	
        /*if (fZminus1) {
            fh2_Califa_theta_phi[1]->Fill(theta, phi); // only with TofD
	    fh2_Califa_NsNf[1]->Fill(hit->GetNs(), hit->GetNf());
	    fh2_Califa_total_energy[1]->Fill(hit->GetEnergy());
	    fh2_Califa_CrystalHits[1]->Fill(hit->GetNbOfCrystalHits());
	    fh2_Califa_ClusterId[1]->Fill(hit->GetClusterId());
	    fh2_Califa_theta_energy[1]->Fill(theta, hit->GetEnergy());

	    TVector3 hitdir;
	    hitdir.SetMagThetaPhi(1., hit->GetTheta(), hit->GetPhi());
	    if (isFootDetect(hitdir)) 
	    {
                fh2_Califa_theta_phi[2]->Fill(theta, phi);
		fh2_Califa_NsNf[2]->Fill(hit->GetNs(), hit->GetNf());
		fh2_Califa_total_energy[2]->Fill(hit->GetEnergy());
		fh2_Califa_CrystalHits[2]->Fill(hit->GetNbOfCrystalHits());
		fh2_Califa_ClusterId[2]->Fill(hit->GetClusterId());
		fh2_Califa_theta_energy[2]->Fill(theta, hit->GetEnergy());
	    }
	}*/

    // Hit data
    /*if (fHitItemsCalifa && fZminus1 && fHitItemsCalifa->GetEntriesFast() > 0)
    {

        Double_t theta = 0., phi = 0.;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hit = (R3BCalifaClusterData*)fHitItemsCalifa->At(ihit);
            if (!hit)
                continue;
            theta = hit->GetTheta() * TMath::RadToDeg();
            phi = hit->GetPhi() * TMath::RadToDeg();
            califa_theta[ihit] = theta;
            califa_phi[ihit] = phi;
            califa_e[ihit] = hit->GetEnergy();
        }

    	TVector3 master[2];
    	Double_t maxEL = 0., maxER = 0.;
   	Double_t NsL=0., NfL=0., NsR=0., NfR=0.;
    	for (Int_t i1 = 0; i1 < nHits; i1++)
    	{
      		if (califa_e[i1] > maxER && TMath::Abs(califa_phi[i1]) > 150.) // wixhausen
      		{
			master[0].SetMagThetaPhi(1., califa_theta[i1] * TMath::DegToRad(), califa_phi[i1] * TMath::DegToRad());
  			maxER = califa_e[i1];
			NsR = califa_Ns[i1];
			NfR = califa_Nf[i1];
      		}
      		if (califa_e[i1] > maxEL && TMath::Abs(califa_phi[i1]) < 60.)
      		{ // messel
        		master[1].SetMagThetaPhi(1., califa_theta[i1] * TMath::DegToRad(), califa_phi[i1] * TMath::DegToRad());
       			maxEL = califa_e[i1];
			NsL = califa_Ns[i1];
			NfL = califa_Nf[i1];
      		}
	}

    	if (maxEL > fMinProtonE && maxER > fMinProtonE)
    	{	
      		fh2_openangle[0]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg());
      		fh2_coinTheta[0]->Fill(master[0].Theta() * TMath::RadToDeg(), master[1].Theta() * TMath::RadToDeg());
      		fh2_coinPhi[0]->Fill(master[0].Phi() * TMath::RadToDeg(), master[1].Phi() * TMath::RadToDeg());
      		fh2_coinE[0]->Fill(maxER,maxEL);
      		fh2_leftE_openangle[0]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg(),maxEL);
      		fh2_rightE_openangle[0]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg(),maxER);

      		// temp mod
      		fh2_Califa_total_energy[0]->Fill(maxEL);
     		 fh2_Califa_total_energy[0]->Fill(maxER);

      		if (fZminus1) {
        		fh2_openangle[1]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg());
        		fh2_coinTheta[1]->Fill(master[0].Theta() * TMath::RadToDeg(), master[1].Theta() * TMath::RadToDeg());
			fh2_coinPhi[1]->Fill(master[0].Phi() * TMath::RadToDeg(), master[1].Phi() * TMath::RadToDeg());
			fh2_coinE[1]->Fill(maxER,maxEL);
			fh2_leftE_openangle[1]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg(),maxEL);
			fh2_rightE_openangle[1]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg(),maxER);

		if (isFootDetect(master[0]) && isFootDetect(master[1])) {
	  fh2_openangle[2]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg());
	  fh2_coinTheta[2]->Fill(master[0].Theta() * TMath::RadToDeg(), master[1].Theta() * TMath::RadToDeg());
	  fh2_coinPhi[2]->Fill(master[0].Phi() * TMath::RadToDeg(), master[1].Phi() * TMath::RadToDeg());
	  fh2_coinE[2]->Fill(maxER,maxEL);
	  fh2_leftE_openangle[2]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg(),maxEL);
	  fh2_rightE_openangle[2]->Fill(master[0].Angle(master[1]) * TMath::RadToDeg(),maxER);
	  
	  double openangle = master[0].Angle(master[1]) * TMath::RadToDeg();
	  double openangle_cut_low = 68.;
	  double openangle_cut_high = 88.;
	  //std::cout << "open angle: " << openangle << std::endl;
	  if ((openangle > openangle_cut_low) && (openangle < openangle_cut_high))
	  {
            //std::cout << "pass here" << std::endl;
            fh2_Califa_NsNf[3]->Fill(NsL,NfL);
	    fh2_Califa_NsNf[3]->Fill(NsR,NfR);
	  }
	}
      }
    }
    }*/

    fNEvents += 1;

    outtree->Fill();
}

void R3BCalifavsTofDOnlineSpectra::FinishEvent()
{
    if (fHitItemsCalifa)
    {
        fHitItemsCalifa->Clear();
    }
    if (fHitItemsTofd)
    {
        fHitItemsTofd->Clear();
    }
    if (fWRItemsCalifa)
    {
	fWRItemsCalifa->Clear();
    }
    if (fWRItemsMaster)
    {
	fWRItemsMaster->Clear();
    }
    if (fMapItemsCalifa)
    {
	fMapItemsCalifa->Clear();
    }
    if (fCalItemsCalifa)
    {
	fCalItemsCalifa->Clear();
    }
    if (fHitItemsLos)
    {
	fHitItemsLos->Clear();
    }
    if (fHitItemsFrs)
    {
	fHitItemsFrs->Clear();
    }
    if (fHitItemsFoot)
    {
	fHitItemsFoot->Clear();
    }
}

void R3BCalifavsTofDOnlineSpectra::FinishTask()
{
    // Write canvas for Hit data
    /*if (fHitItemsCalifa)
    {
        fh2_Califa_theta_phi[0]->Write();
	fh2_Califa_NsNf[0]->Write();
	fh2_Califa_total_energy[0]->Write();
	fh2_Califa_CrystalHits[0]->Write();
	fh2_Califa_ClusterId[0]->Write();
	fh2_Califa_theta_energy[0]->Write();

	fh2_openangle[0]->Write();
	fh2_coinTheta[0]->Write();
        fh2_coinPhi[0]->Write();
	fh2_coinE[0]->Write();
	fh2_leftE_openangle[0]->Write();
	fh2_rightE_openangle[0]->Write();

	if (fHitItemsTofd) 
	{
	    fh2_Q_tof->Write();

            fh2_Califa_theta_phi[1]->Write();
	    fh2_Califa_theta_phi[2]->Write();
	    fh2_Califa_NsNf[1]->Write();
	    fh2_Califa_NsNf[2]->Write();
	    fh2_Califa_total_energy[1]->Write();
	    fh2_Califa_total_energy[2]->Write();
	    fh2_Califa_CrystalHits[1]->Write();
	    fh2_Califa_CrystalHits[2]->Write();
	    fh2_Califa_ClusterId[1]->Write();
	    fh2_Califa_ClusterId[2]->Write();
	    fh2_Califa_theta_energy[1]->Write();
	    fh2_Califa_theta_energy[2]->Write();

    	    fh2_openangle[1]->Write();
	    fh2_openangle[2]->Write();
	    fh2_coinTheta[1]->Write();
	    fh2_coinTheta[2]->Write();
	    fh2_coinPhi[1]->Write();
	    fh2_coinPhi[2]->Write();
	    fh2_coinE[1]->Write();
	    fh2_coinE[2]->Write();
	    fh2_leftE_openangle[1]->Write();
	    fh2_leftE_openangle[2]->Write();
	    fh2_rightE_openangle[1]->Write();
	    fh2_rightE_openangle[2]->Write();
	    fh2_Califa_NsNf[3]->Write();
	}
    }*/

    outtree->SetName("genT");
    outfile->Write();
    outfile->Close();
}

ClassImp(R3BCalifavsTofDOnlineSpectra);
