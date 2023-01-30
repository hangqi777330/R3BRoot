/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVector3.h"
#include "TPolyMarker.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BCalifaCrystalCalPar.h"
#include "R3BCalifaMapped2CrystalCalPar.h"
#include "R3BCalifaMappedData.h"
#include "R3BCalifaMappingPar.h"

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>

R3BCalifaMapped2CrystalCalPar::R3BCalifaMapped2CrystalCalPar()
    : R3BCalifaMapped2CrystalCalPar("R3B CALIFA Calibration Parameters Finder ", 1)
{
}

R3BCalifaMapped2CrystalCalPar::R3BCalifaMapped2CrystalCalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMap_Par(NULL)
    , fCal_Par(NULL)
    , fCalifaMappedDataCA(NULL)
    , fNumCrystals(1)
    , fNumParam(0)
    , fMinStadistics(100)
    , fMapHistos_left(0)
    , fMapHistos_right(0)
    , fMapHistos_bins(0)
    , fMapHistos_leftp(0)
    , fMapHistos_rightp(0)
    , fMapHistos_binsp(0)
    , fNumPeaks(0)
    , fSigma(0)
    , fThreshold(0)
    , fEnergyPeaks(NULL)
    , fDebugMode(0)
    , fSourceName("fitting")
{
}

R3BCalifaMapped2CrystalCalPar::~R3BCalifaMapped2CrystalCalPar()
{
    LOG(INFO) << "R3BCalifaMapped2CrystalCalPar: Delete instance";
    if (fCalifaMappedDataCA)
        delete fCalifaMappedDataCA;
    if (fEnergyPeaks)
        delete fEnergyPeaks;
}

void R3BCalifaMapped2CrystalCalPar::SetParContainers()
{
    // Parameter Container
    // Reading califaMappingPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(ERROR) << "FairRuntimeDb not opened!";
    }

    fMap_Par = (R3BCalifaMappingPar*)rtdb->getContainer("califaMappingPar");
    if (!fMap_Par)
    {
        LOG(ERROR) << "R3BCalifaMapped2CrystalCalPar::Init() Couldn't get handle on califaMappingPar container";
    }
    else
    {
        LOG(INFO) << "R3BCalifaMapped2CrystalCalPar:: califaMappingPar container open";
    }
}

void R3BCalifaMapped2CrystalCalPar::SetParameter()
{
    if (!fMap_Par)
    {
        LOG(WARNING) << "R3BCalifaMapped2CrystalCalPar::Container califaMappingPar not found.";
    }
    //--- Parameter Container ---
    fNumCrystals = fMap_Par->GetNumCrystals(); // Number of crystals x 2
    LOG(INFO) << "R3BCalifaMapped2CrystalCalPar::NumCry " << fNumCrystals;
    // fMap_Par->printParams();
}

InitStatus R3BCalifaMapped2CrystalCalPar::Init()
{
    LOG(INFO) << "R3BCalifaMapped2CrystalCalPar::Init()";

    if (!fEnergyPeaks)
    {
        fEnergyPeaks = new TArrayF;
        fEnergyPeaks->Set(fNumPeaks);
    }

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(ERROR) << "R3BCalifaMapped2CrystalCalPar::Init() FairRootManager not found";
        return kFATAL;
    }

    fCalifaMappedDataCA = (TClonesArray*)rootManager->GetObject("CalifaMappedData");
    if (!fCalifaMappedDataCA)
    {
        LOG(ERROR) << "R3BCalifaMapped2CrystalCalPar::Init() CalifaMappedData not found";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(ERROR) << "R3BCalifaMapped2CrystalCalPar::Init() FairRuntimeDb not found";
        return kFATAL;
    }

    fCal_Par = (R3BCalifaCrystalCalPar*)rtdb->getContainer("califaCrystalCalPar");
    if (!fCal_Par)
    {
        LOG(ERROR) << "R3BCalifaMapped2CrystalCalPar::Init() Couldn't get handle on califaCrystalCalPar container";
        return kFATAL;
    }

    // Initiate output file
    outrootfile = TFile::Open("spectrum.root","UPDATE");
    if (!outrootfile) {
    	outrootfile = new TFile("spectrum.root","RECREATE");
	outrootfile->cd();
	outroottree = new TTree("genT","General Tree");
    }

    // Set container with mapping parameters
    SetParameter();

    // Create histograms for crystal calibration
    char name1[100];
    char name2[100];
    char name3[100];
    char name4[100];
    char name5[100];
    Int_t fright, fleft, fbins;
    if (fSourceName=="60Co" || fSourceName=="22Na" || fSourceName=="AmBe")
    {
        fh_Map_energy_crystal = new TH1F*[fNumCrystals];
    } else if (fSourceName == "fitting")
    {
        fh2_peak_cal = new TH2F*[fNumCrystals];
	fh2_sig_crystal = new TH2F*[3]; // 3 sources
	//fh_Map_nobkg = new TH1F*[3*fNumCrystals];
    }
    for (Int_t i = 0; i < fNumCrystals; i++)
        if (fMap_Par->GetInUse(i + 1) == 1)
        {

            sprintf(name1, "fh_Map_energy_crystal_"+fSourceName+"_%i", i + 1);
	    sprintf(name2, "fh2_peak_cal_%i", i + 1);
	    sprintf(name3, "fh_map_crystal_nobkg_22Na_%i", i + 1);
	    sprintf(name4, "fh_map_crystal_nobkg_60Co_%i", i + 1);
	    sprintf(name5, "fh_map_crystal_nobkg_AmBe_%i", i + 1);
            if (i < fMap_Par->GetNumCrystals() / 2)
            {
                fright = fMapHistos_right;
                fleft = fMapHistos_left;
                fbins = fMapHistos_bins;
            }
            else
            {
                fright = fMapHistos_rightp;
                fleft = fMapHistos_leftp;
                fbins = fMapHistos_binsp;
            }
	    if (fSourceName=="60Co" || fSourceName=="22Na" || fSourceName=="AmBe")
	    {
                fh_Map_energy_crystal[i] = new TH1F(name1, name1, fbins, fleft, fright);
	    } else if (fSourceName=="fitting")
	    {
	        fh2_peak_cal[i] = new TH2F(name2, name2, fbins, fleft, fright, 250, 0, 5000);
		/*fh_Map_nobkg[i] = new TH1F(name3,name3,fbins,fleft,fright);
		fh_Map_nobkg[fNumCrystals+i] = new TH1F(name4,name4,fbins,fleft,fright);
		fh_Map_nobkg[2*fNumCrystals+i] = new TH1F(name5,name5,fbins,fleft,fright);*/
	    }
        }

    if (fSourceName=="60Co" || fSourceName=="22Na" || fSourceName=="AmBe")
    {
        fh2_Map_crystal_gamma = new TH2F("fh2_Map_crystal_gamma","fh2_Map_crystal_gamma;crystal ID;Map Energy",fNumCrystals/2,0,fNumCrystals/2,fMapHistos_bins,fMapHistos_left,fMapHistos_right);
        fh2_Map_crystal_proton = new TH2F("fh2_Map_crystal_proton","fh2_Map_crystal_proton;crystal ID;Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,fMapHistos_binsp,fMapHistos_leftp,fMapHistos_rightp);
        fh_peak_crystal_gamma = new TH2F("fh_peak_crystal_gamma","fh_peak_crystal_gamma;crystal ID;peak Map Energy",fNumCrystals/2,0,fNumCrystals/2,fMapHistos_bins,fMapHistos_left,fMapHistos_right);
        fh_peak_crystal_proton = new TH2F("fh_peak_crystal_proton","fh_peak_crystal_proton;crystal ID;peak Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,fMapHistos_binsp,fMapHistos_leftp,fMapHistos_rightp);
        fh_sigma_crystal_gamma = new TH2F("fh_sigma_crystal_gamma","fh_sigma_crystal_gamma;crystal ID;sigma Map Energy",fNumCrystals/2,0,fNumCrystals/2,10*fSigma,0,10*fSigma);
        fh_sigma_crystal_proton = new TH2F("fh_sigma_crystal_proton","fh_sigma_crystal_proton;crystal ID;sigma Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,10*fSigma,0,10*fSigma);
    } else if (fSourceName == "fitting")
    {
        fh2_slope_crystalID = new TH2F("fh2_slope_crystalID",";crystal ID;slope",fNumCrystals,0,fNumCrystals,40,0,20);
	fh2_sig_crystal[0] = new TH2F("fh2_sig_crystal_22Na","22Na;crystal ID;significance",fNumCrystals,0,fNumCrystals,40,0,400);
	fh2_sig_crystal[1] = new TH2F("fh2_sig_crystal_60Co","60Co;crystal ID;significance",fNumCrystals,0,fNumCrystals,40,0,400);
	fh2_sig_crystal[2] = new TH2F("fh2_sig_crystal_AmBe","AmBe;crystal ID;significance",fNumCrystals,0,fNumCrystals,20,0,100);
	fh2_chi2_crystal = new TH2F("fh2_chi2_crystal",";crystal ID;chi2",fNumCrystals,0,fNumCrystals,400,0,40000);
	fh_numPeak = new TH1F("fh_numPeak","fh_numPeak",10,0,10);
    }

    return kSUCCESS;
}

InitStatus R3BCalifaMapped2CrystalCalPar::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BCalifaMapped2CrystalCalPar::Exec(Option_t* opt)
{
    if (fSourceName == "fitting") return;

    Int_t nHits = fCalifaMappedDataCA->GetEntries();
    if (!nHits)
        return;

    R3BCalifaMappedData** MapHit = new R3BCalifaMappedData*[nHits];
    Int_t crystalId = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        MapHit[i] = (R3BCalifaMappedData*)(fCalifaMappedDataCA->At(i));
        crystalId = MapHit[i]->GetCrystalId();
        // Fill histograms
        if (fMap_Par->GetInUse(crystalId) == 1) {
	    Double_t fleft, fright;
	    if (crystalId<=fNumCrystals/2)
	    {
		fleft = fMapHistos_left;
		fright = fMapHistos_right;
	    } else
	    {
		fleft = fMapHistos_leftp;
		fright = fMapHistos_rightp;
	    }

	    if (MapHit[i]->GetEnergy()>=fleft and MapHit[i]->GetEnergy()<=fright)
	    {
                fh_Map_energy_crystal[crystalId - 1]->Fill(MapHit[i]->GetEnergy());
	        if (crystalId < fNumCrystals/2)
	        {
	            fh2_Map_crystal_gamma->Fill(crystalId-1,MapHit[i]->GetEnergy());
	        } else 
	        {
		    fh2_Map_crystal_proton->Fill(crystalId-1,MapHit[i]->GetEnergy());
		}
	    }
	}
    }

    if (MapHit)
        delete MapHit;
    return;
}

void R3BCalifaMapped2CrystalCalPar::Reset() {}

void R3BCalifaMapped2CrystalCalPar::FinishEvent() {}

void R3BCalifaMapped2CrystalCalPar::FinishTask()
{
    if (fSourceName == "22Na" || fSourceName == "60Co" || fSourceName == "AmBe")
    {
	SearchPeaks();
    }

    if (fSourceName == "fitting")
    {	
	FitPeaks();
    	fCal_Par->printParams();
    }

    outrootfile->Write();
    outrootfile->Close();
}

// find and record
// mapping level peaks
void R3BCalifaMapped2CrystalCalPar::SearchPeaks()
{
    for (Int_t i=0; i<fNumCrystals; i++)
    {
	Int_t nfound = 0;
	TSpectrum* ss = new TSpectrum(fNumPeaks);

	if ((fMap_Par->GetInUse(i+1) == 1) && (fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics))
	{
	    nfound = ss->Search(fh_Map_energy_crystal[i], fSigma, "", fThreshold); // "goff" to turn off drawing
            //fh_Map_energy_crystal[i]->Write();	    
        
            fChannelPeaks = (Double_t*) ss->GetPositionX();
	    Int_t idx[nfound];
            TMath::Sort(nfound, fChannelPeaks, idx, kTRUE);

	    for (Int_t j=0; j<nfound; j++)
	    {
		// gaussian fit
		Double_t posX = fChannelPeaks[idx[nfound-j-1]];
		TF1* gaussfit = new TF1("gaussfit","gaus",posX-fSigma,posX+fSigma);
		TH1F* h_copy = (TH1F*) fh_Map_energy_crystal[i]->Clone("h_copy");
		h_copy->Fit("gaussfit","RQ");
		Double_t mean = gaussfit->GetParameter(1);
		Double_t sigma = gaussfit->GetParameter(2);
	
		Double_t pmX[1] = {(Double_t) i};
		Double_t pmY[1] = {mean};
		Double_t pmZ[1] = {sigma};
		TPolyMarker *pm1 = new TPolyMarker(1,pmX,pmY);
		pm1->SetMarkerStyle(23);
		pm1->SetMarkerColor(kRed);
		TPolyMarker *pm2 = new TPolyMarker(1,pmX,pmZ);
		pm2->SetMarkerStyle(23);
		pm2->SetMarkerColor(kBlack);

		if (i<fNumCrystals/2)
		{
			//fh_peak_crystal_gamma->Fill(i,mean);
			//fh_sigma_crystal_gamma->Fill(i,sigma);
			fh2_Map_crystal_gamma->GetListOfFunctions()->Add(pm1);
			fh2_Map_crystal_gamma->GetListOfFunctions()->Print();
			fh_peak_crystal_gamma->GetListOfFunctions()->Add(pm1);
			fh_peak_crystal_gamma->GetListOfFunctions()->Print();
			fh_sigma_crystal_gamma->GetListOfFunctions()->Add(pm2);
			fh_sigma_crystal_gamma->GetListOfFunctions()->Print();
		} else
		{
			//fh_peak_crystal_proton->Fill(i,mean);
			//fh_sigma_crystal_proton->Fill(i,sigma);
			fh2_Map_crystal_proton->GetListOfFunctions()->Add(pm1);
			fh2_Map_crystal_proton->GetListOfFunctions()->Print();
			fh_peak_crystal_proton->GetListOfFunctions()->Add(pm1);
			fh_peak_crystal_proton->GetListOfFunctions()->Print();
			fh_sigma_crystal_proton->GetListOfFunctions()->Add(pm2);
			fh_sigma_crystal_proton->GetListOfFunctions()->Print();
		}
	    }
	}

	delete ss;
    }

    /*fh2_Map_crystal_gamma->Write("colz");
    fh2_Map_crystal_proton->Write("colz");
    fh_peak_crystal_gamma->Write();
    fh_peak_crystal_proton->Write();
    fh_sigma_crystal_gamma->Write();
    fh_sigma_crystal_gamma->Write();
    fh_sigma_crystal_proton->Write();*/
}

void R3BCalifaMapped2CrystalCalPar::FitPeaks()
{
    Double_t fSigThreshold = 3.0;
    Double_t fChi2Threshold = 20000.0;

    fCal_Par->SetNumCrystals(fNumCrystals);
    fCal_Par->SetNumParametersFit(fNumParam);
    fCal_Par->GetCryCalParams()->Set(fNumParam * fNumCrystals);
    Int_t maxPeaks = 20;
    Double_t X[fNumCrystals][maxPeaks];
    Double_t Y[fNumCrystals][maxPeaks];
    Double_t eX[fNumCrystals][maxPeaks];
    Double_t eY[fNumCrystals][maxPeaks];
    // inference
    TH1F* AmBeSpectra[fNumCrystals];

    for (Int_t i=0; i<fNumCrystals; i++)
    {
	for (Int_t j=0; j<maxPeaks; j++)
	{
	    X[i][j] = -1;
	    Y[i][j] = -1;
	    eX[i][j] = -1;
	    eY[i][j] = -1;
	}
    }

    for (auto k : *outrootfile->GetListOfKeys())
    {
	TKey *key = static_cast<TKey*>(k);
	TClass *cl = gROOT->GetClass(key->GetClassName());
	std::string title = key->GetName();
	std::string begintitle = title.substr(0,13);

        if (begintitle != "fh_Map_energy") continue;
	if (!(cl->InheritsFrom("TH1"))) continue;

	std::string sourceName = title.substr(22,4);
	std::string crystalName = title.substr(27,4);
	Int_t cryId = std::atoi(crystalName.c_str());

	fEnergyPeaks->Reset();

        if (sourceName == "22Na")
	{
	  fNumPeaks = 2;
	  fEnergyPeaks->Set(2);
	  fEnergyPeaks->AddAt(1274.5,0);
	  fEnergyPeaks->AddAt(511.0,1);
	} else if (sourceName == "60Co") 
	{
	  fNumPeaks = 2;
	  fEnergyPeaks->Set(2);
	  fEnergyPeaks->AddAt(1332.5,0);
	  fEnergyPeaks->AddAt(1173.2,1);
	} else if (sourceName == "AmBe")
	{
	  fNumPeaks = 3;
	  fEnergyPeaks->Set(3);
	  fEnergyPeaks->AddAt(4438.0,0);
	  fEnergyPeaks->AddAt(3927.0,1);
	  fEnergyPeaks->AddAt(3416.0,2);

	  AmBeSpectra[cryId-1] = (TH1F*) key->ReadObject<TH1F>();
	}

	Int_t nfound;
	TH1F* h = key->ReadObject<TH1F>();
	TSpectrum* ss = new TSpectrum(fNumPeaks);
	if ((fMap_Par->GetInUse(cryId) == 1) && (h->GetEntries() > fMinStadistics))
	{
	    TH1F *h2 = (TH1F*)h->Clone("h2");
	    nfound = ss->Search(h2,fSigma,"goff",fThreshold);
	    TH1F* hbkg = (TH1F*)ss->Background(h2);
	    fChannelPeaks = (Double_t*) ss->GetPositionX();
	    Int_t idx[nfound];
	    TMath::Sort(nfound, fChannelPeaks, idx, kTRUE);

	    Double_t tempX[fNumPeaks];
	    Double_t tempY[fNumPeaks];
	    Double_t tempS[fNumPeaks];

	    for (Int_t tempi=0; tempi<fNumPeaks; tempi++)
	    {
		tempX[tempi] = -1;
		tempY[tempi] = -1;
		tempS[tempi] = -1;
	    }

	    if (nfound != fNumPeaks) continue;

	    if (fNumPeaks==2)
	    {
		TH1F* h3 = (TH1F*) h2->Clone("h3");
		Double_t posX1 = fChannelPeaks[idx[1]];
		Double_t posX2 = fChannelPeaks[idx[0]];

		TF1* fFit = new TF1("fFit","exp([0]+[1]*x)+[2]*exp(-0.5*((x-[3])/[4])^2)+[5]*exp(-0.5*((x-[6])/[7])^2)",h3->GetXaxis()->GetXmin(),h3->GetXaxis()->GetXmax());
		TF1* fExpo = new TF1("fExpo","expo",hbkg->GetXaxis()->GetXmin(),hbkg->GetXaxis()->GetXmax());
		// axis range within limit
                Double_t gaus1min = TMath::Max(posX1-fSigma,h3->GetXaxis()->GetXmin());
		Double_t gaus1max = TMath::Min(posX1+fSigma,h3->GetXaxis()->GetXmax());
		Double_t gaus2min = TMath::Max(posX2-fSigma,h3->GetXaxis()->GetXmin());
		Double_t gaus2max = TMath::Min(posX2+fSigma,h3->GetXaxis()->GetXmax());
		TF1* fGaus1 = new TF1("fGaus1","gaus",gaus1min,gaus1max);
		TF1* fGaus2 = new TF1("fGaus2","gaus",gaus2min,gaus2max);
		hbkg->Fit("fExpo","RQ");
		h3->Fit("fGaus1","RQ");
		h3->Fit("fGaus2","RQ+");
		fFit->SetParameter(0,fExpo->GetParameter(0));
		fFit->SetParameter(1,fExpo->GetParameter(1));
		fFit->SetParameter(2,fGaus1->GetParameter(0));
		fFit->SetParameter(3,fGaus1->GetParameter(1));
		fFit->SetParameter(4,fGaus1->GetParameter(2));
		fFit->SetParameter(5,fGaus2->GetParameter(0));
		fFit->SetParameter(6,fGaus2->GetParameter(1));
		fFit->SetParameter(7,fGaus2->GetParameter(2));
		h3->Fit("fFit","RQ+");

		tempX[0] = fFit->GetParameter(3);
		tempY[0] = fEnergyPeaks->GetAt(1);
		tempS[0] = fFit->GetParameter(4);
		tempX[1] = fFit->GetParameter(6);
		tempY[1] = fEnergyPeaks->GetAt(0);
		tempS[1] = fFit->GetParameter(7);
	    }
	    
	    int start = 0;
	    while (X[cryId-1][start]>0.) start++;

	    for (Int_t j=0; j<fNumPeaks; j++)
	    {
		X[cryId-1][start+j] = tempX[j];
		Y[cryId-1][start+j] = tempY[j];
		eX[cryId-1][start+j] = tempS[j];
		eY[cryId-1][start+j] = 0.;
	    }
	}
    }

    // Fit with graph
    Int_t fleft, fright;
    for (Int_t i=0; i<fNumCrystals; i++)
    {
	if (fMap_Par->GetInUse(i+1) == 1)
	{
	    if (i<fMap_Par->GetNumCrystals()/2)
	    {
		fright = fMapHistos_right;
		fleft = fMapHistos_left;
	    } else 
	    {
		fright = fMapHistos_rightp;
		fleft = fMapHistos_leftp;
	    }

	    if (fNumParam == 2)
	    {
		f1 = new TF1("f1","[0]+[1]*x",fleft,fright);
	    }

	    Int_t numPeak = 0;

	    while (X[i][numPeak]>0.) 
	    {
		    numPeak++;
	    }

	    if (numPeak >= 2)
	    {
		// fit for function
		auto graph = new TGraphErrors(numPeak,X[i],Y[i],eX[i],eY[i]);
		Int_t crystalID = i+1;
		graph->Fit("f1","WQN");
		graph->Fit("f1","MQ+");

		Double_t chi2 = graph->Chisquare(f1,"");
		if (chi2/numPeak < fChi2Threshold && AmBeSpectra[i])
		{
		  Double_t intercept = f1->GetParameter(0); // y position of f1 at 0
		  Double_t slope = f1->GetParameter(1);
		  
	          Double_t AmBePeaks[3] = {3416.0, 3927.0, 4438.0};

	          TH1F* hPeak = (TH1F*) AmBeSpectra[i]->Clone("hPeak");
                  TSpectrum* ss1 = new TSpectrum(3);
		  TH1F* hbkg = (TH1F*) ss1->Background(hPeak);
		  hPeak->Add(hbkg,-1);

		  for (Int_t t=0; t<3; t++)
		  {
	              Double_t expPeak = (AmBePeaks[t]-intercept)/slope;
		      Double_t gaus3min = TMath::Max(expPeak-2*fSigma,hPeak->GetXaxis()->GetXmin());
		      Double_t gaus3max = TMath::Min(expPeak+2*fSigma,hPeak->GetXaxis()->GetXmax());
		      TF1 *fGaus3 = new TF1("fGaus3","gaus",gaus3min,gaus3max);
		      hPeak->Fit(fGaus3,"RQ+");
		    
		      Int_t posBin = (fGaus3->GetParameter(1) - hPeak->GetXaxis()->GetXmin())/hPeak->GetXaxis()->GetBinWidth(0);
		      Double_t bkgY = hbkg->GetBinContent(posBin);
		      Double_t sigY = hPeak->GetBinContent(posBin);
		      Double_t significance = sigY/TMath::Sqrt(bkgY);
		      if (significance > fSigThreshold)
		      {
			  X[i][numPeak] = fGaus3->GetParameter(1);
			  eX[i][numPeak] = fGaus3->GetParameter(2);
			  Y[i][numPeak] = AmBePeaks[t];
			  eY[i][numPeak] = 0.;
			  numPeak++;
		      }
		  }

		  auto graph1 = new TGraphErrors(numPeak,X[i],Y[i],eX[i],eY[i]);
		  graph1->Fit("f1","WQN");
		  graph1->Fit("f1","MQ+");
		  chi2 = graph1->Chisquare(f1,"");
		}

		if (i > fNumCrystals/2) // reset problematic crystal
		{
		    if (f1->GetParameter(1)<10. or f1->GetParameter(1)>16.)
		    {
			f1->SetParameter(0,0.);
			f1->SetParameter(1,0.);
		    }
		} else 
		{
		    if (f1->GetParameter(1)<1. or f1->GetParameter(1)>1.6)
		    {
			f1->SetParameter(0,0.);
			f1->SetParameter(1,0.);
		    }
		}

		// write histograms
		TPolyMarker *pm1 = new TPolyMarker(numPeak,X[i],Y[i]);
		pm1->SetMarkerStyle(23);
		pm1->SetMarkerColor(kBlack);
		for (Int_t peakIdx=0; peakIdx<numPeak; peakIdx++)
		{
		    //fh2_peak_cal[i]->Fill(X[i][peakIdx],Y[i][peakIdx]);
		    fh2_peak_cal[i]->GetListOfFunctions()->Add(pm1);
		}
		fh2_peak_cal[i]->GetListOfFunctions()->Add(f1);
		fh2_peak_cal[i]->SetStats(0);
		fh2_peak_cal[i]->Write();

		// Fill slope histogram
		Double_t pmX[1] = {(Double_t) crystalID};
		Double_t pmY[1] = {f1->GetParameter(1)};
		TPolyMarker *pm2 = new TPolyMarker(1,pmX,pmY);
		pm2->SetMarkerStyle(23);
		pm2->SetMarkerColor(kBlack);
                //fh2_slope_crystalID->Fill(crystalID,f1->GetParameter(1));
		fh2_slope_crystalID->GetListOfFunctions()->Add(pm2);

		Double_t pm4X[1] = {(Double_t) crystalID};
		Double_t pm4Y[1] = {chi2/numPeak};
		TPolyMarker *pm4 = new TPolyMarker(1,pm4X,pm4Y);
		pm4->SetMarkerStyle(23);
		pm4->SetMarkerColor(kRed);
		fh2_chi2_crystal->GetListOfFunctions()->Add(pm4);
		//fh2_chi2_crystal->Fill(chi2/numPeak,crystalID);

		for (Int_t h=0; h<fNumParam; h++)
		{
		    fCal_Par->SetCryCalParams(f1->GetParameter(h), fNumParam*i+h);
		}
	    }

            fh_numPeak->Fill(numPeak);
	}
    }

    fh2_slope_crystalID->Write();
    fh2_sig_crystal[0]->Write();
    fh2_sig_crystal[1]->Write();
    fh2_sig_crystal[2]->Write();
    fh2_chi2_crystal->Write();
    fh_numPeak->Write();

    fCal_Par->setChanged();
    return;
}

/*void R3BCalifaMapped2CrystalCalPar::SearchPeaks()
{
    Int_t nfound = 0;
    Int_t numPars = 2; // Number of parameters=2 by default
    if (fNumParam)
    {
        numPars = fNumParam;
    }

    fCal_Par->SetNumCrystals(fNumCrystals);
    fCal_Par->SetNumParametersFit(fNumParam);
    fCal_Par->GetCryCalParams()->Set(numPars * fNumCrystals);

    TSpectrum* ss = new TSpectrum(fNumPeaks);

    Int_t fright, fleft;

    for (Int_t i = 0; i < fNumCrystals; i++)
        if (fMap_Par->GetInUse(i + 1) == 1)
        {

            if (fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics)
            {

                if (fDebugMode)
                    nfound = ss->Search(fh_Map_energy_crystal[i], fSigma, "", fThreshold); // number of peaks
                else
                    nfound = ss->Search(fh_Map_energy_crystal[i], fSigma, "goff", fThreshold);
                
		LOG(DEBUG) << "CrystalId=" << i + 1 << " " << nfound << " " << fThreshold;
                fChannelPeaks = (Double_t*)ss->GetPositionX();

                Int_t idx[nfound];
                TMath::Sort(nfound, fChannelPeaks, idx, kTRUE);

                // Calibrated Spectrum
                Double_t X[nfound + 1];
                Double_t Y[nfound + 1];

		if (nfound != 2) {
			std::cout << "crystal Id: " << i+1 << std::endl;
		}

                for (Int_t j = 0; j < nfound; j++)
                {
                    X[j] = fChannelPeaks[idx[nfound - j - 1]];
                    Y[j] = fEnergyPeaks->GetAt(nfound - j - 1);

                    LOG(DEBUG) << "CrystalId=" << i + 1 << " " << j + 1 << " " << X[j + 1];
                }
                X[nfound] = 0.;
                Y[nfound] = 0.;

                if (i < fMap_Par->GetNumCrystals() / 2)
                {
                    fright = fMapHistos_right;
                    fleft = fMapHistos_left;
                }
                else
                {
                    fright = fMapHistos_rightp;
                    fleft = fMapHistos_leftp;
                }

                TF1* f1 = nullptr;
                if (fNumParam)
                {

                    if (fNumParam == 1)
                    {
                        f1 = new TF1("f1", "[0]*x", fleft, fright);
                    }
                    if (fNumParam == 2)
                    {
                        f1 = new TF1("f1", "[0]+[1]*x", fleft, fright);
                    }
                    if (fNumParam == 3)
                    {
                        f1 = new TF1("f1", "[0]+[1]*x+[2]*pow(x,2)", fleft, fright);
                    }
                    if (fNumParam == 4)
                    {
                        f1 = new TF1("f1", "[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)", fleft, fright);
                    }
                    if (fNumParam == 5)
                    {
                        f1 = new TF1("f1", "[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)", fleft, fright);
                    }
                    if (fNumParam > 5)
                    {
                        LOG(WARNING)
                            << "R3BCalifaMapped2CrystalCalPar:: The number of fit parameters can not be higher than 5";
                    }
                }
                else
                {
                    LOG(WARNING)
                        << "R3BCalifaMapped2CrystalCalPar:: No imput number of fit parameters, therefore, by default "
                           "NumberParameters=2";
                    f1 = new TF1("f1", "[0]+[1]*x", fleft, fright);
                }

                TGraph* graph = new TGraph(fNumPeaks + 1, X, Y);
                graph->Fit("f1", "Q"); // Quiet mode (minimum printing)
		graph->Draw("");

		std::cout << "crystalId: " << i+1 << std::endl;
                for (Int_t h = 0; h < numPars; h++)
                {
                    fCal_Par->SetCryCalParams(f1->GetParameter(h), numPars * i + h);
                }
            }
            else
            {
                LOG(WARNING) << "R3BCalifaMapped2CrystalCalPar::Histogram number " << i + 1 << "not Fitted";
            }
        }

    delete ss;
    fCal_Par->setChanged();
    return;
}*/

ClassImp(R3BCalifaMapped2CrystalCalPar)
