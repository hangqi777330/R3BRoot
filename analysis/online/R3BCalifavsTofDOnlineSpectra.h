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

// ------------------------------------------------------------
// -----            R3BCalifavsTofDOnlineSpectra          -----
// -----    Created 21/05/22 by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------

#ifndef R3BCalifavsTofDOnlineSpectra_H
#define R3BCalifavsTofDOnlineSpectra_H 1

#include "R3BFrsData.h"
#include "FairTask.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

class TClonesArray;
class TH2F;
class R3BEventHeader;
class TVector3;

class R3BCalifavsTofDOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BCalifavsTofDOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BCalifavsTofDOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BCalifavsTofDOnlineSpectra();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    /**
     * Method for setting min proton energy (in keV) for opening angle histogram
     */
    inline void SetMinProtonEnergyForOpening(Float_t min) { fMinProtonE = min; }

    /**
     * Method to reset histograms
     */
    void Reset_Histo();

    /**
     * Method for setting the fTpat
     */
    void SetTpat(Int_t tpat) { fTpat = tpat; }

    /**
     * Method for setting the charge of the nuclear residue
     */
    void SetZCharge(Float_t z) { fZselection = z; }

  private:
    bool isFootDetect(TVector3 hit);

    TClonesArray* fWRItemsCalifa;
    TClonesArray* fWRItemsMaster;
    TClonesArray* fMapItemsCalifa;
    TClonesArray* fCalItemsCalifa;
    TClonesArray* fHitItemsCalifa;
    TClonesArray* fHitItemsTofd;
    TClonesArray* fHitItemsLos;
    TClonesArray* fHitItemsFrs;
    TClonesArray* fHitItemsFoot;
    TClonesArray* fCalItemsFoot;

    R3BEventHeader* header;
    Int_t fNEvents;
    Int_t fTpat;
    Float_t fZselection;
    Float_t fMinProtonE; /* Min proton energy (in keV) to calculate the opening angle */

    /*TH2F* fh2_Califa_coinPhi;
    TH2F* fh2_Califa_coinTheta;
    TCanvas* cCalifa_angles;

    TH2F* fh2_Califa_theta_phi[3]; // 0: all, 1: with TofD
    TH2F* fh2_Califa_NsNf[4];
    TH1F* fh2_Califa_total_energy[3];
    TH1F* fh2_Califa_CrystalHits[3];
    TH1F* fh2_Califa_ClusterId[3];
    TH2F* fh2_Califa_theta_energy[3];

    TH1F* fh2_openangle[3];
    TH2F* fh2_coinTheta[3];
    TH2F* fh2_coinPhi[3];
    TH2F* fh2_coinE[3];
    TH2F* fh2_leftE_openangle[3];
    TH2F* fh2_rightE_openangle[3];

    TH2F* fh2_Q_tof;*/

    //TCanvas* cCalifa_angles;

    TFile *outfile;
    TTree *outtree;

    // event variables
    Int_t tpatval;
    Int_t frshits;
    Int_t loshits;
    //Int_t trigger;

    // califa mapped variables
    std::vector<UShort_t> map_crystalId;
    std::vector<int16_t> map_energy;
    std::vector<int16_t> map_Ns;
    std::vector<int16_t> map_Nf;
    std::vector<int16_t> map_febexTime;
    std::vector<int16_t> map_wrts;
    std::vector<int16_t> map_overflow;
    std::vector<int16_t> map_pileup;
    std::vector<int16_t> map_discard;
    
    // califa cal variables
    std::vector<Int_t> cal_crystalId;
    std::vector<Double_t> cal_energy;
    std::vector<Double_t> cal_Ns;
    std::vector<Double_t> cal_Nf;
    std::vector<uint64_t> cal_time; 

    std::vector<Double_t> thetaList;
    std::vector<Double_t> phiList;
    std::vector<Double_t> energyList;
    std::vector<Double_t> NfList;
    std::vector<Double_t> NsList;
    std::vector<Double_t> crystalHitsList;
    std::vector<Double_t> clusterIdList;
    std::vector<ULong64_t> timeList;

    // tofd variables
    std::vector<Double_t> charge;
    std::vector<Double_t> tof;
    std::vector<Int_t> detectorId;

    // foot variables
    std::vector<Int_t> foot_detectorId;
    std::vector<Int_t> foot_nbHits;
    std::vector<Int_t> foot_mulStrip;
    std::vector<Double_t> foot_pos;
    std::vector<TVector3> foot_poslab;
    std::vector<Double_t> foot_theta;
    std::vector<Double_t> foot_phi;
    std::vector<Double_t> foot_energy; 

    // WR items
    std::vector<int64_t> wr;
    int64_t wrm;

  public:
    ClassDef(R3BCalifavsTofDOnlineSpectra, 1)
};

#endif /* R3BCalifavsTofDOnlineSpectra_H */
