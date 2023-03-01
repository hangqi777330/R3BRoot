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

// Created on 23/22/2022 by V.Panin

#ifndef R3BTRACKINGS522
#define R3BTRACKINGS522

#include "FairTask.h"
#include "R3BEventHeader.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFootHitData.h"
#include "R3BTofdHitData.h"
#include "R3BMDFWrapper.h"
#include "R3BTrack.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TRotation.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include "TTree.h"
class R3BTrackingS522 : public FairTask
{
  public:
    R3BTrackingS522();
    R3BTrackingS522(const char* name, Int_t iVerbose = 1);
    virtual ~R3BTrackingS522();
    virtual InitStatus Init();
    virtual void Exec(Option_t* option);
    virtual void FinishEvent();
    virtual void FinishTask();

    // GLAD current which was used for training MDF function
    // Should  be set form the steering macro
    void SetGladReferenceCurrent(Double_t cur) { GladReferenceCurrent = cur; }

    // GLAD current in the run being anlaysed, set from steering macro
    void SetGladCurrent(Double_t cur) { GladCurrent = cur; }

    // Set MDF functions from the steering macro
    void Set_MDF_PoQ(TString name) { MDF_PoQ_filename = name; }
    void Set_MDF_FlightPath(TString name) { MDF_FlightPath_filename = name; }
    void Set_MDF_TX0(TString name) { MDF_TX0_filename = name; }
    void Set_MDF_TY0(TString name) { MDF_TY0_filename = name; }
    void Set_MDF_TX1(TString name) { MDF_TX1_filename = name; }
    void Set_MDF_TY1(TString name) { MDF_TY1_filename = name; }
    void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    void SetTpat(Int_t tpat) { fTpat = tpat; }
    void SetMaxEvent(Int_t nev) { maxevent = nev; }

    // Set lab positions and angles of the detectors from the steering macro
    // LOS
    void SetPositionLos(double x, double y, double z) {  los_position.SetXYZ(x, y, z); }// cm
    void SetEulerAnglesLos(double x, double y, double z) { los_angles.SetXYZ(x, y, z); }// rad

    // MWPCs
    void SetPositionMwpc0(double x, double y, double z) {  m0_position.SetXYZ(x, y, z); }// cm
    void SetPositionMwpc1(double x, double y, double z) {  m1_position.SetXYZ(x, y, z); }// cm
    void SetEulerAnglesMwpc0(double x, double y, double z) {  m0_angles.SetXYZ(x, y, z); }// rad
    void SetEulerAnglesMwpc1(double x, double y, double z) {  m1_angles.SetXYZ(x, y, z); }// rad

    // FOOTs
    void SetPositionFoot1(double x, double y, double z) {  f1_position.SetXYZ(x, y, z); }// cm
    void SetPositionFoot2(double x, double y, double z) {  f2_position.SetXYZ(x, y, z); }// cm
    void SetPositionFoot15(double x, double y, double z) { f15_position.SetXYZ(x, y, z); }// cm
    void SetPositionFoot16(double x, double y, double z) { f16_position.SetXYZ(x, y, z); }// cm
    void SetEulerAnglesFoot1(double x, double y, double z) {  f1_angles.SetXYZ(x, y, z); }// rad
    void SetEulerAnglesFoot2(double x, double y, double z) {  f2_angles.SetXYZ(x, y, z); }// rad
    void SetEulerAnglesFoot15(double x, double y, double z) { f15_angles.SetXYZ(x, y, z); }// rad
    void SetEulerAnglesFoot16(double x, double y, double z) { f16_angles.SetXYZ(x, y, z); }// rad

    //Fibers
    void SetPositionFiber30(double x, double y, double z) { f30_position.SetXYZ(x, y, z); }// cm
    void SetPositionFiber31(double x, double y, double z) { f31_position.SetXYZ(x, y, z); }// cm
    void SetPositionFiber32(double x, double y, double z) { f32_position.SetXYZ(x, y, z); }// cm
    void SetPositionFiber33(double x, double y, double z) { f33_position.SetXYZ(x, y, z); }// cm
    void SetEulerAnglesFiber30(double x, double y, double z) { f30_angles.SetXYZ(x, y, z); }// rad
    void SetEulerAnglesFiber31(double x, double y, double z) { f31_angles.SetXYZ(x, y, z); }// rad
    void SetEulerAnglesFiber32(double x, double y, double z) { f32_angles.SetXYZ(x, y, z); }// rad
    void SetEulerAnglesFiber33(double x, double y, double z) { f33_angles.SetXYZ(x, y, z); }// rad

    //Tofd
    void SetPositionTofd(double x, double y, double z) { tofd_position.SetXYZ(x, y, z); }// cm
    void SetEulerAnglesTofd(double x, double y, double z) { tofd_angles.SetXYZ(x, y, z); }// rad

    // Get lab positions and angles (needed by alignment function)
    inline TVector3 GetPositionLos() { return los_position; }   // cm
    inline TVector3 GetPositionMwpc0() { return m0_position; }   // cm
    inline TVector3 GetPositionMwpc1() { return m1_position; }   // cm
    inline TVector3 GetPositionFoot1() { return f1_position; }   // cm
    inline TVector3 GetPositionFoot2() { return f2_position; }   // cm
    inline TVector3 GetPositionFoot15() { return f15_position; }   // cm
    inline TVector3 GetPositionFoot16() { return f16_position; }   // cm
    inline TVector3 GetPositionFiber30() { return f30_position; }   // cm
    inline TVector3 GetPositionFiber31() { return f31_position; }   // cm
    inline TVector3 GetPositionFiber32() { return f32_position; }   // cm
    inline TVector3 GetPositionFiber33() { return f33_position; }   // cm
    inline TVector3 GetPositionTofd() { return tofd_position; } // cm

    inline TVector3 GetEulerAnglesLos() { return los_angles; }   // cm
    inline TVector3 GetEulerAnglesMwpc0() { return m0_angles; }   // cm
    inline TVector3 GetEulerAnglesMwpc1() { return m1_angles; }   // cm
    inline TVector3 GetEulerAnglesFoot1() { return f1_angles; }   // cm
    inline TVector3 GetEulerAnglesFoot2() { return f2_angles; }   // cm
    inline TVector3 GetEulerAnglesFoot15() { return f15_angles; }   // cm
    inline TVector3 GetEulerAnglesFoot16() { return f16_angles; }   // cm
    inline TVector3 GetEulerAnglesFiber30() { return f30_angles; }   // cm
    inline TVector3 GetEulerAnglesFiber31() { return f31_angles; }   // cm
    inline TVector3 GetEulerAnglesFiber32() { return f32_angles; }   // cm
    inline TVector3 GetEulerAnglesFiber33() { return f33_angles; }   // cm
    inline TVector3 GetEulerAnglesTofd() { return tofd_angles; } // rad

    R3BMDFWrapper* Get_MDF_PoQ() { return MDF_PoQ; }

    void SetTofOffset(double offset) { tof_offset = offset; } // ns

    // Transofrming input detector hit (point) into laboratory system
    void TransformPoint(TVector3& point, TVector3* rotation, TVector3* translation);
    void TransformPoint1(TVector3& point, TVector3 rotation, TVector3 translation);
    // Setup energy cuts in foot and fibers 
    void SetFootEnergyMinMax(double min, double max);
    void SetFiberEnergyMinMax(double min, double max);

    // Setters for the alignment procedure
    void SetDoAlignment(bool flag) { DoAlignment = flag; }
    void SetReferencePoQ(Double_t val) { reference_PoQ = val; }
    Double_t GetReferencePoQ() { return reference_PoQ; }

    void Alignment();
    static double AlignmentErrorS522(const double* par);

    //Storing indices of hits in TCA for potential track candidates
    struct Track
    {
	    double f1_x;
	    double f1_z;
	    double f2_y;
	    double f2_z;
	    double f15_x;
	    double f15_z;
	    double f16_y;
	    double f16_z;
	    double tx0;
	    double ty0;
	    double f32_x;
	    double f32_z;
	    double f30_x;
	    double f30_y;
	    double f30_z;
	    double last_x;
	    double last_z;
    };
    std::vector<Track> tracks_in; //track candidates in FOOT
    std::vector<Track> tracks_out;//track candidates in FOOT

    //Storage of hit indices
    vector<UInt_t> f1_hits;
    vector<UInt_t> f2_hits;
    vector<UInt_t> f15_hits;
    vector<UInt_t> f16_hits;
    TVector3 m0_point, m1_point, f1_point, f2_point, f15_point,f16_point, f30_point, f32_point, flast_point;

    bool IsGoodFiberHit(R3BFiberMAPMTHitData* fhit);
    bool IsGoodFootHit(R3BFootHitData* fhit);
    bool SortFootData();
    bool MakeIncomingTracks();
    bool MakeOutgoingTracks();
    // Data containesr needed only for the Alignment() function
    struct det_points
    {
	    TVector3 f1;
	    TVector3 f2;
	    TVector3 f30;
	    TVector3 f32;
	    TVector3 flast;
    };
    std::vector<det_points> det_points_vec;

    // TVector3 mwpc_point, f10_point, f11_point, f12_point, tofd_point;
    TVector3 f1_ang_offset, f2_ang_offset, f30_ang_offset, f32_ang_offset, flast_ang_offset;
    TVector3 f1_pos_offset, f2_pos_offset, f30_pos_offset, f32_pos_offset, flast_pos_offset;


  private:
    // Input hit data from the TClonesArray
    // do not change the order, add new det in the end
    enum DetectorInstances
    {
	    DET_FI_FIRST,
	    DET_FI30 = DET_FI_FIRST,
	    DET_FI31,
	    DET_FI32,
	    DET_FI33,
	    DET_FI_LAST = DET_FI33,
	    DET_TOFD,
	    FOOT_HITDATA,
	    MWPC0_HITDATA,
	    MWPC1_HITDATA,
	    FRS_DATA,
	    DET_MAX
    };

#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

    // Names of essential branches in the input tree
    // do not change the order! add new data in the end
    const char* fDetectorNames[DET_MAX + 1] = { "Fi30Hit", "Fi31Hit", "Fi32Hit", "Fi33Hit",
	    "TofdHit", "FootHitData", "Mwpc0HitData", "Mwpc1HitData", "FrsData", NULL };

    R3BEventHeader* fHeader;
    std::vector<TClonesArray*> fDataItems; // input data
    TClonesArray* fTrackItems;             // output data

    bool is_good_event;

    TRotation r;
    TVector3 v3_localX;
    TVector3 v3_localZ;

    TVector3 los_position;
    TVector3 los_angles;

    TVector3 m0_position;
    TVector3 m1_position;
    TVector3 f1_position;
    TVector3 f2_position;
    TVector3 f15_position;
    TVector3 f16_position;
    TVector3 f30_position;
    TVector3 f31_position;
    TVector3 f32_position;
    TVector3 f33_position;

    TVector3 m0_angles;
    TVector3 m1_angles;
    TVector3 f1_angles;
    TVector3 f2_angles;
    TVector3 f15_angles;
    TVector3 f16_angles;
    TVector3 f30_angles;
    TVector3 f31_angles;
    TVector3 f32_angles;
    TVector3 f33_angles;

    TVector3 tofd_position;
    TVector3 tofd_angles;

    R3BMDFWrapper* MDF_FlightPath;
    R3BMDFWrapper* MDF_PoQ;
    R3BMDFWrapper* MDF_TX0;
    R3BMDFWrapper* MDF_TX1;
    R3BMDFWrapper* MDF_TY0;
    R3BMDFWrapper* MDF_TY1;

    TString MDF_FlightPath_filename;
    TString MDF_PoQ_filename;
    TString MDF_TX0_filename;
    TString MDF_TY0_filename;
    TString MDF_TX1_filename;
    TString MDF_TY1_filename;

    Double_t mdf_data[8];   // data container for the MDF function
    unsigned long fNEvents; // Event counter
    Int_t fTrigger;
    Int_t fTpat;
    Int_t maxevent;
    Double_t GladCurrent;
    Double_t GladReferenceCurrent;
    Double_t reference_PoQ;
    Bool_t DoAlignment;
    Double_t tof_offset; // ns

    // Energy range in FOOT 
    Double_t FootEnergyMin;
    Double_t FootEnergyMax;

    // Cut on fiber hit energy set by SetFiberEnergyMinMax():
    Double_t FiberEnergyMin;
    Double_t FiberEnergyMax;
    
    TTree tree_out;
    UInt_t N_glob_tracks;
    UInt_t N_in_tracks;
    UInt_t N_out_tracks;
    static constexpr UInt_t N_glob_tracks_max = 100000;
    UInt_t mul_m0;
    UInt_t mul_m1;
    UInt_t mul_f1;
    UInt_t mul_f2;
    UInt_t mul_f15;
    UInt_t mul_f16;
    UInt_t mul_f30;
    UInt_t mul_f31;
    UInt_t mul_f32;
    UInt_t mul_f33;
    UInt_t mul_tofd;
    UInt_t mul_foot;

    Int_t Tpat;
    Double_t PoQ[N_glob_tracks_max];
    Double_t FlightPath[N_glob_tracks_max];
    Double_t ToF[N_glob_tracks_max];
    Double_t TX0[N_glob_tracks_max];
    Double_t TX1[N_glob_tracks_max];
    Double_t TY0[N_glob_tracks_max];
    Double_t TY1[N_glob_tracks_max];
    Double_t Beta[N_glob_tracks_max];
    Double_t Gamma[N_glob_tracks_max];
    Double_t mdf_AoZ[N_glob_tracks_max];

    Float_t  m0_X;
    Float_t  m0_Y;
    Float_t  m0_Z;
    
    Float_t  m1_X;
    Float_t  m1_Y;
    Float_t  m1_Z;

    Float_t vertex_foot_X[N_glob_tracks_max];
    Float_t vertex_foot_Y[N_glob_tracks_max];

    Float_t vertex_mwpc_X[N_glob_tracks_max];
    Float_t vertex_mwpc_Y[N_glob_tracks_max];
    
    Float_t f1_X[N_glob_tracks_max];
    Float_t f1_Y[N_glob_tracks_max];
    Float_t f1_Z[N_glob_tracks_max];
    Float_t f1_Q[N_glob_tracks_max];
    Float_t f1_T[N_glob_tracks_max];

    Float_t f2_X[N_glob_tracks_max];
    Float_t f2_Y[N_glob_tracks_max];
    Float_t f2_Z[N_glob_tracks_max];
    Float_t f2_Q[N_glob_tracks_max];
    Float_t f2_T[N_glob_tracks_max];

    Float_t f15_X[N_glob_tracks_max];
    Float_t f15_Y[N_glob_tracks_max];
    Float_t f15_Z[N_glob_tracks_max];
    Float_t f15_Q[N_glob_tracks_max];
    Float_t f15_T[N_glob_tracks_max];

    Float_t f16_X[N_glob_tracks_max];
    Float_t f16_Y[N_glob_tracks_max];
    Float_t f16_Z[N_glob_tracks_max];
    Float_t f16_Q[N_glob_tracks_max];
    Float_t f16_T[N_glob_tracks_max];

    Float_t f30_X[N_glob_tracks_max];
    Float_t f30_Y[N_glob_tracks_max];
    Float_t f30_Z[N_glob_tracks_max];
    Float_t f30_Q[N_glob_tracks_max];
    Float_t f30_T[N_glob_tracks_max];

    Float_t f32_X[N_glob_tracks_max];
    Float_t f32_Y[N_glob_tracks_max];
    Float_t f32_Z[N_glob_tracks_max];
    Float_t f32_Q[N_glob_tracks_max];
    Float_t f32_T[N_glob_tracks_max];

    Float_t last_X[N_glob_tracks_max];
    Float_t last_Y[N_glob_tracks_max];
    Float_t last_Z[N_glob_tracks_max];
    Float_t last_Q[N_glob_tracks_max];
    Float_t last_T[N_glob_tracks_max];

    Float_t tofd_X[N_glob_tracks_max];
    Float_t tofd_Y[N_glob_tracks_max];
    Float_t tofd_Z[N_glob_tracks_max];
    Float_t tofd_Q[N_glob_tracks_max];
    Float_t tofd_T[N_glob_tracks_max];
    Float_t tofd_ID[N_glob_tracks_max];


    // Essential constants
    const Double_t SPEED_OF_LIGHT = 29.9792458; // cm/ns
    const Double_t AMU = 0.9314940038;          // GeV/c2

    // Private method to fill output track data
    R3BTrack* AddTrackData(TVector3 mw, TVector3 poq, Double_t charge, Double_t aoz);

  public:
    ClassDef(R3BTrackingS522, 1)
};

#endif
