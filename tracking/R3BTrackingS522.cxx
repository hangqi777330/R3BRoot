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

#include "R3BTrackingS522.h"
#include "R3BFootHitData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BEventHeader.h"
#include "R3BLosHitData.h"
#include "R3BTofdHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BFrsData.h"

#include "R3BMCTrack.h"
#include "R3BMDFWrapper.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TCutG.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include <TRandom3.h>
#include <TRandomGen.h>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
R3BTrackingS522* gMDFTrackerS522;

R3BTrackingS522::R3BTrackingS522()
    : R3BTrackingS522("TrackingS522", 1)
{
}

R3BTrackingS522::R3BTrackingS522(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fNEvents(0)
    , maxevent(0)
    , DoAlignment(false)
    , fTrackItems(new TClonesArray("R3BTrack"))
    , reference_PoQ(0.)
    , GladCurrent(-1)
    , GladReferenceCurrent(-1)
    , FootEnergyMin(-1)
    , FootEnergyMax(-1)
    , FiberEnergyMin(-1)
    , FiberEnergyMax(-1)
      , fHeader(nullptr)
{
}

R3BTrackingS522::~R3BTrackingS522()
{
    if (fTrackItems)
        delete fTrackItems;
}

InitStatus R3BTrackingS522::Init()
{
    LOG(info) << "R3BTrackingS522::Init()";
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(fatal) << "FairRootManager not found";
    }
    fHeader = (R3BEventHeader*)mgr->GetObject("EventHeader.");
    //fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!fHeader)
    {
        LOG(warn) << "R3BTrackingS522::Init() EventHeader. not found";
    }
    // Reading all detector branches
    cout << "\nDET_MAX = " << DET_MAX << endl;
    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    LOG(info) << "Reading " << NOF_FIB_DET << " fiber detectors";
    for (int det = 0; det < DET_MAX; det++)
    {
        fDataItems.push_back((TClonesArray*)mgr->GetObject(Form("%s", fDetectorNames[det])));
        if (NULL == fDataItems.at(det))
        {
            R3BLOG(fatal, Form("\n\n Cannot find tree branch %s \n\n", fDetectorNames[det]));
        }
    }
    // check if all cuts are properly set
    if (GladCurrent < 0 || GladReferenceCurrent < 0 || 
            FootEnergyMin < 0 || FootEnergyMax < 0 || 
            FiberEnergyMin < 0 || FiberEnergyMax < 0)
    {
        R3BLOG(fatal, Form(" Some cuts are not set or negative values are used\n\n"));
    }
    // Initializing all MDF functions
    LOG(info) << "Reading MDF function for FlightPath";
    MDF_FlightPath = new R3BMDFWrapper(MDF_FlightPath_filename.Data());

    LOG(info) << "Reading MDF function for PoQ";
    MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_filename.Data());

    LOG(info) << "Reading MDF function for TX0";
    MDF_TX0 = new R3BMDFWrapper(MDF_TX0_filename.Data());

    LOG(info) << "Reading MDF function for TY0";
    MDF_TY0 = new R3BMDFWrapper(MDF_TY0_filename.Data());

    LOG(info) << "Reading MDF function for TX1";
    MDF_TX1 = new R3BMDFWrapper(MDF_TX1_filename.Data());

    LOG(info) << "Reading MDF function for TY1";
    MDF_TY1 = new R3BMDFWrapper(MDF_TY1_filename.Data());

    // linking to global pointer (needed by alignment)
    gMDFTrackerS522 = this;

    tree_out.SetName("tree_out");
    tree_out.Branch("N_glob_tracks", &N_glob_tracks, "N_glob_tracks/i");
    tree_out.Branch("N_in_tracks",   &N_in_tracks,   "N_in_tracks/i");
    tree_out.Branch("N_out_tracks",  &N_out_tracks,  "N_out_tracks/i");

    tree_out.Branch("mul_m0",   &mul_m0,   "mul_m0/i");
    tree_out.Branch("mul_m1",   &mul_m1,   "mul_m1/i");
    tree_out.Branch("mul_foot", &mul_foot, "mul_foot/i");
    tree_out.Branch("mul_f1",   &mul_f1,   "mul_f1/i");
    tree_out.Branch("mul_f2",   &mul_f2,   "mul_f2/i");
    tree_out.Branch("mul_f15",  &mul_f15,  "mul_f15/i");
    tree_out.Branch("mul_f16",  &mul_f16,  "mul_f16/i");
    tree_out.Branch("mul_f30",  &mul_f30,  "mul_f30/i");
    tree_out.Branch("mul_f31",  &mul_f31,  "mul_f31/i");
    tree_out.Branch("mul_f32",  &mul_f32,  "mul_f32/i");
    tree_out.Branch("mul_f33",  &mul_f33,  "mul_f33/i");
    tree_out.Branch("mul_tofd", &mul_tofd, "mul_tofd/i");

    tree_out.Branch("Tpat",  &Tpat,  "Tpat/I");
    tree_out.Branch("PoQ",  PoQ,  "PoQ[N_glob_tracks]/D");
    tree_out.Branch("FlightPath",  FlightPath,  "FlightPath[N_glob_tracks]/D");
    tree_out.Branch("ToF",  ToF,  "ToF[N_glob_tracks]/D");
    tree_out.Branch("TX0",  TX0,  "TX0[N_glob_tracks]/D");
    tree_out.Branch("TX1",  TX1,  "TX1[N_glob_tracks]/D");
    tree_out.Branch("TY0",  TY0,  "TY0[N_glob_tracks]/D");
    tree_out.Branch("TY1",  TY1,  "TY1[N_glob_tracks]/D");
    tree_out.Branch("Beta",  Beta,  "Beta[N_glob_tracks]/D");
    tree_out.Branch("Gamma", Gamma, "Gamma[N_glob_tracks]/D");
    tree_out.Branch("mdf_AoZ", mdf_AoZ, "mdf_AoZ[N_glob_tracks]/D");

    tree_out.Branch("m0_X",    &m0_X,  "m0_X/F");
    tree_out.Branch("m0_Y",    &m0_Y,  "m0_Y/F");
    tree_out.Branch("m0_Z",    &m0_Z,  "m0_Z/F");

    tree_out.Branch("m1_X",    &m1_X,  "m1_X/F");
    tree_out.Branch("m1_Y",    &m1_Y,  "m1_Y/F");
    tree_out.Branch("m1_Z",    &m1_Z,  "m1_Z/F");

    tree_out.Branch("vertex_mwpc_X",    vertex_mwpc_X,  "vertex_mwpc_X[N_in_tracks]/F");
    tree_out.Branch("vertex_mwpc_Y",    vertex_mwpc_Y,  "vertex_mwpc_Y[N_in_tracks]/F");

    tree_out.Branch("vertex_foot_X",    vertex_foot_X,  "vertex_foot_X[N_in_tracks]/F");
    tree_out.Branch("vertex_foot_Y",    vertex_foot_Y,  "vertex_foot_Y[N_in_tracks]/F");

    tree_out.Branch("f1_X",    f1_X,  "f1_X[N_glob_tracks]/F");
    tree_out.Branch("f1_Y",    f1_Y,  "f1_Y[N_glob_tracks]/F");
    tree_out.Branch("f1_Z",    f1_Z,  "f1_Z[N_glob_tracks]/F");

    tree_out.Branch("f2_X",    f2_X,  "f2_X[N_glob_tracks]/F");
    tree_out.Branch("f2_Y",    f2_Y,  "f2_Y[N_glob_tracks]/F");
    tree_out.Branch("f2_Z",    f2_Z,  "f2_Z[N_glob_tracks]/F");

    tree_out.Branch("f15_X",    f15_X,  "f15_X[N_glob_tracks]/F");
    tree_out.Branch("f15_Y",    f15_Y,  "f15_Y[N_glob_tracks]/F");
    tree_out.Branch("f15_Z",    f15_Z,  "f15_Z[N_glob_tracks]/F");

    tree_out.Branch("f16_X",    f16_X,  "f16_X[N_glob_tracks]/F");
    tree_out.Branch("f16_Y",    f16_Y,  "f16_Y[N_glob_tracks]/F");
    tree_out.Branch("f16_Z",    f16_Z,  "f16_Z[N_glob_tracks]/F");

    tree_out.Branch("f30_X",    f30_X,  "f30_X[N_glob_tracks]/F");
    tree_out.Branch("f30_Y",    f30_Y,  "f30_Y[N_glob_tracks]/F");
    tree_out.Branch("f30_Z",    f30_Z,  "f30_Z[N_glob_tracks]/F");

    tree_out.Branch("f32_X",    f32_X,  "f32_X[N_glob_tracks]/F");
    tree_out.Branch("f32_Y",    f32_Y,  "f32_Y[N_glob_tracks]/F");
    tree_out.Branch("f32_Z",    f32_Z,  "f32_Z[N_glob_tracks]/F");

    tree_out.Branch("last_X",    last_X,  "last_X[N_glob_tracks]/F");
    tree_out.Branch("last_Y",    last_Y,  "last_Y[N_glob_tracks]/F");
    tree_out.Branch("last_Z",    last_Z,  "last_Z[N_glob_tracks]/F");

    tree_out.Branch("tofd_Q",    tofd_Q,  "tofd_Q[N_glob_tracks]/F");

    // Request storage of R3BTrack data in the output tree
    mgr->Register("MDFTracks", "MDFTracks data", fTrackItems, kTRUE);

    return kSUCCESS;
}

void R3BTrackingS522::Exec(Option_t* option)
{
    if (fNEvents / 1000. == (int)fNEvents / 1000)
        std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
            << " %) " << std::flush;
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");
    fNEvents += 1;
    N_glob_tracks=0;
    N_in_tracks=0;
    N_out_tracks=0;

    is_good_event = false;

    Tpat = fHeader->GetTpat();//vairable in the output tree
    //if(fTpat>0 && (((fHeader->GetTpat()) & fTpat) != fTpat)) return;
    if(fTpat<0 && fHeader->GetTpat()>64) return;//if Tpat is not set
    //if(fHeader->GetTpat()>64) return;//if Tpat is not set
    mul_m0   = fDataItems[MWPC0_HITDATA]->GetEntriesFast();
    mul_m1   = fDataItems[MWPC1_HITDATA]->GetEntriesFast();
    mul_foot = fDataItems[FOOT_HITDATA]->GetEntriesFast();
    mul_f32  = fDataItems[DET_FI32]->GetEntriesFast();
    mul_f30  = fDataItems[DET_FI30]->GetEntriesFast();
    mul_f31  = fDataItems[DET_FI31]->GetEntriesFast();
    mul_f33  = fDataItems[DET_FI33]->GetEntriesFast();
    mul_tofd = fDataItems[DET_TOFD]->GetEntriesFast();
    
    //cout << "\n\nMul m0: \t" << mul_m0;
    //cout << "\nMul m1: \t" << mul_m1;
    //cout << "\nMul foot: \t" << mul_foot;
    //cout << "\nMul f30: \t" << mul_f30;
    //cout << "\nMul f31: \t" << mul_f31;
    //cout << "\nMul f32: \t" << mul_f32;
    //cout << "\nMul f33: \t" << mul_f33;
    //cout << "\nMul tofd: \t" << mul_tofd;


    if(mul_m0!=1 || mul_m1!=1 || mul_foot<1) return;//for now take only mul=1 in mwpcs
    if(mul_f32<1 || mul_f30<1 || (mul_f31==0 && mul_f33==0) || mul_tofd<1) return;
    if(mul_f32>10 || mul_f30>10 || mul_f31>10 || mul_f33 >10 || mul_tofd>10) return;

    //cout << "\nGood event!\n";
       
    //FRS data
    auto frs_DataItems = fDataItems.at(FRS_DATA);
    if(frs_DataItems->GetEntriesFast() < 1) return; 
    auto frs_data = (R3BFrsData*)frs_DataItems->At(0);
    //if(frs_data->GetBrho()<16 || frs_data->GetBrho()>19) return;
cout << "\nGood event!\n";
        //------ Get TOFD data 
        R3BTofdHitData* tofd_hit{};
    bool is_good_tofd = false;
    for (auto i = 0; i < fDataItems[DET_TOFD]->GetEntriesFast(); ++i)
    {
        tofd_hit = static_cast<R3BTofdHitData*>(fDataItems[DET_TOFD]->At(i));
        if (tofd_hit->GetDetId() == 1) // only hits from first plane
        {
            is_good_tofd = true;
            break;
        }
    }
    if(!is_good_tofd){ 
	    cout<<"Bad_tofd"<<endl; 
	    return;
    }
    if(!MakeIncomingTracks()){ 
	    cout<<"Bad_incoming"<<endl; 
	    return;//at least one good track candidate in FOOT
    }
    if(!MakeOutgoingTracks()){ 
	    cout<<"Bad_outgoing"<<endl; 
	    return;//at least one good track candidate in Fibers
    }
    if(mul_f1>10 || mul_f2>10 || mul_f15>10 || mul_f16 >10) return;

    is_good_event = true;

    if (DoAlignment)
    {
        det_points align_data;
        align_data.f1.SetXYZ(f1_point.X(), 0, f1_point.Z());
        align_data.f2.SetXYZ(0, f2_point.Y(), f2_point.Z());
        align_data.f30.SetXYZ(f30_point.X(), 0, f30_point.Z());
        align_data.f32.SetXYZ(0, f32_point.Y(), f32_point.Z());
        align_data.flast.SetXYZ(flast_point.X(), 0, flast_point.Z());
	det_points_vec.push_back(align_data);
    }
    
    double delta_TY0, delta_TX0;
    for (auto & tin : tracks_in){
        for (auto & tout : tracks_out){
            //preserve the order, it is expected by the MDF function!
            mdf_data[0] = tin.f2_y;
            mdf_data[1] = tin.f2_z;
            mdf_data[2] = tin.f1_x;
            mdf_data[3] = tin.f1_z;
            mdf_data[4] = tout.f32_x;
            mdf_data[5] = tout.f32_z;
            mdf_data[6] = (tout.last_x - tout.f32_x)/(tout.last_z - tout.f32_z);
            mdf_data[7] = (tout.f30_y  - tin.f2_y)/(tout.f30_z - tin.f2_z);

            // Calculate all required MDF values
            PoQ[N_glob_tracks] = MDF_PoQ->MDF(mdf_data) * GladCurrent / GladReferenceCurrent;
            FlightPath[N_glob_tracks] = MDF_FlightPath->MDF(mdf_data);
            TX0[N_glob_tracks] = MDF_TX0->MDF(mdf_data);
            TX1[N_glob_tracks] = MDF_TX1->MDF(mdf_data);
            TY0[N_glob_tracks] = MDF_TY0->MDF(mdf_data);
            TY1[N_glob_tracks] = MDF_TY1->MDF(mdf_data);

            delta_TY0 = TY0[N_glob_tracks] - (tin.f16_y-tin.f2_y)/(tin.f16_z - tin.f2_z);
            if(delta_TY0<(-0.002) || delta_TY0>(0.004))
            //if(delta_TY0<(0.014) || delta_TY0>(0.02))//f31 run 159
                continue;

            delta_TX0 = TX0[N_glob_tracks] - (tin.f15_x - tin.f1_x)/(tin.f15_z-tin.f1_z);
            if(delta_TX0<(-0.029) || delta_TX0>(-0.023))//f31
            //if(delta_TX0<(-0.02) || delta_TX0>(-0.012))//f33
            //if(delta_TX0<(-0.022) || delta_TX0>(-0.012))//f31 run 159
                continue;

            f1_X[N_glob_tracks]   = tin.f1_x;
            f1_Z[N_glob_tracks]   = tin.f1_z;

            f2_Y[N_glob_tracks]   = tin.f2_y;
            f2_Z[N_glob_tracks]   = tin.f2_z;

            f15_X[N_glob_tracks]   = tin.f15_x;
            f15_Z[N_glob_tracks]   = tin.f15_z;

            f16_Y[N_glob_tracks]   = tin.f16_y;
            f16_Z[N_glob_tracks]   = tin.f16_z;

            f30_Y[N_glob_tracks]  = tout.f30_y;
            f30_Z[N_glob_tracks]  = tout.f30_z;

            f32_X[N_glob_tracks]  = tout.f32_x;
            f32_Z[N_glob_tracks]  = tout.f32_z;

            last_X[N_glob_tracks] = tout.last_x;
            last_Z[N_glob_tracks] = tout.last_z;

            tofd_Q[N_glob_tracks] = tofd_hit->GetEloss();

            N_glob_tracks++;

            if(N_glob_tracks == N_glob_tracks_max) return;
        }
    }
    return;
}

void R3BTrackingS522::FinishEvent()
{
    for (auto& DataItem : fDataItems)
    {
        DataItem->Clear();
    }
    fTrackItems->Clear();
    f1_hits.clear();
    f2_hits.clear();
    f15_hits.clear();
    f16_hits.clear();
    tracks_in.clear();
    tracks_out.clear();
    if(is_good_event) {
	cout<<"Good event***"<<endl;	    
	    tree_out.Fill();
    }

    //cout << "\n\nMul m0: " << mul_m0;
    //cout << "\nMul m1: " << mul_m1;
    //cout << "\nMul foot: " << mul_foot;
    //cout << "\nMul f30: " << mul_f30;
    //cout << "\nMul f31: " << mul_f31;
    //cout << "\nMul f32: " << mul_f32;
    //cout << "\nMul f33: " << mul_f33;
    //cout << "\nMul tofd: " << mul_tofd;
}

void R3BTrackingS522::FinishTask()
{
    LOG(info) << "Processed " << fNEvents << " events\n\n";
    if (DoAlignment) Alignment();
    tree_out.Write();
    //cout << "\n\n------- Statisitcs summary --------- ";
}

bool R3BTrackingS522::IsGoodFootHit(R3BFootHitData* fhit)
{
    if(fhit->GetEnergy()>FootEnergyMin && fhit->GetEnergy()<FootEnergyMax)
        return true;
    else 
        return false;
}

bool R3BTrackingS522::IsGoodFiberHit(R3BFiberMAPMTHitData* fhit)
{
    if((fhit->GetEloss() > FiberEnergyMin) && (fhit->GetEloss() < FiberEnergyMax) && 
            (fhit->GetTime() < 20000 && fhit->GetTime()>(-20000) )
      )
        return true;
    else 
        return false;
}

bool R3BTrackingS522::SortFootData()
{
    f1_hits.clear(); f2_hits.clear(); f15_hits.clear(); f16_hits.clear();
    mul_foot = fDataItems[FOOT_HITDATA]->GetEntriesFast();
    if(mul_foot==0) return false;
    for (auto f=0; f<mul_foot; ++f)
    {
        auto foot_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f));
        if(!IsGoodFootHit(foot_hit)) continue;
        switch(foot_hit->GetDetId())
        {
            case 1 ://Second FOOT1 for (Y)
                f1_hits.push_back(f);
                break;
            case 2 ://First FOOT2 (X)
                f2_hits.push_back(f);
                break;
            case 15 ://Last FOOT15 (X) 
                f15_hits.push_back(f);
                break;
            case 16 ://3rd FOOT16 (Y)
                f16_hits.push_back(f);
                break;
        }
    }
    if(f1_hits.empty() || f2_hits.empty() || f15_hits.empty() || f16_hits.empty())     
        return false; 
    mul_f1  = f1_hits.size(); mul_f2  = f2_hits.size();
    mul_f15 = f15_hits.size(); mul_f16 = f16_hits.size();
    return true;
}

bool R3BTrackingS522::MakeIncomingTracks()
{
    if(fDataItems[MWPC0_HITDATA]->GetEntriesFast()==0 || fDataItems[MWPC1_HITDATA]->GetEntriesFast()==0) 
        return false;
    tracks_in.clear();
    N_in_tracks=0;
    if(!SortFootData()) return false;//at least 1 hit in every FOOT
    TVector3 vertex_mwpc, vertex_foot;//projection to the center of the target (0,0,0)
    double tx_in, ty_in, dx_vertex, dy_vertex;
    Track tr;
    //Get MWPC hits, for now only first hit
    auto m0_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC0_HITDATA]->At(0));
    auto m1_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC1_HITDATA]->At(0));
    m0_point.SetXYZ(m0_hit->GetX()*0.1, m0_hit->GetY()*0.1, 0.);//cm
    m1_point.SetXYZ(m1_hit->GetX()*0.1, m1_hit->GetY()*0.1, 0.);//cm
    TransformPoint(m0_point, &m0_angles, &m0_position);//lab
    TransformPoint(m1_point, &m1_angles, &m1_position);//lab
    //Fill output tree variables
    m0_X = m0_point.X();  m0_Y = m0_point.Y();  m0_Z = m0_point.Z();
    m1_X = m1_point.X();  m1_Y = m1_point.Y();  m1_Z = m1_point.Z();
	cout<<"m0_x--> "<<m0_X<<endl;
    //------- Project mwpc track to the center of the target
    tx_in = (m0_point.X() - m1_point.X())/(m0_point.Z() - m1_point.Z());
    ty_in = (m0_point.Y() - m1_point.Y())/(m0_point.Z() - m1_point.Z()); 
    vertex_mwpc.SetX(m1_point.X() - tx_in * m1_point.Z());
    vertex_mwpc.SetY(m1_point.Y() - ty_in * m1_point.Z());
    vertex_mwpc.SetZ(0);

    //----- Make track candidates in FOOT and project them to the center of the target
    for (auto & f2 : f2_hits){
        auto f2_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f2));
        f2_point.SetXYZ(0, f2_hit->GetPosLab().Y() * 0.1, 0); //cm
        TransformPoint(f2_point,  &f2_angles,  &f2_position);

        for (auto & f1 : f1_hits){
            auto f1_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f1));
            f1_point.SetXYZ((-1.) * f1_hit->GetPosLab().X() * 0.1, 0, 0); //cm
            TransformPoint(f1_point,  &f1_angles,  &f1_position);

            for (auto & f16 : f16_hits){
                auto f16_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f16));
                f16_point.SetXYZ(0, f16_hit->GetPosLab().Y() * 0.1, 0); //cm
                TransformPoint(f16_point, &f16_angles, &f16_position);

                for (auto & f15 : f15_hits){
                    auto f15_hit = static_cast<R3BFootHitData*>(fDataItems[FOOT_HITDATA]->At(f15));
                    f15_point.SetXYZ((-1.) * f15_hit->GetPosLab().X() * 0.1, 0, 0); //cm
                    TransformPoint(f15_point, &f15_angles, &f15_position);

                    //project foot track to the center of the target
                    tr.tx0 = (f15_point.X() - f1_point.X())/(f15_point.Z() - f1_point.Z());
                    tr.ty0 = (f16_point.Y() - f2_point.Y())/(f16_point.Z() - f2_point.Z());
                    vertex_foot.SetX(f1_point.X() - tr.tx0 * f1_point.Z());
                    vertex_foot.SetY(f2_point.Y() - tr.ty0 * f2_point.Z());
                    vertex_foot.SetZ(0);
                    //Condition on the vertex matching
                    dx_vertex = vertex_foot.X() - vertex_mwpc.X();
                    dy_vertex = vertex_foot.Y() - vertex_mwpc.Y();
                    
                    //s522
                    if(dx_vertex<1 || dx_vertex>2 || dy_vertex<(-2.5) || dy_vertex>(-1.5)) continue;
                    //s509
                    //if(dx_vertex<1 || dx_vertex>2 || dy_vertex<(-1.9) || dy_vertex>(-1.)) continue;
                    
                    tr.f1_x   = f1_point.X();
                    tr.f1_z   = f1_point.Z();
                    tr.f2_y   = f2_point.Y();
                    tr.f2_z   = f2_point.Z();
                    tr.f15_x  = f15_point.X();
                    tr.f15_z  = f15_point.Z();
                    tr.f16_y  = f16_point.Y();
                    tr.f16_z  = f16_point.Z();
                    tracks_in.push_back(tr);

                    //Fill output tree
                    vertex_mwpc_X[N_in_tracks] = vertex_mwpc.X();
                    vertex_mwpc_Y[N_in_tracks] = vertex_mwpc.Y();
                    vertex_foot_X[N_in_tracks] = vertex_foot.X();
                    vertex_foot_Y[N_in_tracks] = vertex_foot.Y();
                    N_in_tracks++;
                    if(N_in_tracks==N_glob_tracks_max/2) return true;
                }
            }
        }
    }
    if(tracks_in.empty()) return false;
    else return true;
}

bool R3BTrackingS522::MakeOutgoingTracks()
{
    if(fDataItems[DET_FI32]->GetEntriesFast() == 0 || fDataItems[DET_FI30]->GetEntriesFast() == 0 || 
            (fDataItems[DET_FI33]->GetEntriesFast() == 0 && fDataItems[DET_FI31]->GetEntriesFast()==0) ) 
        return false;
    tracks_out.clear();
    Track tr;
    N_out_tracks=0;
    double angle_out, f30_slope, f30_offset, track_slope, track_offset;
    TVector3 f30_edge[2];//to extract z and x in f30
    for (auto i=0; i<fDataItems[DET_FI32]->GetEntriesFast(); ++i)
    {
        auto f32 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI32]->At(i));
        if(!IsGoodFiberHit(f32)) continue;
        f32_point.SetXYZ(f32->GetX(), 0, 0); //cm
        TransformPoint(f32_point, &f32_angles, &f32_position);
        tr.f32_x = f32_point.X();
        tr.f32_z = f32_point.Z();

        for (auto j=0; j<fDataItems[DET_FI30]->GetEntriesFast(); ++j)
        {
            auto f30 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI30]->At(j));
            if(!IsGoodFiberHit(f30)) continue;
            f30_point.SetXYZ(0,f30->GetY(), 0); //cm
            TransformPoint(f30_point, &f30_angles, &f30_position);
            tr.f30_y = f30_point.Y();

            if(fabs(f32->GetTime_ns() - f30->GetTime_ns())>30) continue;

            //make combination with every hit in fibers 33 and 31

            //for (auto k = 0; k<fDataItems[DET_FI33]->GetEntriesFast(); ++k)
            //{
            //    auto f33 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI33]->At(k));
            //    if(!IsGoodFiberHit(f33)) continue; //Messel side

            //    if(fabs(f32->GetTime_ns() - f33->GetTime_ns())>40) continue;

            //    last_point.SetXYZ(f33->GetX(), 0, 0); //cm
            //    TransformPoint(last_point, &f33_angles, &f33_position);
            //    tr.last_x = last_point.X();
            //    tr.last_z = last_point.Z();
            //    angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
            //    if(angle_out>(-10.) || angle_out<(-18.)) continue;
            //    // We need to extrapolate Z position in f30 because it was used for Y measurement
            //    // Define two (X,Z) points on the f30 plane:
            //    //Now track every combination of upstream and downstream tracks 
            //    f30_edge[0].SetXYZ(-1, 0, 0);
            //    f30_edge[1].SetXYZ(1, 0, 0);
            //    TransformPoint(f30_edge[0], &f30_angles, &f30_position);
            //    TransformPoint(f30_edge[1], &f30_angles, &f30_position);
            //    // Parameterize f30 plane
            //    f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
            //    f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
            //    track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
            //    track_offset = (tr.last_x - track_slope * tr.last_z);
            //    // Extrapolate final X and Z position in f30
            //    tr.f30_z = (track_offset - f30_offset) / (f30_slope - track_slope);// extrapolated
            //    tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
            //    N_out_tracks++;
            //    tracks_out.push_back(tr);
            //    if(N_out_tracks==N_glob_tracks/2) return true;
            //}


            //make combination with every hit in fibers 33 and 31
            for (auto k = 0; k<fDataItems[DET_FI31]->GetEntriesFast(); ++k)
            {
                auto f31 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI31]->At(k));
                if(!IsGoodFiberHit(f31)) continue; //Messel side

                if((f32->GetTime_ns() - f31->GetTime_ns())>20 ||
                        (f32->GetTime_ns() - f31->GetTime_ns())<(-10) ) continue;

                flast_point.SetXYZ(f31->GetX(), 0, 0); //cm
                TransformPoint(flast_point, &f31_angles, &f31_position);
                tr.last_x = flast_point.X();
                tr.last_z = flast_point.Z();
                angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
                if(angle_out>(-10.) || angle_out<(-18.)) continue;
                // We need to extrapolate Z position in f30 because it was used for Y measurement
                // Define two (X,Z) points on the f30 plane:
                //Now track every combination of upstream and downstream tracks 
                f30_edge[0].SetXYZ(-1, 0, 0);
                f30_edge[1].SetXYZ(1, 0, 0);
                TransformPoint(f30_edge[0], &f30_angles, &f30_position);
                TransformPoint(f30_edge[1], &f30_angles, &f30_position);
                // Parameterize f30 plane
                f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
                f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
                track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
                track_offset = (tr.last_x - track_slope * tr.last_z);
                // Extrapolate final X and Z position in f30
                tr.f30_z = (track_offset - f30_offset) / (f30_slope - track_slope);// extrapolated
                tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
                N_out_tracks++;
                tracks_out.push_back(tr);
                if(N_out_tracks==N_glob_tracks/2) return true;
            }
        }
    }
    if(tracks_out.empty()) return false;
    return true;
}

void R3BTrackingS522::TransformPoint(TVector3& point, TVector3* rot, TVector3* trans)
{
    r.SetToIdentity();
    // First Euler rotation around Y axis
    r.RotateY(rot->Y());
    // get local X axis after first rotation
    v3_localX.SetMagThetaPhi(1, r.ThetaX(), r.PhiX());
    // Second Euler rotation around local X axis
    r.Rotate(rot->X(), v3_localX);
    // get local Z axis after second rotation
    v3_localZ.SetMagThetaPhi(1, r.ThetaZ(), r.PhiZ());
    // final rotation around local Z axis
    r.Rotate(rot->Z(), v3_localZ);
    point.Transform(r);
    point += (*trans);
    return;
}

void R3BTrackingS522::TransformPoint1(TVector3& point, TVector3 rot, TVector3 trans)
{
    r.SetToIdentity();
    // First Euler rotation around Y axis
    r.RotateY(rot.Y());
    // get local X axis after first rotation
    v3_localX.SetMagThetaPhi(1, r.ThetaX(), r.PhiX());
    // Second Euler rotation around local X axis
    r.Rotate(rot.X(), v3_localX);
    // get local Z axis after second rotation
    v3_localZ.SetMagThetaPhi(1, r.ThetaZ(), r.PhiZ());
    // final rotation around local Z axis
    r.Rotate(rot.Z(), v3_localZ);
    point.Transform(r);
    point += (trans);
    return;
}

R3BTrack* R3BTrackingS522::AddTrackData(TVector3 mw, TVector3 poq, Double_t charge, Double_t aoz)
{
    // Filling output track info
    TClonesArray& clref = *fTrackItems;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTrack(mw.X(), mw.Y(), mw.Z(), poq.X(), poq.Y(), poq.Z(), charge, aoz, 0., 0., 0);
}

// Setup energy cuts in foot and fibers 
void R3BTrackingS522::SetFootEnergyMinMax(double min, double max)
{
    FootEnergyMin = min;
    FootEnergyMax = max;
    return;
}
void R3BTrackingS522::SetFiberEnergyMinMax(double min, double max)
{
    FiberEnergyMin = min;
    FiberEnergyMax = max;
    return;
}
void R3BTrackingS522::Alignment()
{
	int run_flag = 0;
	std::cout << "\n\n----------------------------------------------------";
/*	std::cout << "\n Ready for detector alignment using the following data set:";
	std::cout << "\n\tNumber of reference tracks: " << det_points_vec.size();
	std::cout << "\n\tFRS Brho: min = " << FrsBrhoMin << ", max = " << FrsBrhoMax;
	std::cout << "\n\tFRS PoQ: min = " << FrsBrhoMin / 3.3356 << ", max = " << FrsBrhoMax / 3.3356;
	std::cout << "\n\tP/Z reference: " << reference_PoQ << "GeV/c";
*/	std::cout << "\n\n Do you want to continue? (0 = no, 1 = yes):  ";
	std::cin >> run_flag;
	if (run_flag == 0)
		return;

	// Now define minimizer
	ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
	const UInt_t NVarsFunctor = 25; // number of the alignment offsets: (2 angles + 3 offsets) x 4 detectors
	const double* xs;
	double step[NVarsFunctor];
	double offset[NVarsFunctor];
	Double_t min_offset[NVarsFunctor];
	Double_t max_offset[NVarsFunctor];
	TH1F* hist[NVarsFunctor];
	char hname[100];

	// For every detector: 3 Euler angles and 3 shifts
	for (auto d = 0; d < 4; ++d)
	{
		for (auto a = 0; a < 2; ++a) // angle shifts in rad
		{
			min_offset[d * 5 + a] = -0.04;
			max_offset[d * 5 + a] = 0.04;
			step[d * 5 + a] = 0.001;
		}
		for (auto o = 0; o < 3; ++o) // position shifts in cm
		{
			min_offset[d * 5 + 2 + o] = -4.;
			max_offset[d * 5 + 2 + o] = 4.;
			step[d * 5 + 2 + o] = 0.01;
		}
	}
	for (UInt_t i = 0; i < NVarsFunctor; i++)
	{
		sprintf(hname, "par%d", i);
		hist[i] = new TH1F(hname, hname, 100, min_offset[i], max_offset[i]);
	}
	std::cout << "\n\n-- Perfroming minimization for the detector alignment. Please wait... \n";

	// Setting up Minimizer function parameters
	Double_t precision = 1e-10; // 0 - default precision will be automaticalle determined
	Double_t tolerance = 0.02;
	minimizer->SetMaxFunctionCalls(1000000000); // for Minuit/Minuit2
	minimizer->SetMaxIterations(100);           // for GSL
	minimizer->SetTolerance(tolerance);
	minimizer->SetPrecision(precision);
	minimizer->SetPrintLevel(2);
	minimizer->SetStrategy(0); // to run faster
	ROOT::Math::Functor f(&R3BTrackingS522::AlignmentErrorS522, NVarsFunctor);
	minimizer->SetFunction(f);
	minimizer->SetPrintLevel(0);

	// Getting experimental data and doing minimization
	UInt_t i;
	Int_t bs = 0;
	bool good_minimum = false;
	while (bs < 100) // running several minimizations
	{
		minimizer->Clear();
		for (i = 0; i < NVarsFunctor; i++) // sampling +=50% from limits
		{
			offset[i] = gRandom->Uniform(min_offset[i] * 0.5, max_offset[i] * 0.5);
			minimizer->SetVariable(i, Form("par%d", i), offset[i], step[i]);
			minimizer->SetVariableLimits(i, min_offset[i], max_offset[i]);
		}
		minimizer->Minimize();
		xs = minimizer->X();
		// if(minimizer->Status() !=0) continue; //valid minimum
		// Check if all paramters are "far" from limits
		for (i = 0; i < NVarsFunctor; i++)
		{
			if (fabs((xs[i] - min_offset[i]) / min_offset[i]) < 0.1 ||
					fabs((max_offset[i] - xs[i]) / max_offset[i]) < 0.1)
				break;
		}
		if (i != NVarsFunctor)
		{
			// std::cout << "\n\n -- Parameter is too close to the limit!! ";
			continue;
		}
		else
		{
			std::cout << "\n\n-- Good parameters are found!!";
			good_minimum = true;
		}
		std::cout << "\n-- Minimization No. " << bs << "\n-- Minimum: f(";
		for (i = 0; i < NVarsFunctor; i++)
		{
			std::cout << xs[i];
			hist[i]->Fill(xs[i]);
			if (i != NVarsFunctor - 1)
				std::cout << ",";
		}
		std::cout << ") = " << minimizer->MinValue();
		std::cout << "\n-- Minimizer status: " << minimizer->Status() << std::endl;
		bs++;
	}
	// Outputing final info and histograms
	if (minimizer->MinValue() < tolerance && f(xs) < tolerance)
		std::cout << "-- Minimizer "
			<< "   converged to the right minimum" << std::endl;
	else
	{
		std::cout << "-- Minimizer "
			<< "   failed to converge !!!" << std::endl;
	}
	TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
	c1->Divide(5, 4);
	for (UInt_t v = 0; v < NVarsFunctor; v++)
	{
		c1->cd(v + 1);
		hist[v]->Draw();
	}
	return;
}


double R3BTrackingS522::AlignmentErrorS522(const double* par)
{
	gMDFTrackerS522->f1_ang_offset.SetXYZ(par[0], par[1], par[2]);
	gMDFTrackerS522->f1_pos_offset.SetXYZ(par[3], 0 , par[4]);

	gMDFTrackerS522->f2_ang_offset.SetXYZ(par[5], par[6], par[7]);
	gMDFTrackerS522->f2_pos_offset.SetXYZ(0, par[8], par[9]);

	gMDFTrackerS522->f30_ang_offset.SetXYZ(par[10], par[11], par[12]);
	gMDFTrackerS522->f30_pos_offset.SetXYZ(par[13], 0, par[14]);

	gMDFTrackerS522->f32_ang_offset.SetXYZ(par[15], par[16], par[17]);
	gMDFTrackerS522->f32_pos_offset.SetXYZ(0, par[18], par[19]);

	gMDFTrackerS522->flast_ang_offset.SetXYZ(par[20], par[21], par[22]);
	gMDFTrackerS522->flast_pos_offset.SetXYZ(par[23], 0, par[24]);


	Double_t mdf_input[8]; // data container for the MDF function

	double v2 = 0;
	double v = 0;
	int counter = 0;
	for (auto& d : (gMDFTrackerS522->det_points_vec))
	{
		gMDFTrackerS522->f1_point = d.f1;
		gMDFTrackerS522->f2_point = d.f2;
		gMDFTrackerS522->f30_point = d.f30;
		gMDFTrackerS522->f32_point = d.f32;
		gMDFTrackerS522->flast_point = d.flast;	

		// This will transform "det_point" vectors into lab frame
		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f1_point,
				gMDFTrackerS522->GetEulerAnglesFoot1() + gMDFTrackerS522->f1_ang_offset,
				gMDFTrackerS522->GetPositionFoot1() + gMDFTrackerS522->f1_pos_offset);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f2_point,
				gMDFTrackerS522->GetEulerAnglesFoot2() + gMDFTrackerS522->f2_ang_offset,
				gMDFTrackerS522->GetPositionFoot2() + gMDFTrackerS522->f2_pos_offset);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f30_point,
				gMDFTrackerS522->GetEulerAnglesFiber30() + gMDFTrackerS522->f30_ang_offset,
				gMDFTrackerS522->GetPositionFiber30() + gMDFTrackerS522->f30_pos_offset);

			gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->f32_point,
				gMDFTrackerS522->GetEulerAnglesFiber32() + gMDFTrackerS522->f32_ang_offset,
				gMDFTrackerS522->GetPositionFiber32() + gMDFTrackerS522->f32_pos_offset);

		gMDFTrackerS522->TransformPoint1(gMDFTrackerS522->flast_point,
				gMDFTrackerS522->GetEulerAnglesFiber31() + gMDFTrackerS522->flast_ang_offset,
				gMDFTrackerS522->GetPositionFiber31() + gMDFTrackerS522->flast_pos_offset);



		mdf_input[0] = gMDFTrackerS522->f2_point.Y();
		mdf_input[1] = gMDFTrackerS522->f2_point.Z();
		mdf_input[2] = gMDFTrackerS522->f1_point.X();
		mdf_input[3] = gMDFTrackerS522->f1_point.Z();
		mdf_input[4] = gMDFTrackerS522->f32_point.X();
		mdf_input[5] = gMDFTrackerS522->f32_point.Z();	
		mdf_input[6] = (gMDFTrackerS522->flast_point.X() - mdf_input[4]) / (gMDFTrackerS522->flast_point.Z() - mdf_input[5]);
	mdf_input[7] = (gMDFTrackerS522->f30_point.Y() - mdf_input[0]) / (gMDFTrackerS522->f30_point.Z() - mdf_input[1]);


		v2 += pow((gMDFTrackerS522->Get_MDF_PoQ()->MDF(mdf_input) - gMDFTrackerS522->GetReferencePoQ()), 2);
		counter++;
	}
	v2 /= counter;
	v = sqrt(v2);
	// std::cout << "\nReturning error: " << v;
	return v;

}

ClassImp(R3BTrackingS522);
