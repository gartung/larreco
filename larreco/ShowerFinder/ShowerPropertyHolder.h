//###################################################################
//### Name:        ShowerPropertyHolder                           ### 
//### Author:      Dominic Barker                                 ###
//### Date:        14.05.19                                       ###
//### Description: Class to holder the standard shower property   ###
//###              information. Used in SBNShower_modle and       ###
//###              corresponding tools.                           ###
//###################################################################

#ifndef ShowerPropertyHolder_H
#define ShowerPropertyHolder_H

//LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"

//Root Includes
#include "TVector3.h"

//C++ Inlcudes 
#include <vector>
#include <iostream> 

namespace reco {
  namespace shower {
    class ShowerPropertyHolder;
  }
}

class reco::shower::ShowerPropertyHolder {
 public: 

  //Constructor 
  ShowerPropertyHolder(){
      ShowerDirectionPtr     = 0;
      ShowerStartPositionPtr = 0;
      InitialTrackPtr        = 0;
      InitialTrackHitsPtr    = 0;
      ShowerEnergyPtr        = 0;
      ShowerdEdxPtr          = 0;
      BestPlanePtr           = 0;
      
      ShowerDirectionErrPtr     = 0;
      ShowerStartPositionErrPtr = 0;
      ShowerEnergyErrPtr        = 0;
      ShowerdEdxErrPtr          = 0;
  };

  //Reset the holder
  void clear(){
    ShowerDirectionPtr     = 0;
    ShowerStartPositionPtr = 0;
    InitialTrackPtr        = 0;
    InitialTrackHitsPtr    = 0;
    ShowerEnergyPtr        = 0;
    ShowerdEdxPtr          = 0;
    BestPlanePtr           = 0;

    ShowerDirectionErrPtr     = 0;
    ShowerStartPositionErrPtr = 0;
    ShowerEnergyErrPtr        = 0;
    ShowerdEdxErrPtr          = 0;
  };

  //Set Functions 
  void SetShowerDirection    (TVector3 &showerdirection)              {ShowerDirection     = showerdirection;     ShowerDirectionPtr     = &ShowerDirection;} 
  void SetShowerStartPosition(TVector3 &showerstartposition)          {ShowerStartPosition = showerstartposition; ShowerStartPositionPtr = &ShowerStartPosition;} 
  void SetInitialTrack       (recob::Track &track)                    {InitialTrack        = track;               InitialTrackPtr        = &InitialTrack;}
  void SetInitialTrackHits   (std::vector<art::Ptr<recob::Hit> > hits){InitialTrackHits    = hits;                InitialTrackHitsPtr    = &InitialTrackHits;}      
  void SetShowerEnergy       (std::vector<double> energyvec)          {ShowerEnergy        = energyvec;           ShowerEnergyPtr        = &energyvec;}
  void SetShowerdEdx         (std::vector<double> dedxvec)            {ShowerdEdx          = dedxvec;             ShowerdEdxPtr          = &ShowerdEdx;}
  void SetBestPlane          (int plane)                              {BestPlane           = plane;               BestPlanePtr           = &plane;}
  

  void SetShowerDirectionErr    (TVector3 &showerdirectionerr)    {ShowerDirectionErr     = showerdirectionerr;     ShowerDirectionErrPtr     = &ShowerDirectionErr;} 
  void SetShowerStartPositionErr(TVector3 &showerstartpositionerr){ShowerStartPositionErr = showerstartpositionerr; ShowerStartPositionErrPtr = &ShowerStartPositionErr;} 
  void SetShowerEnergyErr       (std::vector<double> energyvecerr){ShowerEnergyErr        = energyvecerr;           ShowerEnergyErrPtr        = &ShowerEnergyErr;}
  void SetShowerdEdxErr         (std::vector<double> dedxvecerr)  {ShowerdEdxErr          = dedxvecerr;             ShowerdEdxErrPtr          = &ShowerdEdxErr;}


  //Get Functions
  TVector3 GetShowerDirection(){
    if(ShowerDirectionPtr) return ShowerDirection;
    TVector3 EmptyVec = {-999,-999,-999};
    return EmptyVec;
  }
    
  TVector3 GetShowerStartPosition(){
    if(ShowerStartPositionPtr) return ShowerStartPosition;
    TVector3 EmptyVec = {-999,-999,-999};
    return EmptyVec;
  }
  
  recob::Track GetInitialTrack(){
    if(InitialTrackPtr) return InitialTrack;
    return InitialTrack;
  } 

  std::vector<art::Ptr<recob::Hit> > GetInitialTrackHits(){
    if(InitialTrackHitsPtr) return InitialTrackHits;
    std::vector<art::Ptr<recob::Hit> > EmptyVec;
    return  EmptyVec;
  }

  std::vector<double> GetShowerEnergy(){
    if(ShowerEnergyPtr) return ShowerEnergy;
    std::vector<double> EmptyVec = {-999,-999,-999};
    return EmptyVec;
  }

  std::vector<double> GetShowerdEdx(){
    if(ShowerdEdxPtr)return ShowerdEdx;
    std::vector<double> EmptyVec = {-999,-999,-999};
    return EmptyVec;

  }

  int GetBestPlane(){
    if(BestPlanePtr)return BestPlane;
    int Empty = -999;
    return Empty;

  }

  TVector3 GetShowerDirectionErr(){
    if(ShowerDirectionErrPtr) return ShowerDirectionErr;
    TVector3 EmptyVec = {-999,-999,-999};
    return EmptyVec;
  }
    
  TVector3 GetShowerStartPositionErr(){
    if(ShowerStartPositionErrPtr) return ShowerStartPositionErr;
    TVector3 EmptyVec = {-999,-999,-999};
    return EmptyVec;
  }
  
  std::vector<double> GetShowerEnergyErr(){
    if(ShowerEnergyErrPtr) return ShowerEnergyErr;
    std::vector<double> EmptyVec = {-999,-999,-999};
    return EmptyVec;
  }

  std::vector<double> GetShowerdEdxErr(){
    if(ShowerdEdxErrPtr)return ShowerdEdxErr;
    std::vector<double> EmptyVec = {-999,-999,-999};
    return EmptyVec;
  }

  ///Check Functions 
  bool CheckShowerDirection(){
    if(ShowerDirectionPtr)  return true; 
    else return false;
  } 
  bool CheckShowerStartPosition(){
    if(ShowerStartPositionPtr) return true; 
    else return false;
  } 
  
  bool CheckInitialTrack(){
    if(InitialTrackPtr) return true; 
    else return false;
  } 

  bool CheckInitialTrackHits(){
    if(InitialTrackHitsPtr) return true;
    else return false;
  }
  
  bool CheckShowerEnergy(){
      if(ShowerEnergyPtr) return true; 
    else return false;
  } 

  bool CheckShowerdEdx(){
    if(ShowerdEdxPtr) return true; 
    else return false;
  } 

  bool CheckBestPlane(){
    if(BestPlanePtr) return true; 
    else return false;
  } 

  bool CheckShowerDirectionErr(){
    if(ShowerDirectionErrPtr)  return true; 
    else return false;
  } 
  bool CheckShowerStartPositionErr(){
    if(ShowerStartPositionErrPtr) return true; 
    else return false;
  } 
  
  bool CheckShowerEnergyErr(){
    if(ShowerEnergyErrPtr) return true; 
    else return false;
  } 

  bool CheckShowerdEdxErr(){
    if(ShowerdEdxErrPtr) return true; 
    else return false;
  } 

 private:

  //Shower Parameter Ptrs   
  TVector3*                           ShowerDirectionPtr;
  TVector3*                           ShowerStartPositionPtr;
  recob::Track*                       InitialTrackPtr; 
  std::vector<art::Ptr<recob::Hit> >* InitialTrackHitsPtr;
  std::vector<double>*                ShowerEnergyPtr;
  std::vector<double>*                ShowerdEdxPtr;
  int*                                BestPlanePtr;

  //Shower Parameters 
  TVector3                           ShowerDirection;
  TVector3                           ShowerStartPosition;
  recob::Track                       InitialTrack; 
  std::vector<art::Ptr<recob::Hit> > InitialTrackHits;
  std::vector<double>                ShowerEnergy;
  std::vector<double>                ShowerdEdx;
  int                                BestPlane;

  //Shower Parameter Error Ptrs   
  TVector3*            ShowerDirectionErrPtr;
  TVector3*            ShowerStartPositionErrPtr;
  std::vector<double>* ShowerEnergyErrPtr;
  std::vector<double>* ShowerdEdxErrPtr;

  //Shower Parameter Errors
  TVector3            ShowerDirectionErr;
  TVector3            ShowerStartPositionErr;
  std::vector<double> ShowerEnergyErr;
  std::vector<double> ShowerdEdxErr;
  
};


#endif
