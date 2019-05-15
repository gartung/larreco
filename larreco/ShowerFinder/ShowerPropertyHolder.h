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
#include "messagefacility/MessageLogger/MessageLogger.h"

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
      ShowerDirection = 0;
      ShowerStartPosition = 0;
      InitialTrack = 0;
      ShowerEnergy = 0;
      ShowerdEdx = 0;
      err = 0;
  };

  //Reset the holder
  void clear(){
    ShowerDirection = 0;
    ShowerStartPosition = 0;
    InitialTrack = 0;
    ShowerEnergy = 0;
    ShowerdEdx = 0;
    err = 0;
  };

  //Set Functions
  void SetShowerDirection    (TVector3 &showerdirection)    {ShowerDirection     = &showerdirection;} 
  void SetShowerStartPosition(TVector3 &showerstartposition){ShowerStartPosition = &showerstartposition;} 
  void SetInitialTrack       (recob::Track &track)          {InitialTrack        = &track;}
  void SetShowerEnergy       (std::vector<double> energyvec){ShowerEnergy        = &energyvec;}
  void SetShowerdEdx         (std::vector<double> dedxvec)  {ShowerdEdx          = &dedxvec;}

  //Get Functions
  TVector3 GetShowerDirection(){
    if(ShowerDirection) return *ShowerDirection;
    else throw cet::exception("ShowerPropertyHolder") << "Shower Direction not set" << std::endl;      
  }
    
  TVector3 GetShowerStartPosition(){
    if(ShowerStartPosition) return *ShowerStartPosition;
    else throw cet::exception("ShowerPropertyHolder") << "Shower start position not set" << std::endl;
  }
  
  recob::Track GetInitialTrack(){
    if(InitialTrack) return *InitialTrack;
    else throw cet::exception("ShowerPropertyHolder") << "Shower Initial track not set" << std::endl;
  } 

  std::vector<double> GetShowerEnergy(){
    if(ShowerEnergy) return *ShowerEnergy;
    else throw cet::exception("ShowerPropertyHolder") << "Shower Energy not set" << std::endl;
  }

  std::vector<double> GetShowerdEdx(){
    if(ShowerdEdx)return *ShowerdEdx;
    else throw cet::exception("ShowerPropertyHolder") << "Shower dEdx not set" << std::endl;
  }

  //Error 
  int GetErrorStatus(){return err;}

 private:

  //Shower Parameters 
  TVector3*            ShowerDirection;
  TVector3*            ShowerStartPosition;
  recob::Track*        InitialTrack; 
  std::vector<double>* ShowerEnergy;
  std::vector<double>* ShowerdEdx;

  //ErrorCode 
  int err;
  
};


#endif
