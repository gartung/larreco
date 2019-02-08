#include "ShowerUtils.h"

std::pair<int,double> ShowerUtils::TrueParticleIDFromTrueChain(std::map<int,std::vector<int>> &ShowersMothers,const std::vector<art::Ptr<recob::Hit> >& hits, int planeid) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  //Find the energy for each track ID.
  std::map<int,double> trackIDToEDepMap;
  std::map<int,double> trackIDTo3EDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;

    //Get the plane ID                                                                                
    geo::WireID wireid = (*hitIt)->WireID();
    int PlaneID = wireid.Plane;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      trackIDTo3EDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;
      if(PlaneID == planeid){trackIDToEDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;}
    }
  }

  //Find the energy for each showermother.
  std::map<int,double> MotherIDtoEMap; 
  std::map<int,double> MotherIDto3EMap;
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
    for(std::vector<int>::iterator daughter=(showermother->second).begin(); daughter!=(showermother->second).end(); ++daughter){
      MotherIDtoEMap[showermother->first]  +=  trackIDToEDepMap[*daughter];
      MotherIDto3EMap[showermother->first] +=  trackIDTo3EDepMap[*daughter];
    }
  }

  //Loop over the mothers to find the most like candidate by identifying which true shower deposited the most energy in the hits. 
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator mapIt = MotherIDto3EMap.begin(); mapIt != MotherIDto3EMap.end(); mapIt++){
    double energy_three = mapIt->second;
    double trackid = mapIt->first;
    if (energy_three > maxenergy){
      maxenergy = energy_three;
      objectTrack = trackid;
    }
  }

  //If the none of the shower mother deposited no energy then we cannot match this.
  if(maxenergy == 0){
    return std::make_pair(-99999,-99999);
  }

  return std::make_pair(objectTrack,MotherIDtoEMap[objectTrack]);
}


std::map<geo::PlaneID,int> ShowerUtils::NumberofWiresHitByShower(std::vector<int> &TrackIDs, const std::vector<art::Ptr<recob::Hit> >& hits){ 
  
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  std::vector<geo::WireID> WiresUsed;
  std::map<geo::PlaneID,int> HitWirePlaneMap;

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;

    //Find the wire and  plane id                                                                                
    geo::WireID wireid = hit->WireID();
    geo::PlaneID  PlaneID = wireid.planeID();

    //Check to see if the wire is already been continued.
    if(std::find(WiresUsed.begin(),WiresUsed.end(),wireid) != WiresUsed.end()){continue;}

    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      if(std::find(TrackIDs.begin(), TrackIDs.end(),TMath::Abs(trackIDs.at(idIt).trackID)) != TrackIDs.end()){ 
	WiresUsed.push_back(wireid);
	++HitWirePlaneMap[PlaneID];
	break;
      }
    }
  }
  return HitWirePlaneMap;
}
