#include "ShowerUtils.h"

std::pair<int,double> ShowerUtils::TrueParticleIDFromTrueChain(std::map<int,std::vector<int>> &ShowersMothers,const std::vector<art::Ptr<recob::Hit> >& hits, int planeid) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  //Find the energy for each track ID.
  std::map<int,double> trackIDToEDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;

    //Get the plane ID                                                                                                                                                           
    geo::WireID wireid = (*hitIt)->WireID();
    int PlaneID = wireid.Plane;
    if(PlaneID != planeid){continue;}
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      trackIDToEDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;
    }
  }

  //Find the energy for each showermother.
  std::map<int,double> MotherIDtoEMap; 
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
    for(std::vector<int>::iterator daughter=(showermother->second).begin(); daughter!=(showermother->second).end(); ++daughter){
      MotherIDtoEMap[showermother->first] +=  trackIDToEDepMap[*daughter];
    }
  }

  //Loop over the mothers to find the most like candiate. 
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator mapIt = MotherIDtoEMap.begin(); mapIt != MotherIDtoEMap.end(); mapIt++){
    double energy = mapIt->second;
    double trackid = mapIt->first;
    if (energy > maxenergy){
      maxenergy = energy;
      objectTrack = trackid;
      
    }
  }
  return std::make_pair(objectTrack,maxenergy);
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
