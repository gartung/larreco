#include "larreco/RecoAlg/SBNShowerAlg.h"

shower::SBNShowerAlg::SBNShowerAlg(const fhicl::ParameterSet& pset) {
  fUseCollectionOnly = pset.get<bool>("UseCollectionOnly");
}


void shower::SBNShowerAlg::OrderShowerHits(int test){
  std::cout<<"Testing SBNShowerAlg: "<<test<<std::endl;
  return;
}

void shower::SBNShowerAlg::OrderShowerSpacePoints( std::vector<art::Ptr<recob::SpacePoint> >& 
						   showersps, TVector3& vertex, 
						   TVector3& direction){

  std::map<double,art::Ptr<recob::SpacePoint> > OrderedSpacePoints;

  //Loop over the spacepoints and get the pojected distance from the vertex.                       
  for(auto const& sp: showersps){

    //Get the position of the spacepoint                                                           
    TVector3 pos = shower::SBNShowerAlg::SpacePointPosition(sp) - vertex;

    //Get the the projected length                                                                 
    double len = pos.Dot(direction);

    //Add to the list                                                                              
    OrderedSpacePoints[len] = sp;
  }

  //Return an ordered list.                                                                        
  showersps.clear();
  for(auto const& sp: OrderedSpacePoints){
    showersps.push_back(sp.second);
  }
  return;
}

TVector3 shower::SBNShowerAlg::ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& 
					    showerspcs, art::FindManyP<recob::Hit>& fmh){
  float totalCharge=0;
  TVector3 centre =  shower::SBNShowerAlg::ShowerCentre(showerspcs,fmh,totalCharge);
  return centre;

}


TVector3 shower::SBNShowerAlg::ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& showersps,
					    art::FindManyP<recob::Hit>& fmh, float& totalCharge){

  TVector3 pos, chargePoint = TVector3(0,0,0);

  //Loop over the spacepoints and get the charge weighted center.
  for(auto const& sp: showersps){

    //Get the position of the spacepoint
    pos = SpacePointPosition(sp);

    //Get the associated hits
    std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());

    //Average the charge unless sepcified.
    float charge  = 0;
    float charge2 = 0;
    for(auto const& hit: hits){

      if(fUseCollectionOnly){
	if(hit->SignalType() == geo::kCollection){
	  charge = hit->Integral();
	  break;
	}
      }

      //Check if any of the points are not withing 2 sigma.
      if(!fUseCollectionOnly){
	charge += hit->Integral();
	charge2 += hit->Integral();
      }
    }

    if(!fUseCollectionOnly){
      //Calculate the unbiased standard deviation and mean.
      float mean = charge/((float) hits.size());
      float rms  = TMath::Sqrt((charge2 - charge*charge)/((float)(hits.size()-1)));

      charge = 0;
      for(auto const& hit: hits){
	if(hit->Integral() > (mean - 2*rms) && hit->Integral() < (mean + 2*rms))
	  charge += hit->Integral();
      }
    }

    chargePoint += charge * pos;
    totalCharge += charge;

    if(charge == 0){
      mf::LogWarning("ShowerStartPosition") << 
	"Averaged charge, within 2 sigma, for a spacepoint is zero, Maybe this not a good method";
    }
  }
    
  double intotalcharge = 1/totalCharge;
  TVector3 centre = chargePoint *  intotalcharge;
  return centre;
 
}


TVector3 shower::SBNShowerAlg::SpacePointPosition(const art::Ptr<recob::SpacePoint>& sp){

  const Double32_t* sp_xyz = sp->XYZ();
  TVector3 sp_postiion = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
  return sp_postiion;
}

double shower::SBNShowerAlg::SpacePointCharge(art::Ptr<recob::SpacePoint> sp,
					      art::FindManyP<recob::Hit>& fmh){

  return 1;
}


