////////////////////////////////////////////////////////////////////
// Implementation of the Blurred Clustering algorithm
//
// Converts a hit map into a 2D image of the hits before convoling
// with a Gaussian function to introduce a weighted blurring.
// Clustering proceeds on this blurred image to create more
// complete clusters.
//
// M Wallbank (m.wallbank@sheffield.ac.uk), May 2015
////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/BlurredClusterAlg.h"

cluster::BlurredClusterAlg::BlurredClusterAlg(fhicl::ParameterSet const& pset) {

  this->reconfigure(pset);
  this->MakeKernels();

  // For the debug PDF
  fDebugCanvas = NULL;
  fDebugPDFName = "";

}

cluster::BlurredClusterAlg::~BlurredClusterAlg() {
  if (fDebugCanvas) {
    std::string closeName = fDebugPDFName;
    closeName.append("]");
    fDebugCanvas->Print(closeName.c_str());
    delete fDebugCanvas;
  }
}

void cluster::BlurredClusterAlg::reconfigure(fhicl::ParameterSet const& p) {
  fBlurWire             = p.get<int>   ("BlurWire");
  fBlurTick             = p.get<int>   ("BlurTick");
  fSigmaWire            = p.get<double>("SigmaWire");
  fSigmaTick            = p.get<double>("SigmaTick");
  fMaxTickWidthBlur     = p.get<int>   ("MaxTickWidthBlur");
  fShowerDirectionWidth = p.get<double>("ShowerDirectionWidth");
  fClusterWireDistance  = p.get<int>   ("ClusterWireDistance");
  fClusterTickDistance  = p.get<int>   ("ClusterTickDistance");
  fNeighboursThreshold  = p.get<int>   ("NeighboursThreshold");
  fMinNeighbours        = p.get<int>   ("MinNeighbours");
  fMinSize              = p.get<int>   ("MinSize");
  fMinSeed              = p.get<double>("MinSeed");
  fChargeThreshold      = p.get<double>("ChargeThreshold");
  fDebug                = p.get<bool>  ("Debug",false);
  fDetector             = p.get<std::string>("Detector");

  fKernelWidth = 2 * fBlurWire + 1;
  fKernelHeight = 2 * fBlurTick*fMaxTickWidthBlur + 1;

  fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

int cluster::BlurredClusterAlg::BinKey(const std::vector<std::vector<double> >& image, int bin) {

  art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, bin);

  return hit.isNull() ? -999 : hit.key();

}

void cluster::BlurredClusterAlg::CreateDebugPDF(int run, int subrun, int event) {

  if (!fDebugCanvas) {

    // Create the grayscale palette for the Z axis
    Double_t Red[2] = { 1.00, 0.00 };
    Double_t Green[2] = { 1.00, 0.00 };
    Double_t Blue[2] = { 1.00, 0.00 };
    Double_t Length[2] = { 0.00, 1.00 };
    TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, 1000);
    gStyle->SetOptStat(110000);

    // Decide what to call this PDF
    std::ostringstream oss;
    oss << "BlurredImages_Run" << run << "_Subrun" << subrun;
    fDebugPDFName = oss.str();
    fDebugCanvas = new TCanvas(fDebugPDFName.c_str(), "Image canvas", 1000, 500);
    fDebugPDFName.append(".pdf");

    std::string openName = fDebugPDFName;
    openName.append("[");
    fDebugCanvas->Print(openName.c_str());
    fDebugCanvas->Divide(2, 2);
    fDebugCanvas->SetGrid();
  }

  // Clear the pads on the canvas
  for (int i = 1; i <= 4; i++) {
    fDebugCanvas->GetPad(i)->Clear();
  }

  std::ostringstream oss;
  oss << "Event " << event;
  fDebugCanvas->cd(1);
  TLatex l;
  l.SetTextSize(0.15);
  l.DrawLatex(0.1, 0.1, oss.str().c_str());
  fDebugCanvas->Print(fDebugPDFName.c_str());

}

TVector2 cluster::BlurredClusterAlg::Convert3DPointToPlaneBins(const TVector3& point, int plane) {

  // Project the 3D point onto the correct plane
  TVector2 wireTickPos = Project3DPointOntoPlane(point, plane);

  wireTickPos -= TVector2(fLowerWire, fLowerTick);

  return wireTickPos;

  // // Now convert wire/tick to bin number
  // int bin = ConvertWireTickBinsToBin(image, wireTickPos.X()-fLowerWire, wireTickPos.Y()-fLowerTick);

  // return bin;

}

TVector2 cluster::BlurredClusterAlg::ConvertBinTo2DPosition(const std::vector<std::vector<double > >& image, int bin) {

  art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, bin);
  return HitPosition(hit);

}

art::PtrVector<recob::Hit> cluster::BlurredClusterAlg::ConvertBinsToRecobHits(const std::vector<std::vector<double> >& image, std::vector<int> const& bins) {

  // Create the vector of hits to output
  art::PtrVector<recob::Hit> hits;

  // Look through the hits in the cluster
  for (std::vector<int>::const_iterator binIt = bins.begin(); binIt != bins.end(); binIt++) {

    // Take each hit and convert it to a recob::Hit
    art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, *binIt);

    // If this hit was a real hit put it in the hit selection
    if (!hit.isNull())
      hits.push_back(hit);
  }

  // Return the vector of hits to make cluster
  return hits;
}

art::Ptr<recob::Hit> cluster::BlurredClusterAlg::ConvertBinToRecobHit(const std::vector<std::vector<double> >& image, int bin) {

  int wire = bin % image.size();
  int tick = bin / image.size();

  return fHitMap[wire][tick];

}

void cluster::BlurredClusterAlg::ConvertBinsToClusters(std::vector<std::vector<double> > const& image,
						       std::vector<std::vector<int> > const& allClusterBins,
						       std::vector<art::PtrVector<recob::Hit> >& clusters) {

  // Loop through the clusters (each a vector of bins)
  for (std::vector<std::vector<int> >::const_iterator clustIt = allClusterBins.begin(); clustIt != allClusterBins.end(); clustIt++) {
    std::vector<int> bins = *clustIt;

    // Convert the clusters (vectors of bins) to hits in a vector of recob::Hits
    art::PtrVector<recob::Hit> clusHits = ConvertBinsToRecobHits(image, bins);

    mf::LogInfo("BlurredCluster") << "Cluster made from " << bins.size() << " bins, of which " << clusHits.size() << " were real hits";

    // Make sure the clusters are above the minimum cluster size
    if (clusHits.size() < fMinSize) {
      mf::LogVerbatim("BlurredCluster") << "Cluster of size " << clusHits.size() << " not saved since it is smaller than the minimum cluster size, set to " << fMinSize;
      continue;
    }

    clusters.push_back(clusHits);

  }

  return;

}

void cluster::BlurredClusterAlg::ConvertBinToWireTickBins(int bin, const std::vector<std::vector<double> >& blurred, int& xbin, int& ybin) {

  int nbinsx = blurred.size();
  int nbinsy = blurred.at(0).size();

  xbin = bin % nbinsx;
  ybin = ((bin - xbin) / nbinsx) % nbinsy;

  return;

}

std::vector<std::vector<double> > cluster::BlurredClusterAlg::ConvertRecobHitsToVector(std::vector<art::Ptr<recob::Hit> > const& hits) {

  // Define the size of this particular plane -- dynamically to avoid huge histograms
  int lowerTick = fDetProp->ReadOutWindowSize(), upperTick = 0, lowerWire = fGeom->MaxWires(), upperWire = 0;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    TVector2 hitCoords = HitCoordinates(*hitIt);
    if (hitCoords.Y() < lowerTick) lowerTick = hitCoords.Y();
    if (hitCoords.Y() > upperTick) upperTick = hitCoords.Y();
    if (hitCoords.X() < lowerWire) lowerWire = hitCoords.X();
    if (hitCoords.X() > upperWire) upperWire = hitCoords.X();
    //std::cout << "Hit at " << (*hitIt)->WireID() << " has peak time " << (*hitIt)->PeakTime() << " and RMS " << (*hitIt)->RMS() << ": lower and upper tick is " << (*hitIt)->StartTick() << " and " << (*hitIt)->EndTick() << ")" << std::endl;
  }
  fLowerTick = lowerTick-20;
  fUpperTick = upperTick+20;
  fLowerWire = lowerWire-20;
  fUpperWire = upperWire+20;

  // Use a map to keep a track of the real hits and their wire/ticks
  fHitMap.clear();
  fHitMap.resize(fUpperWire-fLowerWire, std::vector<art::Ptr<recob::Hit> >(fUpperTick-fLowerTick, art::Ptr<recob::Hit>()));

  // Create a 2D vector
  std::vector<std::vector<double> > image(fUpperWire-fLowerWire, std::vector<double>(fUpperTick-fLowerTick, 0));

  // Look through the hits
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    int wire = GlobalWire((*hitIt)->WireID());
    int tick = (int)(*hitIt)->PeakTime();
    float charge = (*hitIt)->Integral();

    // Fill hit map and keep a note of all real hits for later
    if (charge > image.at(wire-fLowerWire).at(tick-fLowerTick)) {
      image.at(wire-fLowerWire).at(tick-fLowerTick) = charge;
      fHitMap[wire-fLowerWire][tick-fLowerTick] = (*hitIt);
      fHits.push_back(std::make_pair(hitIt->key(), std::make_pair(HitPosition(*hitIt), (*hitIt)->Integral())));
    }
  }

  // Keep a note of dead wires
  fDeadWires.clear();
  fDeadWires.resize(fUpperWire-fLowerWire,false);
  geo::PlaneID planeID = hits.front()->WireID().planeID();

  for (int wire = fLowerWire; wire < fUpperWire; ++wire) {
    raw::ChannelID_t channel = fGeom->PlaneWireToChannel(planeID.Plane,wire,planeID.TPC,planeID.Cryostat);
    fDeadWires[wire-fLowerWire] = !fChanStatus.IsGood(channel);
  }

  return image;

}

int cluster::BlurredClusterAlg::ConvertWireTickBinsToBin(std::vector<std::vector<double> > const& image, int xbin, int ybin) {

  return (ybin * image.size()) + xbin;

}

double cluster::BlurredClusterAlg::ConvertBinToCharge(std::vector<std::vector<double> > const& image, int bin) {

  int x = bin % image.size();
  int y = bin / image.size();

  return image.at(x).at(y);

}

std::pair<int,int> cluster::BlurredClusterAlg::DeadWireCount(int wire_bin, int width) {

  std::pair<int,int> deadWires = std::make_pair<int,int>(0,0);

  int lower_bin = width / 2;
  int upper_bin = (width+1) / 2;

  for (int wire = TMath::Max(wire_bin + fLowerWire - lower_bin, fLowerWire); wire < TMath::Min(wire_bin + fLowerWire + upper_bin, fUpperWire); ++wire) {
    if (fDeadWires[wire-fLowerWire]) {
      if (wire < wire_bin + fLowerWire)
	++deadWires.first;
      else if (wire > wire_bin + fLowerWire)
	++deadWires.second;
    }
  }

  return deadWires;

}

BlurringParameters cluster::BlurredClusterAlg::FindBlurringParameters(const std::vector<std::pair<int,std::pair<TVector2,double> > >& hits,
								      const TVector2& point,
								      const std::vector<bool>& used) {

  // Define a unit vector for the rough shower direction
  const TVector2 unit = TVector2(0,1).Unit();
  TVector2 direction;

  // Reduce the hit container to just ununsed hits
  std::vector<std::pair<TVector2,double> > unusedHits;
  for (std::vector<std::pair<int,std::pair<TVector2,double> > >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
    if (!used[hitIt->first])
      unusedHits.push_back(hitIt->second);

  // Look over a range of angles and save the direction vector for the angle
  // which contains most of the charge within a certain rectangle width
  TVector2 pos, proj;
  double max_charge = 0;
  for (int angle = 0; angle < 360; angle+=5) {
    TVector2 rotated = unit.Rotate(angle*TMath::Pi()/180);
    int charge = 0;
    for (std::vector<std::pair<TVector2,double> >::const_iterator hitIt = unusedHits.begin(); hitIt != unusedHits.end(); ++hitIt) {
      pos = hitIt->first - point;
      proj = pos.Proj(rotated);
      if ((pos-proj).Mod() < fShowerDirectionWidth)
	charge += hitIt->second;
    }
    if (charge > max_charge) {
      max_charge = charge;
      direction = rotated;
    }
  }

  return FindBlurringParameters(direction);

}

BlurringParameters cluster::BlurredClusterAlg::FindBlurringParameters() {

  // Calculate least squares slope
  int x, y;
  double nhits=0, sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (unsigned int wireIt = 0; wireIt < fHitMap.size(); ++wireIt) {
    for (unsigned int tickIt = 0; tickIt < fHitMap.at(wireIt).size(); ++tickIt) {
      if (fHitMap[wireIt][tickIt].isNull())
	continue;
      ++nhits;
      x = wireIt + fLowerWire;
      y = tickIt + fLowerTick;
      sumx += x;
      sumy += y;
      sumx2 += x*x;
      sumxy += x*y;
    }
  }
  double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);

  // Get the rough unit vector for the trajectories
  TVector2 unit = TVector2(1,gradient).Unit();

  // Catch vertical gradients
  if (std::isnan(gradient))
    unit = TVector2(0,1);

  return FindBlurringParameters(unit);

}

BlurringParameters cluster::BlurredClusterAlg::FindBlurringParameters(const TVector2& direction) {

  // Use this direction to scale the blurring radii and Gaussian sigma
  int blur_wire = std::max(std::abs(std::round(fBlurWire * direction.X())),1.);
  int blur_tick = std::max(std::abs(std::round(fBlurTick * direction.Y())),1.);

  int sigma_wire = std::max(std::abs(std::round(fSigmaWire * direction.X())),1.);
  int sigma_tick = std::max(std::abs(std::round(fSigmaTick * direction.Y())),1.);

  // std::cout << "  Blurring: direction (" << direction.X() << ", " << direction.Y() << ")"
  // 	    << "; wire " << blur_wire << " and tick " << blur_tick << "; sigma: wire " << sigma_wire << " and tick " << sigma_tick << std::endl;
  
  return BlurringParameters(blur_wire, blur_tick, sigma_wire, sigma_tick);

}

std::vector<std::vector<int> > cluster::BlurredClusterAlg::FindClusters(const std::vector<std::vector<double> >& hit_map, const std::vector<TVector2>& vertices, bool reblur) {

  // Containers to hold output clusters
  std::vector<int> cluster;
  std::vector<std::vector<int> > allClusters;

  // Size of image in x and y
  const int nbinsx = hit_map.size();
  const int nbinsy = hit_map.at(0).size();

  // Keep an note of which bins have been clustered
  std::vector<std::vector<bool> > used(nbinsx, std::vector<bool>(nbinsy, false));
  std::vector<bool> usedHits(nbinsx*nbinsy, false);

  // Get list of bins and charges
  std::vector<std::pair<double,std::pair<int,int> > > bin_charges;
  for (int xbin = 0; xbin < nbinsx; ++xbin)
    for (int ybin = 0; ybin < nbinsy; ++ybin)
      bin_charges.push_back(std::make_pair(hit_map[xbin][ybin], std::make_pair(xbin, ybin)));
  std::sort(bin_charges.rbegin(), bin_charges.rend());

  // Hit map to hold the blurred image
  std::vector<std::vector<double> > blurred, reblurred;

  int niter = 0;

  // Clustering loops
  // First loop - considers highest charge hits in decreasing order, and puts them in a new cluster if they aren't already clustered (makes new cluster every iteration)
  // Second loop - looks at the direct neighbours of this seed and clusters to this if above charge thresholds. Runs recursively over all hits in cluster (inc. new ones)
  bool createNewClusters = true;
  while (createNewClusters) {

    // Start a new cluster each time loop is executed
    cluster.clear();

    // Find the highest charged unclustered bin
    double bin_charge = bin_charges[niter].first;
    if (bin_charge < fMinSeed) {
      createNewClusters = false;
      break;
    }
    std::pair<int,int> bins = bin_charges[niter++].second;
    int xbin = bins.first, ybin = bins.second;
    if (used[xbin][ybin])
      continue;

    // Reblur the hit map around this bin if necessary
    int bin = ConvertWireTickBinsToBin(hit_map, xbin, ybin);
    if (reblur)
      reblurred = GaussianBlur(hit_map, bin, usedHits);
    const std::vector<std::vector<double> >& blurred = reblur ? reblurred : hit_map;

    // Start a new cluster
    cluster.push_back(bin);
    used[xbin][ybin] = true;
    if (reblur) {
      int key = BinKey(blurred, bin);
      if (key >= 0)
	usedHits[key] = true;
    }

    // Now cluster neighbouring hits to this seed
    bool growCluster = true;
    while (growCluster) {

      int nadded = 0;

      for (unsigned int clusBin = 0; clusBin < cluster.size(); ++clusBin) {

	// Get x and y bins
	ConvertBinToWireTickBins(cluster[clusBin], blurred, xbin, ybin);

	// Look for hits in the neighbouring x/y bins
        for (int x = xbin - fClusterWireDistance; x <= xbin + fClusterWireDistance; ++x) {
          for (int y = ybin - fClusterTickDistance; y <= ybin + fClusterTickDistance; ++y) {
            if (used[x][y] or (x == xbin and y == ybin) or (x >= nbinsx or y >= nbinsy) or (x < 0 or y < 0))
              continue;

	    // Check we're not clustering across a vertex
	    bool crossVertex = false;
	    for (std::vector<TVector2>::const_iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
	      if (((vertexIt->X() - xbin) * (vertexIt->X() - x)) <= 0 or
		  ((vertexIt->Y() - ybin) * (vertexIt->Y() - y)) <= 0)
		crossVertex = true;
	    if (crossVertex)
	      continue;

	    // Cluster bin if above charge threshold
	    bin = ConvertWireTickBinsToBin(blurred, x, y);
	    // if (bin >= nbinsx * nbinsy or bin < 0)
            //   continue;
            if (ConvertBinToCharge(blurred, bin) > fChargeThreshold) {
              cluster.push_back(bin);
              used[x][y] = true;
	      if (reblur) {
		int key = BinKey(blurred, bin);
		if (key >= 0)
		  usedHits[key] = true;
	      }
              ++nadded;
	    }

          }
	} // End of looking at directly neighbouring bins

      } // End of looping over bins already in this cluster

      if (nadded == 0)
        growCluster = false;

    } // End of adding hits to this cluster

    // Check this cluster is above minimum size
    if (cluster.size() < fMinSize) {
      for (unsigned int clusterBin = 0; clusterBin < cluster.size(); clusterBin++) {
	ConvertBinToWireTickBins(cluster[clusterBin], blurred, xbin, ybin);
        used[xbin][ybin] = false;
      }
      continue;
    }

    // Fill in holes in the cluster
    for (unsigned int clusBin = 0; clusBin < cluster.size(); clusBin++) {

      // Looks at directly neighbouring bins (and not itself)
      for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
          if (x == 0 and y == 0) continue;

	  // Look at neighbouring bins to the clustered bin which are inside the cluster
          int neighbouringBin = cluster[clusBin] + x + (y * nbinsx);
          if (neighbouringBin < nbinsx or neighbouringBin % nbinsx == 0 or neighbouringBin % nbinsx == nbinsx - 1 or neighbouringBin >= nbinsx * (nbinsy - 1))
            continue;

	  // If not already clustered and passes neighbour thresholds, add to cluster
	  ConvertBinToWireTickBins(neighbouringBin, blurred, xbin, ybin);
          if (!used[xbin][ybin] and (NumNeighbours(xbin, ybin, used) > fNeighboursThreshold) and RealHit(blurred, neighbouringBin)) {
            used[xbin][ybin] = true;
	    if (reblur) {
	      int key = BinKey(blurred, neighbouringBin);
	      if (key >= 0)
		usedHits[key] = true;
	    }
            cluster.push_back(neighbouringBin);
          }

        }
      } // End of looping over neighbouring bins

    } // End of looping over bins already in cluster

    mf::LogVerbatim("Blurred Clustering") << "Size of cluster after filling in holes: " << cluster.size();


    // Remove peninsulas - hits which have too few neighbouring hits in the cluster (defined by fMinNeighbours)
    bool peninsulasToRemove = true;
    while (peninsulasToRemove) {

      int nremoved = 0;

      // Loop over all the bins in the cluster
      for (int clusBin = cluster.size() - 1; clusBin >= 0; clusBin--) {
        bin = cluster[clusBin];

	// If bin is in cluster ignore
        if (bin < nbinsx || bin % nbinsx == 0 || bin % nbinsx == nbinsx - 1 || bin >= nbinsx * (nbinsy - 1))
	  continue;

	// Remove hit if it has too few neighbouring hits
	ConvertBinToWireTickBins(bin, blurred, xbin, ybin);
        if ((int) NumNeighbours(xbin, ybin, used) < fMinNeighbours) {
          used[xbin][ybin] = false;
          nremoved++;
          cluster.erase(cluster.begin() + clusBin);
        }
      }

      if (nremoved == 0)
	peninsulasToRemove = false;

    }

    mf::LogVerbatim("Blurred Clustering") << "Size of cluster after removing peninsulas: " << cluster.size();

    // Disregard cluster if not of minimum size
    if (cluster.size() < fMinSize) {
      for (unsigned int clusterBin = 0; clusterBin < cluster.size(); clusterBin++) {
	ConvertBinToWireTickBins(cluster[clusterBin], blurred, xbin, ybin);
        used[xbin][ybin] = false;
      }
      continue;
    }

    // Put this cluster in the vector of clusters
    allClusters.push_back(cluster);

  } // End loop over this cluster

  return allClusters;

}

TVector2 cluster::BlurredClusterAlg::HitCoordinates(const art::Ptr<recob::Hit>& hit) {

  return TVector2(GlobalWire(hit->WireID()), hit->PeakTime());

}

TVector2 cluster::BlurredClusterAlg::BinCoordinates(int bin, const std::vector<std::vector<double> >& image) {

  int wire_bin, tick_bin;
  ConvertBinToWireTickBins(bin, image, wire_bin, tick_bin);
  return TVector2(wire_bin+fLowerWire, tick_bin+fLowerTick);

}

TVector2 cluster::BlurredClusterAlg::HitPosition(const art::Ptr<recob::Hit>& hit) {

  geo::PlaneID planeID = hit->WireID().planeID();

  return HitPosition(HitCoordinates(hit), planeID);

}

TVector2 cluster::BlurredClusterAlg::HitPosition(const TVector2& pos, geo::PlaneID planeID) {

  return TVector2(pos.X() * fGeom->WirePitch(planeID),
		  fDetProp->ConvertTicksToX(pos.Y(), planeID));

}

TVector2 cluster::BlurredClusterAlg::BinPosition(int bin, const std::vector<std::vector<double> >& image) {

  TVector2 binCoords = BinCoordinates(bin, image);
  return HitPosition(binCoords, fDefaultPlaneID);

}

std::vector<std::vector<double> > cluster::BlurredClusterAlg::GaussianBlur(std::vector<std::vector<double> > const& image, int bin, const std::vector<bool>& used) {

  if (fSigmaWire == 0 and fSigmaTick == 0)
    return image;

  int wire_bin, tick_bin;
  ConvertBinToWireTickBins(bin, image, wire_bin, tick_bin);
  //std::cout << "About to reblur around wire " << wire_bin+fLowerWire << ", tick " << tick_bin+fLowerTick << std::endl;

  // Find the blurring parameters
  BlurringParameters blurring_parameters;
  if (bin == -1)
    blurring_parameters = FindBlurringParameters();
  else
    blurring_parameters = FindBlurringParameters(fHits, ConvertBinTo2DPosition(image, bin), used);

  // Convolve the Gaussian
  int width = 2 * blurring_parameters.blur_wire + 1;
  int height = 2 * blurring_parameters.blur_tick + 1;
  int nbinsx = image.size();
  int nbinsy = image.at(0).size();

  // Blurred histogram and normalisation for each bin
  std::vector<std::vector<double> > blurred(nbinsx, std::vector<double>(nbinsy, 0));

  // Loop through all the bins in the histogram to blur
  for (int x = 0; x < nbinsx; ++x) {
    for (int y = 0; y < nbinsy; ++y) {

      if (image[x][y] == 0)
      	continue;

      // Scale the tick blurring based on the width of the hit
      int tick_scale = TMath::Sqrt(TMath::Power(fHitMap[x][y]->RMS(),2) + TMath::Power(blurring_parameters.sigma_tick,2)) / (double)blurring_parameters.sigma_tick;
      tick_scale = TMath::Max(TMath::Min(tick_scale,fMaxTickWidthBlur),1);
      std::vector<double> correct_kernel = fAllKernels[blurring_parameters.sigma_wire][blurring_parameters.sigma_tick*tick_scale];

      // Find any dead wires in the potential blurring region
      std::pair<int,int> num_deadwires = DeadWireCount(x, width);
      //num_deadwires = std::make_pair<int,int>(0,0);

      // Note of how many dead wires we have passed whilst blurring in the wire direction
      // If blurring below the seed hit, need to keep a note of how many dead wires to come
      // If blurring above, need to keep a note of how many dead wires have passed
      int dead_wires_passed = num_deadwires.first;

      // bool dead = false;
      // if (num_deadwires.first != 0 or num_deadwires.second != 0)
      // 	dead = true;

      // if (dead) {
      // 	std::cout << "Wire is " << x+fLowerWire << std::endl;
      // 	std::cout << "Width is " << width << std::endl;
      // }

      // Loop over the blurring region around this hit
      for (int blurx = -(width/2+num_deadwires.first); blurx < (width+1)/2+num_deadwires.second; ++blurx) {
  	for (int blury = -height/2*tick_scale; blury < ((((height+1)/2)-1)*tick_scale)+1; ++blury) {

	  // if (dead)
	  //   std::cout << "Start... dead_wires_passed is " << dead_wires_passed << " and blurx is " << blurx << std::endl;

	  if (blurx < 0 and fDeadWires[x+blurx])
	    dead_wires_passed -= 1;

  	  // Smear the charge of this hit
	  double weight = correct_kernel[fKernelWidth * (fKernelHeight / 2 + blury) + (fKernelWidth / 2 + (blurx - dead_wires_passed))];
  	  if (x + blurx >= 0 and x + blurx < nbinsx and y + blury >= 0 and y + blury < nbinsy)
  	    blurred[x+blurx][y+blury] += weight * image[x][y];

	  if (blurx > 0 and fDeadWires[x+blurx])
	    dead_wires_passed += 1;

	  // if (dead)
	  //   std::cout << "Start... dead_wires_passed is " << dead_wires_passed << " and blurx is " << blurx << std::endl;

  	}
      } // blurring region

    }
  } // hits to blur

  // HAVE REMOVED NOMALISATION CODE
  // WHEN USING DIFFERENT KERNELS, THERE'S NO EASY WAY OF DOING THIS...
  // RECONSIDER...

  // Return the blurred histogram
  return blurred;

}

int cluster::BlurredClusterAlg::GlobalWire(const geo::WireID& wireID) {

  double globalWire = -999;

  // Induction
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    double wireCentre[3];
    fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }

  // Collection
  else {
    // FOR COLLECTION WIRES, HARD CODE THE GEOMETRY FOR GIVEN DETECTORS
    // THIS _SHOULD_ BE TEMPORARY. GLOBAL WIRE SUPPORT IS BEING ADDED TO THE LARSOFT GEOMETRY AND SHOULD BE AVAILABLE SOON
    if (fDetector == "dune35t") {
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      if (wireID.TPC == 0 or wireID.TPC == 1) globalWire = wireID.Wire;
      else if (wireID.TPC == 2 or wireID.TPC == 3 or wireID.TPC == 4 or wireID.TPC == 5) globalWire = nwires + wireID.Wire;
      else if (wireID.TPC == 6 or wireID.TPC == 7) globalWire = (2*nwires) + wireID.Wire;
      else mf::LogError("BlurredClusterAlg") << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC << " (geometry " << fDetector << ")";
    }
    else if (fDetector == "dune10kt") {
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      // Detector geometry has four TPCs, two on top of each other, repeated along z...
      int block = wireID.TPC / 4;
      globalWire = (nwires*block) + wireID.Wire;
    }
    else {
      double wireCentre[3];
      fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);
      if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
      else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
    }
  }

  return std::round(globalWire);

}

void cluster::BlurredClusterAlg::MakeKernels() {

  // Kernel size is the largest possible given the hit width rescaling
  fAllKernels.clear();
  fAllKernels.resize(fSigmaWire+1,std::vector<std::vector<double> >(fSigmaTick*fMaxTickWidthBlur+1,std::vector<double>(fKernelWidth*fKernelHeight,0.)));

  // Ranges of kernels to make
  // Complete range of sigmas possible after dynamic fixing and hit width convolution
  for (int sigma_wire = 1; sigma_wire <= fSigmaWire; ++sigma_wire) {
    for (int sigma_tick = 1; sigma_tick <= fSigmaTick*fMaxTickWidthBlur; ++sigma_tick) {

      // New kernel
      std::vector<double> kernel(fKernelWidth*fKernelHeight,0);

      // Smear out according to the blur radii in each direction
      for (int i = -fBlurWire; i <= fBlurWire; i++) {
	for (int j = -fBlurTick*fMaxTickWidthBlur; j <= fBlurTick*fMaxTickWidthBlur; j++) {

	  // Fill kernel
	  double sig2i = 2. * sigma_wire * sigma_wire;
	  double sig2j = 2. * sigma_tick * sigma_tick;

	  int key = (fKernelWidth * (j + fBlurTick*fMaxTickWidthBlur)) + (i + fBlurWire);
	  double value = 1. / sqrt(sig2i * M_PI) * exp(-i * i / sig2i) * 1. / sqrt(sig2j * M_PI) * exp(-j * j / sig2j);
	  kernel.at(key) = value;

	}
      } // End loop over blurring region

      fAllKernels[sigma_wire][sigma_tick] = kernel;

    }
  }

  return;

}

TH2F* cluster::BlurredClusterAlg::MakeHistogram(std::vector<std::vector<double> > const& image, TString name) {

  TH2F* hist = new TH2F(name,name,fUpperWire-fLowerWire,fLowerWire-0.5,fUpperWire-0.5,fUpperTick-fLowerTick,fLowerTick-0.5,fUpperTick-0.5);
  hist->Clear();
  hist->SetXTitle("Wire number");
  hist->SetYTitle("Tick number");
  hist->SetZTitle("Charge");

  for (unsigned int imageWireIt = 0; imageWireIt < image.size(); ++imageWireIt) {
    int wire = imageWireIt + fLowerWire;
    for (unsigned int imageTickIt = 0; imageTickIt < image.at(imageWireIt).size(); ++imageTickIt) {
      int tick = imageTickIt + fLowerTick;
      hist->Fill(wire, tick, image.at(imageWireIt).at(imageTickIt));
    }
  }

  return hist;

}

unsigned int cluster::BlurredClusterAlg::NumNeighbours(int xbin, int ybin, std::vector<std::vector<bool> > const& used) {

  unsigned int neighbours = 0;

  // Loop over all directly neighbouring hits (not itself)
  for (int x = -1; x <= 1; ++x) {
    for (int y = -1; y <= 1; ++y) {
      if (x == 0 and y == 0)
	continue;
      if (used[xbin+x][ybin+y])
	neighbours++;
    }
  }

  // Return the number of neighbours in the cluster of a particular hit
  return neighbours;

}

TVector2 cluster::BlurredClusterAlg::Project3DPointOntoPlane(TVector3 const& point, int plane, int cryostat) {

  TVector2 wireTickPos = TVector2(-999., -999.);

  double pointPosition[3] = {point.X(), point.Y(), point.Z()};

  geo::TPCID tpcID = fGeom->FindTPCAtPosition(pointPosition);
  int tpc = 0;
  if (tpcID.isValid)
    tpc = tpcID.TPC;
  else
    return wireTickPos;

  // Construct wire ID for this point projected onto the plane
  geo::PlaneID planeID = geo::PlaneID(cryostat, tpc, plane);
  geo::WireID wireID = fGeom->NearestWireID(point, planeID);

  wireTickPos = TVector2(GlobalWire(wireID),
                         fDetProp->ConvertXToTicks(point.X(), planeID));

  return wireTickPos;

}

bool cluster::BlurredClusterAlg::RealHit(std::vector<std::vector<double> > const& image, int bin) {

  art::Ptr<recob::Hit> hit = ConvertBinToRecobHit(image, bin);

  return hit.isNonnull();

}

void cluster::BlurredClusterAlg::SaveImage(TH2F* image, std::vector<art::PtrVector<recob::Hit> > const& allClusters, int pad, int tpc, int plane) {

  // Make a vector of clusters
  std::vector<std::vector<int> > allClusterBins;

  for (std::vector<art::PtrVector<recob::Hit> >::const_iterator clusterIt = allClusters.begin(); clusterIt != allClusters.end(); clusterIt++) {
    art::PtrVector<recob::Hit> cluster = *clusterIt;

    if (!cluster.size())
      continue;

    std::vector<int> clusterBins;

    for (art::PtrVector<recob::Hit>::iterator hitIt = cluster.begin(); hitIt != cluster.end(); hitIt++) {
      art::Ptr<recob::Hit> hit = *hitIt;
      unsigned int wire = GlobalWire(hit->WireID());
      float tick = hit->PeakTime();
      int bin = image->GetBin((wire-fLowerWire)+1,(tick-fLowerTick)+1);
      if (cluster.size() < fMinSize)
        bin *= -1;

      clusterBins.push_back(bin);
    }

    allClusterBins.push_back(clusterBins);
  }

  SaveImage(image, allClusterBins, pad, tpc, plane);
}

void cluster::BlurredClusterAlg::SaveImage(TH2F* image, int pad, int tpc, int plane) {
  std::vector<std::vector<int> > allClusterBins;
  SaveImage(image, allClusterBins, pad, tpc, plane);
}

void cluster::BlurredClusterAlg::SaveImage(TH2F* image, std::vector<std::vector<int> > const& allClusterBins, int pad, int tpc, int plane) {

  fDebugCanvas->cd(pad);
  std::string stage;

  switch (pad) {
    case 1:
      stage = "Stage 1: Unblurred";
      break;
    case 2:
      stage = "Stage 2: Blurred";
      break;
    case 3:
      stage = "Stage 3: Blurred with clusters overlaid";
      break;
    case 4:
      stage = "Stage 4: Output clusters";
      break;
    default:
      stage = "Unknown stage";
      break;
  }

  std::stringstream title;
  title << stage << " -- TPC " << tpc << ", Plane " << plane;// << " (Event " << fEvent << ")";

  image->SetName(title.str().c_str());
  image->SetTitle(title.str().c_str());
  image->DrawCopy("colz");

  // Draw the clustered hits on the histograms
  int clusterNum = 2;
  for (std::vector<std::vector<int> >::const_iterator it = allClusterBins.begin(); it != allClusterBins.end(); it++, clusterNum++) {
    std::vector<int> bins = *it;
    TMarker mark(0, 0, 20);
    mark.SetMarkerColor(clusterNum);
    mark.SetMarkerSize(0.1);

    for (std::vector<int>::iterator binIt = bins.begin(); binIt != bins.end(); binIt++) {
      int bin = *binIt;
      int wire, tick, z;

      // Hit from a cluster that we aren't going to save
      if (bin < 0) {
        bin *= -1;
        mark.SetMarkerStyle(24);
      }

      image->GetBinXYZ(bin,wire,tick,z);
      mark.DrawMarker(wire+fLowerWire-1, tick+fLowerTick-1);
      mark.SetMarkerStyle(20);
    }
  }

  if (pad == 4) {
    fDebugCanvas->Print(fDebugPDFName.c_str());
    fDebugCanvas->Clear("D");
  }

}
