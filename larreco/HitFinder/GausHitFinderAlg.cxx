#ifndef GAUS_HIT_FINDER_ALG_CXX
#define GAUS_HIT_FINDER_ALG_CXX


#include "GausHitFinderAlg.h"

namespace hit {

void GausHitFinderAlg::beginJob() {
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  // ======================================
  // === Hit Information for Histograms ===
  fFirstChi2  = tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 5000);
  fChi2           = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 5000);
}

GausHitFinderAlg::GausHitFinderAlg(fhicl::ParameterSet const& p ) {

  bool const doHitFiltering = p.get<bool>("FilterHits", false);
  if (doHitFiltering) {
    if (fHitFilterAlg) { // create a new algorithm instance
      fHitFilterAlg->reconfigure(p.get<fhicl::ParameterSet>("HitFilterAlg"));
    }
    else { // reconfigure the existing instance
      fHitFilterAlg = std::make_unique<HitFilterAlg>
                      (p.get<fhicl::ParameterSet>("HitFilterAlg"));
    }
  }

  FillOutHitParameterVector(p.get< std::vector<double> >("MinSig"),         fMinSig);
  FillOutHitParameterVector(p.get< std::vector<double> >("InitWidth"),      fInitWidth);
  FillOutHitParameterVector(p.get< std::vector<double> >("MinWidth"),       fMinWidth);
  FillOutHitParameterVector(p.get< std::vector<double> >("AreaNorms"),      fAreaNorms);

  fLongMaxHits      = p.get< std::vector<int>>("LongMaxHits",    std::vector<int>() = {25, 25, 25});
  fLongPulseWidth   = p.get< std::vector<int>>("LongPulseWidth", std::vector<int>() = {16, 16, 16});
  fMaxMultiHit      = p.get< int             >("MaxMultiHit");
  fAreaMethod       = p.get< int             >("AreaMethod");
  fTryNplus1Fits    = p.get< bool            >("TryNplus1Fits");
  fChi2NDFRetry     = p.get< double          >("Chi2NDFRetry");
  fChi2NDF          = p.get< double          >("Chi2NDF");
  fNumBinsToAverage = p.get< size_t          >("NumBinsToAverage", 0);
}


void GausHitFinderAlg::RunOnWire(art::Ptr<recob::Wire> const& wire,
                                 std::vector<recob::Hit>& hit_vector,
                                 geo::Geometry const& geom)
{

  // We'll use the view as well...
  geo::View_t view = wire->View();


  // ----------------------------------------------------------
  // -- Setting the appropriate signal widths and thresholds --
  // --    for the right plane.      --
  // ----------------------------------------------------------

  threshold = fMinSig.at(view);
  fitWidth  = fInitWidth.at(view);
  minWidth  = fMinWidth.at(view);

  // #################################################
  // ### Set up to loop over ROI's for this wire   ###
  // #################################################
  const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();
  for (const auto& range : signalROI.get_ranges())
  {

    // Find the hits for this ROI:
    RunOnROI(range, hit_vector, wire, geom);

  }
}

void GausHitFinderAlg::RunOnROI(const lar::sparse_vector<float>::datarange_t & range,
                                std::vector<recob::Hit> & hit_vector,
                                art::Ptr<recob::Wire> wire,
                                geo::Geometry const& geom)
{

  // --- Setting Channel Number and Signal type ---
  raw::ChannelID_t channel =  wire->Channel();


  // get the WireID for this hit
  std::vector<geo::WireID> wids = geom.ChannelToWire(channel);
  // for now, just take the first option returned from ChannelToWire
  geo::WireID wid  = wids[0];
  // We'll use the view as well...
  geo::View_t view = wire->View();

  // #################################################
  // ### Getting a vector of signals for this wire ###
  // #################################################

  const std::vector<float>& signal = range.data();

  // ##########################################################
  // ### Making an iterator for the time ticks of this wire ###
  // ##########################################################
  std::vector<float>::const_iterator timeIter;        // iterator for time bins

  // ROI start time
  raw::TDCtick_t roiFirstBinTick = range.begin_index();

  MergedTimeWidVec mergedVec;
  float            roiThreshold(threshold);


  //#################################################
  //###    Set the charge determination method    ###
  //### Default is to compute the normalized area ###
  //#################################################
  std::function<double (double, double, double, double, int, int)> chargeFunc = [](double peakMean, double peakAmp, double peakWidth, double areaNorm, int low, int hi) {return std::sqrt(2 * TMath::Pi()) * peakAmp * peakWidth / areaNorm;};

  //##############################################
  //### Alternative is to integrate over pulse ###
  //##############################################
  if (fAreaMethod == 0)
    chargeFunc = [](double peakMean, double peakAmp, double peakWidth, double areaNorm, int low, int hi)
  {
    double charge(0);
    for (int sigPos = low; sigPos < hi; sigPos++)
      charge += peakAmp * TMath::Gaus(sigPos, peakMean, peakWidth);
    return charge;
  };



  // ###########################################################
  // ### If option set do bin averaging before finding peaks ###
  // ###########################################################

  if (fNumBinsToAverage > 1)
  {
    std::vector<float> timeAve;

    doBinAverage(signal, timeAve, fNumBinsToAverage);

    // ###################################################################
    // ### Search current averaged ROI for candidate peaks and widths  ###
    // ###################################################################

    TimeValsVec timeValsVec;
    findCandidatePeaks(timeAve.begin(), timeAve.end(), timeValsVec, roiThreshold, 0);

    // ####################################################
    // ### If no startTime hit was found skip this wire ###
    // ####################################################
    if (timeValsVec.empty()) return;

    // #############################################################
    // ### Merge potentially overlapping peaks and do multi fit  ###
    // #############################################################

    mergeCandidatePeaks(timeAve, timeValsVec, mergedVec);
  }

  // ###########################################################
  // ### Otherwise, operate directonly on signal vector      ###
  // ###########################################################

  else
  {
    // ##########################################################
    // ### Search current ROI for candidate peaks and widths  ###
    // ##########################################################

    TimeValsVec timeValsVec;
    findCandidatePeaks(signal.begin(), signal.end(), timeValsVec, roiThreshold, 0);

    // ####################################################
    // ### If no startTime hit was found skip this wire ###
    // ####################################################
    if (timeValsVec.empty()) return;

    // #############################################################
    // ### Merge potentially overlapping peaks and do multi fit  ###
    // #############################################################

    mergeCandidatePeaks(signal, timeValsVec, mergedVec);
  }

  // #######################################################
  // ### Lets loop over the pulses we found on this wire ###
  // #######################################################

  for (auto& mergedCands : mergedVec)
  {
    int             startT   = std::get<0>(mergedCands);
    int             endT     = std::get<1>(mergedCands);
    PeakTimeWidVec& peakVals = std::get<2>(mergedCands);

    // ### Putting in a protection in case things went wrong ###
    // ### In the end, this primarily catches the case where ###
    // ### a fake pulse is at the start of the ROI           ###
    if (endT - startT < 5) continue;

    // #######################################################
    // ### Clearing the parameter vector for the new pulse ###
    // #######################################################

    // === Setting the number of Gaussians to try ===
    int nGausForFit = peakVals.size();

    // ##################################################
    // ### Calling the function for fitting Gaussians ###
    // ##################################################
    double       chi2PerNDF(0.);
    int          NDF(0);
    ParameterVec paramVec;

    // #######################################################
    // ### If # requested Gaussians is too large then punt ###
    // #######################################################
    if (peakVals.size() <= fMaxMultiHit)
    {
      FitGaussians(signal, peakVals, startT, endT, 1.0, paramVec, chi2PerNDF, NDF);

      // If the chi2 is infinite then there is a real problem so we bail
      if (!(chi2PerNDF < std::numeric_limits<double>::infinity())) continue;

      fFirstChi2->Fill(chi2PerNDF);

      // #######################################################
      // ### Clearing the parameter vector for the new pulse ###
      // #######################################################
      double       chi2PerNDF2(0.);
      int          NDF2(0);
      ParameterVec paramVec2;

      // #####################################################
      // ### Trying extra gaussians for an initial bad fit ###
      // #####################################################
      if ( (chi2PerNDF > (2 * fChi2NDFRetry) && fTryNplus1Fits && nGausForFit == 1) ||
           (chi2PerNDF > (fChi2NDFRetry)   && fTryNplus1Fits && nGausForFit >  1))
      {
        // ############################################################
        // ### Modify input parameters for re-fitting n+1 Gaussians ###
        // ############################################################
        int newPeakTime = peakVals[0].first + 5 * nGausForFit;

        // We need to make sure we are not out of range and new peak amplitude is non-negative
        if (newPeakTime < endT - 1 && signal[newPeakTime] > 0.)
        {
          peakVals.emplace_back(newPeakTime, 2. * peakVals[0].second);

          // #########################################################
          // ### Calling the function for re-fitting n+1 Gaussians ###
          // #########################################################
          FitGaussians(signal, peakVals, startT, endT, 0.5, paramVec2, chi2PerNDF2, NDF2);

          // #########################################################
          // ### Getting the appropriate parameter into the vector ###
          // #########################################################
          if (chi2PerNDF2 < chi2PerNDF)
          {
            nGausForFit = peakVals.size();
            chi2PerNDF  = chi2PerNDF2;
            NDF         = NDF2;
            paramVec    = paramVec2;
          }
        }
      }
    }

    // #######################################################
    // ### If too large then force alternate solution      ###
    // ### - Make n hits from pulse train where n will     ###
    // ###   depend on the fhicl parameter fLongPulseWidth ###
    // ### Also do this if chi^2 is too large              ###
    // #######################################################
    if (peakVals.size() > fMaxMultiHit || chi2PerNDF > fChi2NDF)
    {
      int longPulseWidth = fLongPulseWidth.at(view);
      int nHitsThisPulse = (endT - startT) / longPulseWidth;

      if (nHitsThisPulse > fLongMaxHits.at(view))
      {
        nHitsThisPulse = fLongMaxHits.at(view);
        longPulseWidth = (endT - startT) / nHitsThisPulse;
      }

      if (nHitsThisPulse * longPulseWidth < endT - startT) nHitsThisPulse++;

      int firstTick = startT;
      int lastTick  = firstTick + std::min(endT, longPulseWidth);

      paramVec.clear();
      nGausForFit = nHitsThisPulse;
      NDF         = 1.;
      chi2PerNDF  =  chi2PerNDF > fChi2NDF ? chi2PerNDF : -1.;

      for (int hitIdx = 0; hitIdx < nHitsThisPulse; hitIdx++)
      {
        // This hit parameters
        double sumADC    = std::accumulate(signal.begin() + firstTick, signal.begin() + lastTick, 0.);
        double peakSigma = (lastTick - firstTick) / 3.;  // Set the width...
        double peakAmp   = 0.3989 * sumADC / peakSigma;  // Use gaussian formulation
        double peakMean  = (firstTick + lastTick) / 2.;

        // Store hit params
        paramVec.emplace_back(peakAmp,   0.1 * peakAmp);
        paramVec.emplace_back(peakMean,  0.1 * peakMean);
        paramVec.emplace_back(peakSigma, 0.1 * peakSigma);

        // set for next loop
        firstTick = lastTick;
        lastTick  = std::min(lastTick  + longPulseWidth, endT);
      }
    }

    // #######################################################
    // ### Loop through returned peaks and make recob hits ###
    // #######################################################

    int numHits(0);

    for (int hitIdx = 0; hitIdx < 3 * nGausForFit; hitIdx += 3)
    {
      // Extract values for this hit
      double peakAmp   = paramVec[hitIdx    ].first;
      double peakMean  = paramVec[hitIdx + 1].first;
      double peakWidth = paramVec[hitIdx + 2].first;

      // Extract errors
      double peakAmpErr   = paramVec[hitIdx    ].second;
      double peakMeanErr  = paramVec[hitIdx + 1].second;
      double peakWidthErr = paramVec[hitIdx + 2].second;

      // ### Charge ###
      double charge    = chargeFunc(peakMean, peakAmp, peakWidth, fAreaNorms[view], startT, endT);;
      double chargeErr = std::sqrt(TMath::Pi()) * (peakAmpErr * peakWidthErr + peakWidthErr * peakAmpErr);

      // ### limits for getting sums
      std::vector<float>::const_iterator sumStartItr = signal.begin() + startT;
      std::vector<float>::const_iterator sumEndItr   = signal.begin() + endT;

      // ### Sum of ADC counts
      double sumADC = std::accumulate(sumStartItr, sumEndItr, 0.);

      // ok, now create the hit
      recob::HitCreator hitcreator(*wire,                            // wire reference
                                   wid,                              // wire ID
                                   startT + roiFirstBinTick,         // start_tick TODO check
                                   endT + roiFirstBinTick,           // end_tick TODO check
                                   peakWidth,                        // rms
                                   peakMean + roiFirstBinTick,       // peak_time
                                   peakMeanErr,                      // sigma_peak_time
                                   peakAmp,                          // peak_amplitude
                                   peakAmpErr,                       // sigma_peak_amplitude
                                   charge,                           // hit_integral
                                   chargeErr,                        // hit_sigma_integral
                                   sumADC,                           // summedADC FIXME
                                   nGausForFit,                      // multiplicity
                                   numHits,                          // local_index TODO check that the order is correct
                                   chi2PerNDF,                       // goodness_of_fit
                                   NDF                               // dof
                                  );

      const recob::Hit hit(hitcreator.move());

      if (!fHitFilterAlg || fHitFilterAlg->IsGoodHit(hit)) {
        hit_vector.emplace_back(std::move(hit));
        numHits++;
      }
    } // <---End loop over gaussians

    fChi2->Fill(chi2PerNDF);

  }//<---End loop over merged candidate hits


  return;
}


//-------------------------------------------------
//-------------------------------------------------
void GausHitFinderAlg::FillOutHitParameterVector(
  const std::vector<double>& input,
  std::vector<double>&       output)
{
  if (input.size() == 0)
    throw std::runtime_error("GausHitFinderAlg::FillOutHitParameterVector ERROR! Input config vector has zero size.");

  art::ServiceHandle<geo::Geometry> geom;
  const unsigned int N_PLANES = geom->Nplanes();

  if (input.size() == 1)
    output.resize(N_PLANES, input[0]);
  else if (input.size() == N_PLANES)
    output = input;
  else
    throw std::runtime_error("GausHitFinderAlg::FillOutHitParameterVector ERROR! Input config vector size !=1 and !=N_PLANES.");

}



// --------------------------------------------------------------------------------------------
// Initial finding of candidate peaks
// --------------------------------------------------------------------------------------------
void hit::GausHitFinderAlg::findCandidatePeaks(
  std::vector<float>::const_iterator    startItr,
  std::vector<float>::const_iterator    stopItr,
  std::vector<std::tuple<int, int, int>>& timeValsVec,
  float&                                roiThreshold,
  int                                   firstTick) const
{
  // Need a minimum number of ticks to do any work here
  if (std::distance(startItr, stopItr) > 4)
  {
    // Find the highest peak in the range given
    auto maxItr = std::max_element(startItr, stopItr);

    float maxValue = *maxItr;
    int   maxTime  = std::distance(startItr, maxItr);

    if (maxValue > roiThreshold)
    {
      // backwards to find first bin for this candidate hit
      auto firstItr = std::distance(startItr, maxItr) > 2 ? maxItr - 1 : startItr;

      while (firstItr != startItr)
      {
        // Check for pathology where waveform goes too negative
        if (*firstItr < -roiThreshold) break;

        // Check both sides of firstItr and look for min/inflection point
        if (*firstItr < * (firstItr + 1) && *firstItr <= *(firstItr - 1)) break;

        firstItr--;
      }

      int firstTime = std::distance(startItr, firstItr);

      // Recursive call to find all candidate hits earlier than this peak
      findCandidatePeaks(startItr, firstItr + 1, timeValsVec, roiThreshold, firstTick);

      // forwards to find last bin for this candidate hit
      auto lastItr = std::distance(maxItr, stopItr) > 2 ? maxItr + 1 : stopItr - 1;

      while (lastItr != stopItr - 1)
      {
        // Check for pathology where waveform goes too negative
        if (*lastItr < -roiThreshold) break;

        // Check both sides of firstItr and look for min/inflection point
        if (*lastItr <= *(lastItr + 1) && *lastItr < * (lastItr - 1)) break;

        lastItr++;
      }

      int lastTime = std::distance(startItr, lastItr);

      // Now save this candidate's start and max time info
      timeValsVec.push_back(std::make_tuple(firstTick + firstTime, firstTick + maxTime, firstTick + lastTime));

      // Recursive call to find all candidate hits later than this peak
      findCandidatePeaks(lastItr + 1, stopItr, timeValsVec, roiThreshold, firstTick + std::distance(startItr, lastItr + 1));
    }
  }

  return;
} // findCandidatePeaks



// --------------------------------------------------------------------------------------------
// Merging of nearby candidate peaks
// --------------------------------------------------------------------------------------------

void hit::GausHitFinderAlg::mergeCandidatePeaks(
  const std::vector<float>& signalVec,
  TimeValsVec& timeValsVec,
  MergedTimeWidVec& mergedVec) const
{
  // ################################################################
  // ### Lets loop over the candidate pulses we found in this ROI ###
  // ################################################################
  auto timeValsVecItr = timeValsVec.begin();

  while (timeValsVecItr != timeValsVec.end())
  {
    PeakTimeWidVec peakVals;

    // Setting the start, peak, and end time of the pulse
    auto& timeVal = *timeValsVecItr++;

    int startT = std::get<0>(timeVal);
    int maxT   = std::get<1>(timeVal);
    int endT   = std::get<2>(timeVal);
    int widT   = std::max(2, (endT - startT) / 6);

    peakVals.emplace_back(maxT, widT);

    // See if we want to merge pulses together
    // First check if we have more than one pulse on the wire
    bool checkNextHit = timeValsVecItr != timeValsVec.end();

    // Loop until no more merged pulses (or candidates in this ROI)
    while (checkNextHit)
    {
      // If the start time of the next pulse is the end time of the current pulse then merge
      // Alternatively, if the start time of the next pulses is one tick away then
      // merge if the intervening signal is above 0.5 * threshold
      int nextStartT = std::get<0>(*timeValsVecItr);

      if ((nextStartT == endT) || (nextStartT - endT  < 2 && signalVec[endT + 1] > threshold / 2))
      {
        timeVal = *timeValsVecItr++;
        maxT    = std::get<1>(timeVal);
        endT    = std::get<2>(timeVal);
        widT    = std::max(2, (endT - nextStartT) / 6);

        peakVals.emplace_back(maxT, widT);

        checkNextHit = timeValsVecItr != timeValsVec.end();
      }//<---Checking adjacent pulses
      else checkNextHit = false;

    }//<---End checking if there is more than one pulse on the wire

    // Add these to our merged vector
    mergedVec.emplace_back(startT, endT, peakVals);
  }

  return;
}



// --------------------------------------------------------------------------------------------
// Fit Gaussians
// --------------------------------------------------------------------------------------------
void hit::GausHitFinderAlg::FitGaussians(const std::vector<float>& SignalVector,
    const PeakTimeWidVec&     PeakVals,
    int                       StartTime,
    int                       EndTime,
    double                    ampScaleFctr,
    ParameterVec&             paramVec,
    double&                   chi2PerNDF,
    int&                      NDF)
{
  int size = EndTime - StartTime;
  // #############################################
  // ### If size < 0 then set the size to zero ###
  // #############################################
  if (EndTime - StartTime < 0) {size = 0;}

  // --- TH1D HitSignal ---
  TH1F hitSignal("hitSignal", "", std::max(size, 1), StartTime, EndTime);
  hitSignal.Sumw2();

  // #############################
  // ### Filling the histogram ###
  // #############################
  for (int aa = StartTime; aa < EndTime; aa++)
  {
    hitSignal.Fill(aa, SignalVector[aa]);

    if (EndTime > 10000) {break;} // FIXME why?
  }//<---End aa loop

  // ############################################
  // ### Building TFormula for basic Gaussian ###
  // ############################################
  std::string eqn = "gaus(0)";  // string for equation for gaus fit
  std::stringstream numConv;

  for (size_t i = 3; i < PeakVals.size() * 3; i += 3)
  {
    eqn.append("+gaus(");
    numConv.str("");
    numConv << i;
    eqn.append(numConv.str());
    eqn.append(")");
  }

  // ---------------------------------
  // --- TF1 function for GausHit  ---
  // ---------------------------------
  TF1 Gaus("Gaus", eqn.c_str(), 0, std::max(size, 1));

  // ### Setting the parameters for the Gaussian Fit ###
  int parIdx(0);
  for (auto& peakVal : PeakVals)
  {
    double peakMean   = peakVal.first;
    double peakWidth  = peakVal.second;
    double amplitude  = ampScaleFctr * SignalVector[peakMean];
    double meanLowLim = std::max(peakMean - 2.*peakWidth, double(StartTime));
    double meanHiLim  = std::min(peakMean + 2.*peakWidth, double(EndTime));

    Gaus.SetParameter(  parIdx, amplitude);
    Gaus.SetParameter(1 + parIdx, peakMean);
    Gaus.SetParameter(2 + parIdx, peakWidth);
    Gaus.SetParLimits(  parIdx, 0.0,        1.5 * amplitude);
    Gaus.SetParLimits(1 + parIdx, meanLowLim, meanHiLim);
    Gaus.SetParLimits(2 + parIdx, minWidth,   10.*peakWidth);

    parIdx += 3;
  }//<---End bb loop

  // ####################################################
  // ### PERFORMING THE TOTAL GAUSSIAN FIT OF THE HIT ###
  // ####################################################

  //hitSignal.Fit(&Gaus,"QNRWIB","", StartTime, EndTime);
  //hitSignal.Fit(&Gaus,"QNRWB","", StartTime, EndTime);
  //hitGraph.Fit(&Gaus,"QNB","",StartTime, EndTime);

  try
  { hitSignal.Fit(&Gaus, "QNRWB", "", StartTime, EndTime);}
  catch (...)
  {mf::LogWarning("GausHitFinderAlg") << "Fitter failed finding a hit";}

  // ##################################################
  // ### Getting the fitted parameters from the fit ###
  // ##################################################
  chi2PerNDF = (Gaus.GetChisquare() / Gaus.GetNDF());
  NDF        = Gaus.GetNDF();

  for (size_t ipar = 0; ipar < (3 * PeakVals.size()); ++ipar)
    paramVec.emplace_back(Gaus.GetParameter(ipar), Gaus.GetParError(ipar));

  Gaus.Delete();
  hitSignal.Delete();
}//<----End FitGaussians



void hit::GausHitFinderAlg::doBinAverage(const std::vector<float>& inputVec,
    std::vector<float>&       outputVec,
    size_t                    binsToAverage) const
{
  size_t halfBinsToAverage(binsToAverage / 2);

  float runningSum(0.);

  for (size_t idx = 0; idx < halfBinsToAverage; idx++) runningSum += inputVec[idx];

  outputVec.resize(inputVec.size());
  std::vector<float>::iterator outputVecItr = outputVec.begin();

  // First pass through to build the erosion vector
  for (std::vector<float>::const_iterator inputItr = inputVec.begin(); inputItr != inputVec.end(); inputItr++)
  {
    size_t startOffset = std::distance(inputVec.begin(), inputItr);
    size_t stopOffset  = std::distance(inputItr, inputVec.end());
    size_t count       = std::min(2 * halfBinsToAverage, std::min(startOffset + halfBinsToAverage + 1, halfBinsToAverage + stopOffset - 1));

    if (startOffset >= halfBinsToAverage) runningSum -= *(inputItr - halfBinsToAverage);
    if (stopOffset  >  halfBinsToAverage) runningSum += *(inputItr + halfBinsToAverage);

    *outputVecItr++ = runningSum / float(count);
  }

  return;
}


void hit::GausHitFinderAlg::reBin(const std::vector<float>& inputVec,
                                  std::vector<float>&       outputVec,
                                  size_t                    nBinsToCombine) const
{
  size_t nNewBins = inputVec.size() / nBinsToCombine;

  if (inputVec.size() % nBinsToCombine > 0) nNewBins++;

  outputVec.resize(nNewBins, 0.);

  size_t outputBin = 0;

  for (size_t inputIdx = 0; inputIdx < inputVec.size();)
  {
    outputVec[outputBin] += inputVec[inputIdx++];

    if (inputIdx % nBinsToCombine == 0) outputBin++;

    if (outputBin > outputVec.size())
    {
      std::cout << "***** DISASTER!!! ****** outputBin: " << outputBin << ", inputIdx = " << inputIdx << std::endl;
      break;
    }
  }

  return;
}

} // namespace hit

#endif // GAUS_HIT_FINDER_ALG_CXX
