#ifndef GAUSHITFINDER_H
#define GAUSHITFINDER_H

////////////////////////////////////////////////////////////////////////
//
// GaussHitFinder class
//
// jaasaadi@syr.edu
//
//  This algorithm is designed to find hits on wires after deconvolution.
// -----------------------------------
// This algorithm is based on the FFTHitFinder written by Brian Page,
// Michigan State University, for the ArgoNeuT experiment.
//
//
// The algorithm walks along the wire and looks for pulses above threshold
// The algorithm then attempts to fit n-gaussians to these pulses where n
// is set by the number of peaks found in the pulse
// If the Chi2/NDF returned is "bad" it attempts to fit n+1 gaussians to
// the pulse. If this is a better fit it then uses the parameters of the
// Gaussian fit to characterize the "hit" object
//
// To use this simply include the following in your producers:
// gaushit:     @local::microboone_gaushitfinder
// gaushit:	@local::argoneut_gaushitfinder
////////////////////////////////////////////////////////////////////////


// C/C++ standard library
#include <algorithm> // std::accumulate()
#include <vector>
#include <string>
#include <memory> // std::unique_ptr()
#include <utility> // std::move()
//#include <Rcpp.h> // nls()
#include <sstream>

#if defined(GCC_VERSION)
#undef GCC_VERSION
#endif

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"


// LArSoft Includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoBaseArt/HitCreator.h"
#include "HitFilterAlg.h"

#include <iostream>
#include <pthread.h>
#include <omp.h>
#include <time.h>

// ROOT Includes
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TThread.h"

#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "HFitInterface.h"
#include "Fit/Fitter.h"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

namespace hit{

  /******************************************************************************
   * The fit function for thread safe applications - Gaus_TS
  */
  class Gaus_TS : public ROOT::Math::IParamFunction {

    private:
      //---------------------------------------->
      //
      double DoEvalPar(double x, const double *par) const{
        double fitval = par[0];
        double arg = 0;

        for(int igaus = 0; igaus < npar/4; igaus++) {
          if(par[3+4*igaus] != 0){
            arg = (x - par[2+4*igaus])/par[3+4*igaus];
          }
          fitval += par[1+4*igaus]*TMath::Exp(-0.5*arg*arg);
        }

        return fitval;
      }

      //---------------------------------------->
      //
      std::vector<double> fp;
      int                 npar;
    int Counter;
    public:
      //---------------------------------------->
      //
      virtual ~Gaus_TS () {
	//	std::cout<<"Deleting: "<<Counter<<std::endl;
      }

      Gaus_TS(int p) : npar(p) {
	fp.reserve(NPar());
        //if (npar > 4) std::cout << "Gaus_TS constructor: Reserved fp with " << NPar() << " parameters." << std::endl;

	static int sCounter = 0;
	Counter = sCounter++;
	//	std::cout<<"Created: "<<Counter<<std::endl;
      }

      //---------------------------------------->
      //
      void SetParameters (const double *p) {
        fp.clear();
        fp.insert(fp.begin(), p, p+NPar());
      }

      //---------------------------------------->
      //
      const double *Parameters() const {
        return fp.data();
      }

      //---------------------------------------->
      //
      ROOT::Math::IGenFunction *Clone() const {
      //  if (npar>4) std::cout << "Gaus_TS_Clone(): npar, fp.size() " << npar << ", " << fp.size() << std::endl;
        Gaus_TS *gs = new Gaus_TS(npar);
        gs->SetParameters(fp.data());
        return gs;
      }

      //---------------------------------------->
      //
      unsigned int NPar() const {
        const int const_npar(npar);
        return const_npar;
      }

  };//==>end of class Gaus_TS


  class GausHitFinder : public art::EDProducer {

  public:

    explicit GausHitFinder(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt) override;
    void beginJob() override;
    void endJob() override;
    void reconfigure(fhicl::ParameterSet const& p) override;


  private:

    using TimeValsVec      = std::vector<std::tuple<int,int,int>>;
    using PeakTimeWidVec   = std::vector<std::pair<int,int>>;
    using MergedTimeWidVec = std::vector<std::tuple<int,int,PeakTimeWidVec>>;

    void findCandidatePeaks(std::vector<float>::const_iterator startItr,
                            std::vector<float>::const_iterator stopItr,
                            TimeValsVec&                       timeValsVec,
                            float&                             roiThreshold,
                            int                                firstTick) const;

    void mergeCandidatePeaks(const std::vector<float>&, TimeValsVec&, MergedTimeWidVec&) const;

    // ### This function will fit N-Gaussians to at TH1D where N is set ###
    // ###            by the number of peaks found in the pulse         ###

    using ParameterVec = std::vector<std::pair<double,double>>;  //< parameter/error vec

    void FillOutHitParameterVector(const std::vector<double>& input,
				   std::vector<double>& output);

    void doBinAverage(const std::vector<float>& inputVec,
                      std::vector<float>&       outputVec,
                      size_t                    binsToAverage) const;

    void reBin(const std::vector<float>& inputVec,
               std::vector<float>&       outputVec,
               size_t                    nBinsToCombine) const;

  public:

    std::map <unsigned int, double > fchi2pndf ;         /////// maps of results from fitter, indexed by threadID
    std::map <unsigned int, int> fndf;
    std::map < unsigned int, std::vector<std::pair<double,double> > > fparamvec;
    std::map < unsigned int, bool > fthreaddone;

    // FitGaussians() must now be Public if it's to be accessed by TThread external function call.
    void FitGaussians(const std::vector<float>& SignalVector,
                      const PeakTimeWidVec&     PeakVals,
                      int                       StartTime,
                      int                       EndTime,
                      double                    ampScaleFctr,
                      ParameterVec&             paramVec,
                      double&                   chi2PerNDF,
                      int&                      NDF,
		      const int& t
		      );

    void FitGaussians_gsl(const std::vector<float>& SignalVector,
                      const PeakTimeWidVec&     PeakVals,
                      int                       StartTime,
                      int                       EndTime,
                      double                    ampScaleFctr,
                      ParameterVec&             paramVec,
                      double&                   chi2PerNDF,
                      int&                      NDF,
		      const int& t
		      );

  private:
    double              threshold           = 0.;  // minimum signal size for id'ing a hit
    double              fitWidth            = 0.;  // hit fit width initial value
    double              minWidth		    = 0 ;  // hit minimum width
    std::string         fCalDataModuleLabel;

    std::vector<double> fMinSig;                   ///<signal height threshold
    std::vector<double> fInitWidth;                ///<Initial width for fit
    std::vector<double> fMinWidth;                 ///<Minimum hit width

    size_t              fMaxMultiHit;              ///<maximum hits for multi fit
    int                 fAreaMethod;               ///<Type of area calculation
    std::vector<double> fAreaNorms;                ///<factors for converting area to same units as peak height
    int		            fTryNplus1Fits;            ///<whether we will (0) or won't (1) try n+1 fits
    double	            fChi2NDFRetry;             ///<Value at which a second n+1 Fit will be tried
    double	            fChi2NDF;                  ///maximum Chisquared / NDF allowed for a hit to be saved
    size_t              fNumBinsToAverage;         ///< If bin averaging for peak finding, number bins to average

    std::unique_ptr<HitFilterAlg> fHitFilterAlg;   ///algorithm used to filter out noise hits

    TH1F* fFirstChi2;
    TH1F* fChi2;
    double              ftoolongusec = 3.0E6;  /// to kill long-running fit threads
    //    std::map <unsigned int, TH1F> fhitSignal;
    //    std::map <unsigned int, TF1* > fGaus;

    int histCnt;

  protected:


  }; // class GausHitFinder

  //---------------------------------------->
  //
  typedef struct{
    std::vector<float>                      SigVec;
    std::vector<std::pair<int, int>>        PeakV;
    int                                     StartT;
    int                                     EndT;
    double                                  AmpScaleF;
    std::vector<std::pair<double, double>>  ParamV;
    double                                  Chi2PerNDF;
    int                                     NDF;
  } args_t;
  //---------------------------------------->
  //

  /******************************************************************************
   * This function encapsulates the information transfered into the thread
   ** EC: appropriated 10-May-2016 from online ROOT chatter, Andreas Zoglauer (andreas@megalibtoolkit.com)
   */
  class ThreadCaller {

    private:
      //---------------------------------------->
      // store the calling class for retrieval
      GausHitFinder *fFitter;
      //---------------------------------------->
      //ID of the worker threadID
      unsigned int fThreadID;
      //extra arguments
      std::shared_ptr<args_t> fArgs;

    public:
      //---------------------------------------->
      // Standard constructor
      ThreadCaller(GausHitFinder *M, unsigned int ThreadID, std::shared_ptr<args_t> args){
        fFitter   = M;
        fThreadID = ThreadID;
        fArgs     = std::move(args);
      }

      //---------------------------------------->
      // Return the calling class
      GausHitFinder *GetThreadCaller(){
        return fFitter;
      }

      //---------------------------------------->
      // Return the thread ID
      unsigned int GetThreadID(){
        return fThreadID;
      }

      //---------------------------------------->
      //
      std::shared_ptr<args_t> GetArgs(){
        return fArgs;
      }

  };//==>end of class ThreadCaller



//-------------------------------------------------
//-------------------------------------------------
GausHitFinder::GausHitFinder(fhicl::ParameterSet const& pset)
{
    this->reconfigure(pset);

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);

} // GausHitFinder::GausHitFinder()


//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::FillOutHitParameterVector(const std::vector<double>& input,
                                              std::vector<double>&       output)
{
    if(input.size()==0)
        throw std::runtime_error("GausHitFinder::FillOutHitParameterVector ERROR! Input config vector has zero size.");

    art::ServiceHandle<geo::Geometry> geom;
    const unsigned int N_PLANES = geom->Nplanes();

    if(input.size()==1)
        output.resize(N_PLANES,input[0]);
    else if(input.size()==N_PLANES)
        output = input;
    else
        throw std::runtime_error("GausHitFinder::FillOutHitParameterVector ERROR! Input config vector size !=1 and !=N_PLANES.");

}


//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::reconfigure(fhicl::ParameterSet const& p)
{
    // Implementation of optional member function here.
    fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");


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

    FillOutHitParameterVector(p.get< std::vector<double> >("MinSig"),fMinSig);
    FillOutHitParameterVector(p.get< std::vector<double> >("InitWidth"),fInitWidth);
    FillOutHitParameterVector(p.get< std::vector<double> >("MinWidth"),fMinWidth);
    FillOutHitParameterVector(p.get< std::vector<double> >("AreaNorms"),fAreaNorms);

    fMaxMultiHit      = p.get< int          >("MaxMultiHit");
    fAreaMethod       = p.get< int          >("AreaMethod");
    fTryNplus1Fits    = p.get< int          >("TryNplus1Fits");
    fChi2NDFRetry     = p.get< double       >("Chi2NDFRetry");
    fChi2NDF          = p.get< double       >("Chi2NDF");
    fNumBinsToAverage = p.get< size_t       >("NumBinsToAverage", 0);
}

//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::beginJob()
{
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;


    // ======================================
    // === Hit Information for Histograms ===
    fFirstChi2	= tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 5000);
    fChi2	        = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 5000);
}

//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::endJob()
{

}

//  This algorithm uses the fact that deconvolved signals are very smooth
//  and looks for hits as areas between local minima that have signal above
//  threshold.
//-------------------------------------------------
void GausHitFinder::produce(art::Event& evt)
{
    //==================================================================================================
  clock_t str1, stp1;
    str1 = clock();

    TH1::AddDirectory(kFALSE);

    // Instantiate and Reset a stop watch
    //TStopwatch StopWatch;
    //StopWatch.Reset();

    // ################################
    // ### Calling Geometry service ###
    // ################################
    art::ServiceHandle<geo::Geometry> geom;

    // ###############################################
    // ### Making a ptr vector to put on the event ###
    // ###############################################
    // this contains the hit collection
    // and its associations to wires and raw digits
    recob::HitCollectionCreator hcol(*this, evt);

    // ##########################################
    // ### Reading in the Wire List object(s) ###
    // ##########################################
    art::Handle< std::vector<recob::Wire> > wireVecHandle;
    evt.getByLabel(fCalDataModuleLabel,wireVecHandle);

    // #################################################################
    // ### Reading in the RawDigit associated with these wires, too  ###
    // #################################################################
    art::FindOneP<raw::RawDigit> RawDigits
        (wireVecHandle, evt, fCalDataModuleLabel);

    // Channel Number
    raw::ChannelID_t channel = raw::InvalidChannelID;

    std::map < uint16_t, std::vector<recob::Hit>> hitsthreads; //hitsthreads.reserve((int)wireVecHandle->size());
    std::map < uint16_t, std:: vector <art::Ptr<recob::Wire>> > wiresthreads;// wiresthreads.reserve((int)wireVecHandle->size());
    std::map < uint16_t, std:: vector <art::Ptr<raw::RawDigit>> > digitsthreads;// digitsthreads.reserve((int)wireVecHandle->size());
    std::vector <int> t_ID(wireVecHandle->size());
    std::iota (std::begin(t_ID), std::end(t_ID), 0); // fill it with 0,1,2, .... 8255

    //##############################
    //### Looping over the wires ###
    //##############################

    //    std::cout << "\t Wire count: " << wireVecHandle->size() << std::endl;
    //#pragma omp parallel for schedule(static)
    std::vector<std::pair<double, double>> tmp1;
    tmp1.push_back(std::make_pair(0.0, 0.0));

    //initializing the keys for the maps. Must be done for thread safety later.

    for(int ii=0; ii<(int)wireVecHandle->size(); ii++){
      hitsthreads[ii] = std::vector<recob::Hit>();
      wiresthreads[ii] = std::vector<art::Ptr<recob::Wire>>();
      digitsthreads[ii] = std::vector<art::Ptr<raw::RawDigit>>();
      fthreaddone[ii] = false;
      fparamvec[ii] = tmp1;
      fndf[ii] = 0.0;
      fchi2pndf[ii] = 0.0;

      //      fhitSignal[ii] = TH1F();
      //      fGaus[ii] = NULL;
    }


    #pragma omp parallel
    {
    //iterating over all wires
    #pragma omp for schedule(dynamic)
    for(size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter++)
    {

      int tid = (int)wireIter;
      std:: vector<recob::Hit> onehitperthread;
      std:: vector<art::Ptr<recob::Wire>> onewireperthread;
      std:: vector<art::Ptr<raw::RawDigit>> onedigitperthread;
      //      std::cout << "[" << (int)omp_get_thread_num() << "] wire: " << wireIter << std::endl;

        // ####################################
        // ### Getting this particular wire ###
        // ####################################
        art::Ptr<recob::Wire>   wire(wireVecHandle, wireIter);
        art::Ptr<raw::RawDigit> rawdigits = RawDigits.at(wireIter);

        // --- Setting Channel Number and Signal type ---
        channel = wire->Channel();

        // get the WireID for this hit
        std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        // for now, just take the first option returned from ChannelToWire
        geo::WireID wid = wids[0];

        // ----------------------------------------------------------
        // -- Setting the appropriate signal widths and thresholds --
        // --    for the right plane.      --
        // ----------------------------------------------------------

        threshold = fMinSig.at(wire->View());
        fitWidth  = fInitWidth.at(wire->View());
        minWidth  = fMinWidth.at(wire->View());

//            if (wid.Plane == geo::kV)
//                roiThreshold = std::max(threshold,std::min(2.*threshold,*std::max_element(signal.begin(),signal.end())/3.));

        // #################################################
        // ### Set up to loop over ROI's for this wire   ###
        // #################################################
        const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();

        for(const auto& range : signalROI.get_ranges())
        {
            // #################################################
            // ### Getting a vector of signals for this wire ###
            // #################################################
            //std::vector<float> signal(wire->Signal());

            const std::vector<float>& signal = range.data();

            // ##########################################################
            // ### Making an iterator for the time ticks of this wire ###
            // ##########################################################
            std::vector<float>::const_iterator timeIter;  	    // iterator for time bins

            // ROI start time
            raw::TDCtick_t roiFirstBinTick = range.begin_index();

            MergedTimeWidVec mergedVec;
            float       roiThreshold(threshold);

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
                findCandidatePeaks(timeAve.begin(),timeAve.end(),timeValsVec,roiThreshold,0);

                // ####################################################
                // ### If no startTime hit was found skip this wire ###
                // ####################################################
                if (timeValsVec.empty()) continue;

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
                findCandidatePeaks(signal.begin(),signal.end(),timeValsVec,roiThreshold,0);

                // ####################################################
                // ### If no startTime hit was found skip this wire ###
                // ####################################################
                if (timeValsVec.empty()) continue;

                // #############################################################
                // ### Merge potentially overlapping peaks and do multi fit  ###
                // #############################################################

                mergeCandidatePeaks(signal, timeValsVec, mergedVec);
            }

            // #######################################################
            // ### Lets loop over the pulses we found on this wire ###
            // #######################################################

	    size_t ntt(0);
	    size_t ntt2(0);
            for(auto& mergedCands : mergedVec)
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
		int tid = (int)wireIter;
                if (peakVals.size() <= fMaxMultiHit)
                {
		  ntt++;
		  FitGaussians_gsl(signal, peakVals, startT, endT, 1.0, paramVec, chi2PerNDF, NDF, ntt);

		  //#####################################
		  //delete this later:
		  /*
		  std::ofstream outfile;
		  std::stringstream fname;
		  fname << "test_" << (int)wireIter << ".txt";
		  outfile.open(fname.str(), std::ios_base::app);
		  for(auto& aa : signal){
		    outfile << aa << "\t" << signal[aa] << std::endl;
		  }
		  */
		  //outfile << "wire id: ##########################" << std::endl;
		  //#####################################
      
      

		  //std::cout << "paramVec.size(), chi2PerNDF, NDF are: " << paramVec.size() << ", " << chi2PerNDF <<  ", " << NDF << std::endl;

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
                    if( (chi2PerNDF > (2*fChi2NDFRetry) && fTryNplus1Fits == 0 && nGausForFit == 1)||
                        (chi2PerNDF > (fChi2NDFRetry)   && fTryNplus1Fits == 0 && nGausForFit >  1))
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

			    ntt2++;
			    FitGaussians_gsl(signal, peakVals, startT, endT, 0.5, paramVec2, chi2PerNDF2, NDF2, tid);

                            // #########################################################
                            // ### Getting the appropriate parameter into the vector ###
                            // #########################################################
                            if (chi2PerNDF2 < chi2PerNDF && chi2PerNDF2 > 1E-12)
                            {
                                nGausForFit = peakVals.size();

				chi2PerNDF  = chi2PerNDF2;
				NDF         = NDF2;
				paramVec    = paramVec2;
                            } // end of chi2PerNDF2 if
			    ///thread2->Kill();
                        } // end if newPeakTime
                    } // end compound if
                } // end if peakVals.size()

                // ############################################
                // ### If too large then make one large hit ###
                // ### Also do this if chi^2 is too large   ###
                // ############################################
                if (peakVals.size() > fMaxMultiHit || chi2PerNDF > fChi2NDF)
                {
                    double sumADC    = std::accumulate(signal.begin() + startT, signal.begin() + endT,0.);
                    double peakAmp   = 1.5 * sumADC / (endT - startT);  // hedge between triangle and a box
                    double peakMean  = (startT + endT) / 2.;
                    double peakWidth = (endT - startT) / 4.;

                    nGausForFit =  1;
                    chi2PerNDF  =  chi2PerNDF > fChi2NDF ? chi2PerNDF : -1.;
                    NDF         =  1;

                    paramVec.clear();
                    paramVec.emplace_back(peakAmp,   0.1 * peakAmp);
                    paramVec.emplace_back(peakMean,  0.1 * peakMean);
                    paramVec.emplace_back(peakWidth, 0.1 * peakWidth);
                }

                // #######################################################
                // ### Loop through returned peaks and make recob hits ###
                // #######################################################

                int numHits(0);
		numHits = tid*10;
		if( !paramVec.size() ) {
		  std::cout << "paramVec.size() is 0..." << std::endl;
		  continue;
		}

                for(int hitIdx = 0; hitIdx < 3*nGausForFit; hitIdx+=3)
                {
                    // Extract values for this hit
                    double peakAmp   = paramVec[hitIdx    ].first;
                    double peakMean  = paramVec[hitIdx + 1].first;
                    double peakWidth = paramVec[hitIdx + 2].first;

                    // Selection cut
                    if (nGausForFit == 1 && peakAmp < threshold) continue;

                    // Extract errors
                    double peakAmpErr   = paramVec[hitIdx    ].second;
                    double peakMeanErr  = paramVec[hitIdx + 1].second;
                    double peakWidthErr = paramVec[hitIdx + 2].second;

                    // ### Charge ###
                    double totSig(0.);

                    // ######################################################
                    // ### Getting the total charge using the area method ###
                    // ######################################################
                    if(fAreaMethod)
                    {
                        totSig = std::sqrt(2*TMath::Pi())*peakAmp*peakWidth/fAreaNorms[(size_t)(wire->View())];
                    }//<---End Area Method

                    // ##################################
                    // ### Integral Method for charge ###
                    // ##################################
                    else
                    {
                        for(int sigPos = startT; sigPos < endT; sigPos++)
                            totSig += peakAmp * TMath::Gaus(sigPos,peakMean,peakWidth);
                    }

                    double charge(totSig);
                    double chargeErr = std::sqrt(TMath::Pi()) * (peakAmpErr*peakWidthErr + peakWidthErr*peakAmpErr);

                    // ### limits for getting sums
                    std::vector<float>::const_iterator sumStartItr = signal.begin() + startT;
                    std::vector<float>::const_iterator sumEndItr   = signal.begin() + endT;

                    // ### Sum of ADC counts
                    double sumADC = std::accumulate(sumStartItr, sumEndItr, 0.);

                    // ok, now create the hit
		    //		    std::cout << "\t create hit - OMP process: " << (int)omp_get_thread_num() << std::endl;
		    //		    std::cout << " nGausForFit, numHits, chi2pndf, signalType, ndf " << nGausForFit << "," << numHits << "," << chi2PerNDF << ", " << geom->SignalType(wire->Channel()) << "," << NDF << "." << std::endl;
		    //		    std::cout << " peakWidth,  peakMeanproi, peakMeanErr, peakAmp, charge, chargeErr, sumADC " << peakWidth << ", " <<peakMean+roiFirstBinTick<< ", " << peakMeanErr << ", " << peakAmp << ", " << charge << ", " << chargeErr << ", " << sumADC << ". " << std::endl;


                    recob::HitCreator hitcreator(*wire,                            // wire reference
                                        	 wid,                              // wire ID
                                        	 startT+roiFirstBinTick,           // start_tick TODO check
                                        	 endT+roiFirstBinTick,             // end_tick TODO check
                                        	 peakWidth,                        // rms
                                        	 peakMean+roiFirstBinTick,         // peak_time
                                        	 peakMeanErr,                      // sigma_peak_time
                                        	 peakAmp,                          // peak_amplitude
                                        	 peakAmpErr,                       // sigma_peak_amplitude
                                        	 charge,                           // hit_integral
                                        	 chargeErr,                        // hit_sigma_integral
                                        	 sumADC,                           // summedADC FIXME
                                        	 nGausForFit,                      // multiplicity
                                        	 numHits,                          // local_index TODO check that the order is correct
                                        	 chi2PerNDF,                       // goodness_of_fit
						 geom->SignalType(wire->Channel()),// EC added to avoid thread problem of invoking geom services down inside HitCreator().
                                        	 NDF                               // dof
                                        	 );

		    // const recob::Hit hit(hitcreator.move());

		    // hcol is global and can't be filled here in multi-threaded code.
		    // Instead, gather the arguments for later use.
		    /*
		    if (!fHitFilterAlg || fHitFilterAlg->IsGoodHit(hit)) {
                      hcol.emplace_back(std::move(hit), wire, rawdigits);
                      numHits++;
		    }
		    */

		    numHits++;
		    onehitperthread.emplace_back(hitcreator.move());
		    onewireperthread.emplace_back(wire);
		    onedigitperthread.emplace_back(rawdigits);



                } // <---End loop over gaussians

                fChi2->Fill(chi2PerNDF);

           }//<---End loop over merged candidate hits

       } //<---End looping over ROI's

      hitsthreads[tid]   = onehitperthread;
      wiresthreads[tid]  = onewireperthread;
      digitsthreads[tid] = onedigitperthread;

      //      std::cout << "finished with OMP thread " << tid << std::endl;

    }//<---End looping over all the wires
    } // end of omp parallel region


    // Do not let any rando threads touch anything after this. Lock goes out of scope at end of produce().
    // Entirely placebo. Not obvious this should, or in fact, does do anything.
    std::mutex protect;
    std::lock_guard<std::mutex> lock(protect);

    for(auto const& jj : t_ID ){
      //    std::vector <int> t_IDtmp(0);
      //std::iota (std::begin(t_IDtmp), std::end(t_IDtmp), 0); // fill it with 0,1,2, .... 8255
      //for(auto const& jj : t_IDtmp ){

      auto  hind = hitsthreads.find(jj);
      auto  wind = wiresthreads.find(jj);
      auto  dind = digitsthreads.find(jj);
      if (hind == hitsthreads.end()) break;
      if (wind == wiresthreads.end()) break;
      if (dind == digitsthreads.end()) break;
      for (uint16_t kk=0 ; kk < (hind->second).size(); kk++ ) {
	//	std::cout << "GausHitFinder: emplace_back'ing hit " << kk << " on wire " << jj << std::endl;
	hcol.emplace_back(hind->second.at(kk),wind->second.at(kk),dind->second.at(kk));
	}
    }

    //    std::cout << "\n\t Placing " << hcol.size() << " hits onto the collection. " << std::endl;

    //==================================================================================================
    // End of the event

    // move the hit collection and the associations into the event
    hcol.put_into(evt);


    stp1 = clock();
    //    FILE *fd;
    //    fd = fopen("timeOut.txt","w");
    //    fprintf(fd, "RunTime:%f\n", (double)(stp1-str1)/CLOCKS_PER_SEC);
    std::cout << "This event's run time is " <<  (double)(stp1-str1)/CLOCKS_PER_SEC  << std::endl;

} // End of produce()

// --------------------------------------------------------------------------------------------
// Initial finding of candidate peaks
// --------------------------------------------------------------------------------------------
void hit::GausHitFinder::findCandidatePeaks(std::vector<float>::const_iterator    startItr,
                                            std::vector<float>::const_iterator    stopItr,
                                            std::vector<std::tuple<int,int,int>>& timeValsVec,
                                            float&                                roiThreshold,
                                            int                                   firstTick) const
{
    // Need a minimum number of ticks to do any work here
    if (std::distance(startItr,stopItr) > 4)
    {
        // Find the highest peak in the range given
        auto maxItr = std::max_element(startItr, stopItr);

        float maxValue = *maxItr;
        int   maxTime  = std::distance(startItr,maxItr);

        if (maxValue > roiThreshold)
        {
            // backwards to find first bin for this candidate hit
            auto firstItr = std::distance(startItr,maxItr) > 2 ? maxItr - 1 : startItr;

            while(firstItr != startItr)
            {
                // Check both sides of firstItr and look for min/inflection point
                if (*firstItr < *(firstItr+1) && *firstItr <= *(firstItr-1)) break;

                firstItr--;
            }

            int firstTime = std::distance(startItr,firstItr);

            // Recursive call to find all candidate hits earlier than this peak
            findCandidatePeaks(startItr, firstItr + 1, timeValsVec, roiThreshold, firstTick);

            // forwards to find last bin for this candidate hit
            auto lastItr = std::distance(maxItr,stopItr) > 2 ? maxItr + 1 : stopItr - 1;

            while(lastItr != stopItr - 1)
            {
                // Check both sides of firstItr and look for min/inflection point
                if (*lastItr <= *(lastItr+1) && *lastItr < *(lastItr-1)) break;

                lastItr++;
            }

            int lastTime = std::distance(startItr,lastItr);

            // Now save this candidate's start and max time info
            timeValsVec.push_back(std::make_tuple(firstTick+firstTime,firstTick+maxTime,firstTick+lastTime));

            // Recursive call to find all candidate hits later than this peak
            findCandidatePeaks(lastItr + 1, stopItr, timeValsVec, roiThreshold, firstTick + std::distance(startItr,lastItr + 1));
        }
    }

    return;
}

// --------------------------------------------------------------------------------------------
// Merging of nearby candidate peaks
// --------------------------------------------------------------------------------------------

void hit::GausHitFinder::mergeCandidatePeaks(const std::vector<float>& signalVec, TimeValsVec& timeValsVec, MergedTimeWidVec& mergedVec) const
{
    // ################################################################
    // ### Lets loop over the candidate pulses we found in this ROI ###
    // ################################################################
    auto timeValsVecItr = timeValsVec.begin();

    while(timeValsVecItr != timeValsVec.end())
    {
        PeakTimeWidVec peakVals;

        // Setting the start, peak, and end time of the pulse
        auto& timeVal = *timeValsVecItr++;

        int startT = std::get<0>(timeVal);
        int maxT   = std::get<1>(timeVal);
        int endT   = std::get<2>(timeVal);
        int widT   = std::max(2,(endT - startT) / 6);

        peakVals.emplace_back(maxT,widT);

        // See if we want to merge pulses together
        // First check if we have more than one pulse on the wire
        bool checkNextHit = timeValsVecItr != timeValsVec.end();

        // Loop until no more merged pulses (or candidates in this ROI)
        while(checkNextHit)
        {
            // If the start time of the next pulse is the end time of the current pulse then merge
            // Alternatively, if the start time of the next pulses is one tick away then
            // merge if the intervening signal is above 0.5 * threshold
            int nextStartT = std::get<0>(*timeValsVecItr);

            if ((nextStartT == endT) ||(nextStartT - endT  < 2 && signalVec[endT+1] > threshold/2))
            {
                timeVal = *timeValsVecItr++;
                maxT    = std::get<1>(timeVal);
                endT    = std::get<2>(timeVal);
                widT    = std::max(2,(endT - nextStartT) / 6);

                peakVals.emplace_back(maxT,widT);

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
// Free functions f and df/dx for Gasussians_gsl minimizer.
// --------------------------------------------------------------------------------------------

struct data {
  const size_t npts;
  const int ngauss;
  double * x;
  double * y;
};

int GaussN_f(const gsl_vector * p, void *d, 
        gsl_vector * f)
{

  for (size_t i = 0; i < ((struct data *)d)->npts; i++)
    {
      /* Model Yi = SUM {N exp [(x-mu)^2/2sig^2] }  */

      double Yi(0);
      for (int jj=0; jj<((struct data *)d)->ngauss; jj++)
	{
	  double expg;
	  expg = std::exp( -std::pow( ((struct data *)d)->x[i] - gsl_vector_get(p,1+jj*3) ,2.)/(2*std::pow(gsl_vector_get(p,2+jj*3),2.)) );
	  Yi+=gsl_vector_get(p,0+jj*3)*expg;
	}
      
      gsl_vector_set (f, i, Yi - ((struct data *)d)->y[i]);
      //      std::cout << "GaussN_f: func, data at point "  << ((struct data *)d)->x[i] << " are " << Yi << ", " << ((struct data *)d)->y[i] << std::endl;
    }

  return GSL_SUCCESS;
}

int GaussN_df (const gsl_vector * p, void *d, 
         gsl_matrix * J)
{

  for (size_t i = 0; i < ((struct data *)d)->npts; i++)
    {
      double dfidP0(0);
      double dfidP1(0);
      double dfidP2(0);

      for (int jj=0; jj<((struct data *)d)->ngauss; jj++)
	{
	  double expg;
	  expg = std::exp(  -std::pow( ((struct data *)d)->x[i] - gsl_vector_get(p,1+jj*3) ,2.)/(2*std::pow(gsl_vector_get(p,2+jj*3),2.)) );

	  dfidP0 = expg;
	  dfidP1 = gsl_vector_get(p,0+jj*3)/std::pow( gsl_vector_get(p,2+jj*3), 2.) * (((struct data *)d)->x[i] - gsl_vector_get(p,1+jj*3)) * expg;
	  dfidP2 = gsl_vector_get(p,0+jj*3)/std::pow( gsl_vector_get(p,2+jj*3), 3.) * std::pow(((struct data *)d)->x[i] - gsl_vector_get(p,1+jj*3), 2.) * expg;

	  gsl_matrix_set (J, i, 0+jj*3, dfidP0);
	  gsl_matrix_set (J, i, 1+jj*3, dfidP1);
	  gsl_matrix_set (J, i, 2+jj*3, dfidP2);
	  //	  std::cout << "GaussN_df: Jacobian elements at point "  << ((struct data *)d)->x[i] << " are " << dfidP0 << ", " << dfidP1 << ", " << dfidP2 << std::endl;
	}
    }

  return GSL_SUCCESS;
}


// --------------------------------------------------------------------------------------------
// Fit Gaussians
// --------------------------------------------------------------------------------------------
void hit::GausHitFinder::FitGaussians_gsl(const std::vector<float>& SignalVector,
                                      const PeakTimeWidVec&     PeakVals,
                                      int                       StartTime,
                                      int                       EndTime,
                                      double                    ampScaleFctr,
                                      ParameterVec&             paramVec,
                                      double&                   chi2PerNDF,
                                      int&                      NDF,
				      const int&                t
				      )
{
    int size = EndTime - StartTime;
    // #############################################
    // ### If size < 0 then set the size to zero ###
    // #############################################
    if(EndTime - StartTime < 0){size = 0;}

    // --- TH1D HitSignal ---
    TH1F hitSignal("hitSignal","",std::max(size,1),StartTime,EndTime);
    /////fhitSignal[t] = new TH1F("hitSignal","",std::max(size,1),StartTime,EndTime);
    ///fhitSignal[t].SetNameTitle(("hitSignal"+std::to_string(t)).c_str(),"");
    ///fhitSignal[t].SetBins(std::max(size,1),StartTime,EndTime);
    hitSignal.Sumw2();

    // #############################
    // ### Filling the histogram ###
    // #############################

    gsl_vector * start;
    const int N(PeakVals.size());
    const size_t Ndatapoints(size);
    std::vector<double> wts ; wts.reserve(Ndatapoints);
    std::vector<double> xval; xval.reserve(Ndatapoints);
    std::vector<double> yval; yval.reserve(Ndatapoints);

    start = gsl_vector_alloc (N*3);

    for(int aa = StartTime; aa < EndTime; aa++)
    {
      //      std::cout<<"Time Tick = "<<aa<<" (of " << size << "ticks), ADC = "<<SignalVector[aa]<<std::endl;
      hitSignal.Fill(aa,SignalVector[aa]);
      wts.push_back ((double) std::sqrt(std::abs(SignalVector[aa])));
      xval.push_back((double) aa);
      yval.push_back((double) SignalVector[aa]);
      //      std::cout<<"\t "<<" wts, xval at " << aa << " are "<<wts.at(aa) <<", " << xval.at(aa) <<std::endl;
      
      if(EndTime > 10000){break;} // FIXME why?
    }//<---End aa loop


    int status, info;
    int max_iter = 100;


    // ### Setting the parameters for the Gaussian Fit ###
    int parIdx(0);
    for(auto& peakVal : PeakVals)
    {
        double peakMean   = peakVal.first;
        double peakWidth  = peakVal.second;
        double amplitude  = ampScaleFctr * SignalVector[peakMean];
	///        double meanLowLim = std::max(peakMean - 2.*peakWidth, double(StartTime));
        /// double meanHiLim  = std::min(peakMean + 2.*peakWidth, double(EndTime));


	gsl_vector_set (start, parIdx, amplitude);
	gsl_vector_set (start,1+parIdx, peakMean);
	gsl_vector_set (start,2+parIdx, peakWidth);

	/// Not obvious one can set limits with gsl.

        parIdx += 3;

    }//<---End bb loop

    // ####################################################
    // ### PERFORMING THE TOTAL GAUSSIAN FIT OF THE HIT ###
    // ####################################################


    struct data d = { Ndatapoints, N, (double *) xval.data(), (double *)(yval.data()) };

    gsl_multifit_function_fdf GaussN;

    GaussN.f = &(hit::GaussN_f);
    GaussN.df = &(hit::GaussN_df);
    GaussN.n = Ndatapoints;
    GaussN.p = N*3;
    GaussN.params = &d;

    const gsl_multifit_fdfsolver_type * T 
      = gsl_multifit_fdfsolver_lmder;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;
    
    gsl_multifit_fdfsolver * s = gsl_multifit_fdfsolver_alloc (T, Ndatapoints, 3*N);

    // I give the wts, the fdf, and the initial guess x.
    gsl_vector_view w = gsl_vector_view_array(wts.data(), Ndatapoints);
    status = gsl_multifit_fdfsolver_wset (s, &GaussN, start, &w.vector);



    gsl_matrix * J = gsl_matrix_alloc(Ndatapoints, 3*N);
    gsl_matrix *covar = gsl_matrix_alloc (3*N, 3*N);

    bool ReturnCode = false;

    try
      {
	///hitSignal.Fit(fGaus[t],"QNRWB","", StartTime, EndTime);
	//	std::cout << "calling gsl_multifit_fdsolver ..." << std::endl;
	
	status = gsl_multifit_fdfsolver_driver(s, max_iter, xtol, gtol, ftol, &info);
	if (!status) ReturnCode = true;
	//	std::cout << "fit status is " << status << std::endl;
	
	gsl_multifit_fdfsolver_jac(s, J);
	gsl_multifit_covar (J, 0.0, covar);
	
	// if (!ReturnCode)
	//	  std::cout << "Fitter failed." << std::endl;
      }
    catch(...)
      {mf::LogWarning("GausHitFinder") << "Fitter failed finding a hit";}

    // ##################################################
    // ### Getting the fitted parameters from the fit ###
    // ##################################################

    NDF = Ndatapoints - 3*N;
    chi2PerNDF = std::pow(gsl_blas_dnrm2(gsl_multifit_fdfsolver_residual(s)),2.)/NDF;
    /*
    fprintf(stderr, "summary from method '%s'\n",
	    gsl_multifit_fdfsolver_name(s));
    fprintf(stderr, "number of iterations: %zu\n",
	    gsl_multifit_fdfsolver_niter(s));
    fprintf(stderr, "function evaluations: %zu\n", GaussN.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", GaussN.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
          (info == 1) ? "small step size" : "small gradient");
    */

    //    if (ReturnCode)  std::cout << "GSL fit NDF,chi2PerNDF: " << NDF << ", " << chi2PerNDF  <<  std::endl;

    for(size_t ipar = 0; ipar < (3 * PeakVals.size()); ++ipar)
      {
	if (ReturnCode)
	  {

	    //	    std::cout << "GSL fit ipar param, err: " << ipar << ", " << gsl_vector_get (s->x, 0+ipar) << ", " << std::sqrt(gsl_matrix_get (covar, 0+ipar, 0+ipar))<< std::endl;

	    paramVec.emplace_back(gsl_vector_get (s->x, 0+ipar),std::sqrt(gsl_matrix_get(covar,0+ipar,0+ipar)));
	  }
	else
	  paramVec.emplace_back(0.,0.);
	  //std::cout << PeakVals.size() << " Gaussian-printout: Parameter number/val/error " << ipar << "/" << TheFitter.Result().Parameter(ipar)<< "/" << TheFitter.Result().Error(ipar)  << std::endl;
      }

    gsl_multifit_fdfsolver_free (s);
    gsl_vector_free (start);
    gsl_matrix_free (covar);
    gsl_matrix_free (J);

    ///fGaus[t]->Delete();
    hitSignal.Delete();

}//<----End FitGaussians_gsl


void hit::GausHitFinder::doBinAverage(const std::vector<float>& inputVec,
                                      std::vector<float>&       outputVec,
                                      size_t                    binsToAverage) const
{
    size_t halfBinsToAverage(binsToAverage/2);

    float runningSum(0.);

    for(size_t idx = 0; idx < halfBinsToAverage; idx++) runningSum += inputVec[idx];

    outputVec.resize(inputVec.size());
    std::vector<float>::iterator outputVecItr = outputVec.begin();

    // First pass through to build the erosion vector
    for(std::vector<float>::const_iterator inputItr = inputVec.begin(); inputItr != inputVec.end(); inputItr++)
    {
        size_t startOffset = std::distance(inputVec.begin(),inputItr);
        size_t stopOffset  = std::distance(inputItr,inputVec.end());
        size_t count       = std::min(2 * halfBinsToAverage, std::min(startOffset + halfBinsToAverage + 1, halfBinsToAverage + stopOffset - 1));

        if (startOffset >= halfBinsToAverage) runningSum -= *(inputItr - halfBinsToAverage);
        if (stopOffset  >  halfBinsToAverage) runningSum += *(inputItr + halfBinsToAverage);

        *outputVecItr++ = runningSum / float(count);
    }

    return;
}


void hit::GausHitFinder::reBin(const std::vector<float>& inputVec,
                               std::vector<float>&       outputVec,
                               size_t                    nBinsToCombine) const
{
    size_t nNewBins = inputVec.size() / nBinsToCombine;

    if (inputVec.size() % nBinsToCombine > 0) nNewBins++;

    outputVec.resize(nNewBins, 0.);

    size_t outputBin = 0;

    for(size_t inputIdx = 0; inputIdx < inputVec.size();)
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


  DEFINE_ART_MODULE(GausHitFinder)

} // end of hit namespace
#endif // GAUSHITFINDER_H
