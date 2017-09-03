/**
 *  @file   Cluster3D_module.cc
 *
 *          Class:       Cluster3D
 *          Module Type: Producer
 * 
 *  @brief  Producer module to create 3D clusters from input recob::Hit objects
 *
 *          This producer module will drive the 3D association of recob::Hit objects
 *          to form 3D clusters. This information will be output as:
 *          1) a PFParticle to anchor all the other objects (as associations)
 *          2) three recob::Cluster objects representing 2D hit clusterswith associations 
 *             to the 2D hits comprising each of them
 *          3) One or more recob::PCAxis objects representing the Principal Components
 *             Analysis output of the space points associated to the 3D objects
 *          4) recob::SpacePoints representing the accepted 3D points for each PFParticle
 *          5) recob::Seed objects and associated seed hits representing candidate straight
 *             line segments in the space point collection. 
 *
 *          The module has two main sections
 *          1) Find the 3D clusters of 3D hits
 *          2) Get the output objects for each of the 3D clusters
 *          See the code below for more detail on these steps
 *
 *          Note that the general 3D cluster suite of algorithms make extensive use of a set of data objects 
 *          which contain volatile data members. At the end of the routine these are used to make the output
 *          LArSoft data products described above. See LarData/RecoObjects/Cluster3D.h
 *   
 *          Configuration parameters:
 *          HitFinderModuleLabel:         the producer module responsible for making the recob:Hits to use
 *          EnableMonitoring:             if true then basic monitoring of the module performed
 *          DBScanAlg:                    Parameter block required by the 3D clustering algorithm
 *          PrincipalComponentsAlg:       Parameter block required by the Principal Components Analysis Algorithm
 *          SkeletonAlg:                  Parameter block required by the 3D skeletonization algorithm
 *          SeedFinderAlg:                Parameter block required by the Hough Seed Finder algorithm
 *          PCASeedFinderAlg:             Parameter block required by the PCA Seed Finder algorithm
 *          ParrallelHitsAlg:             Parameter block required by the parallel hits algorithm
 *
 *          The current producer module does not try to analyze or break apart PFParticles
 *          so, for example, all tracks emanating from a common vertex will be associated
 *          to a single PFParticle
 *
 *  @author usher@slac.stanford.edu 
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardata/RecoObjects/Cluster3D.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

#include "larreco/RecoAlg/Cluster3DAlgs/HoughSeedFinderAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PCASeedFinderAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/ParallelHitsSeedFinderAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/SkeletonAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/DBScanAlg.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterRecoUtil/OverriddenClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "larreco/ClusterFinder/ClusterCreator.h"

// ROOT includes
#include "TTree.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d
{
//------------------------------------------------------------------------------------------------------------------------------------------
// Start with typedefs and definitions of some utility classes

typedef std::map<const recob::Hit*, art::Ptr<recob::Hit> > RecobHitToPtrMap;
typedef std::vector<art::Ptr<recob::Hit> >                 RecobHitVector;
    
typedef std::pair<reco::PrincipalComponents, reco::HitPairClusterMap::iterator> PCAHitPairClusterMapPair;
    
/**
 *  @brief A utility class used in construction of 3D clusters
 */
class RecobClusterParameters
{
public:
    RecobClusterParameters() : m_startTime(999999.),
        m_sigmaStartTime(1.),
        m_endTime(0.),
        m_sigmaEndTime(1.),
        m_totalCharge(0.),
        m_startWire(9999999),
        m_endWire(0),
        m_view(geo::kUnknown)
    {
        m_hitVector.clear();
    }
        
    void UpdateParameters(const reco::ClusterHit2D* hit);
        
    double         m_startTime;
    double         m_sigmaStartTime;
    double         m_endTime;
    double         m_sigmaEndTime;
    double         m_totalCharge;
    unsigned int   m_startWire;
    unsigned int   m_endWire;
    geo::View_t    m_view;
    HitVectorConst m_hitVector;
};
    
typedef std::map<geo::View_t, RecobClusterParameters> ViewToClusterParamsMap;

/**
 *  @brief Class wrapping the above and containing volatile information to characterize the cluster
 */
class ClusterParameters
{
public:
    ClusterParameters(reco::HitPairClusterMap::iterator& mapItr) : m_hitPairListPtr(mapItr->second)
    {
        m_clusterParams.clear();
    }
    
    ClusterParameters(reco::HitPairListPtr& hitList) : m_hitPairListPtr(hitList)
    {
        m_clusterParams.clear();
    }
    
    void UpdateParameters(const reco::ClusterHit2D* hit)
    {
        m_clusterParams[hit->getHit().View()].UpdateParameters(hit);
    }
    
    ViewToClusterParamsMap    m_clusterParams;
    reco::HitPairListPtr&     m_hitPairListPtr;
    reco::PrincipalComponents m_fullPCA;
    reco::PrincipalComponents m_skeletonPCA;
};
    
typedef std::list<ClusterParameters> ClusterParametersList;

//------------------------------------------------------------------------------------------------------------------------------------------
// Definition of the producer module here
    
/**
 *  @brief  Definition of the Cluster3D class
 */
class Cluster3D : public art::EDProducer
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset - reference to the parameters used by this module and its algorithms
     */
    Cluster3D(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~Cluster3D();

    /**
     *  @brief declare the standard art functions that we'll implement in this producer module
     */
    void beginJob();
    void endJob();
    void produce(art::Event &evt);
    void reconfigure(fhicl::ParameterSet const &pset);

private:

    /**
     *  @brief  Event Preparation 
     * 
     *  @param  evt  the ART event 
     */
    void PrepareEvent(const art::Event &evt);  

    /**
     *  @brief  Extract the ART hits and the ART hit-particle relationships
     * 
     *  @param  evt                   the ART event
     *  @param  hit2DVector           A container for the internal Cluster3D 2D hit objects
     *  @param  viewToHitVectorMap    A map between view and the internal Cluster3D 2D hit objects
     *  @param  viewToWireToHitSetMap This maps 2D hits to wires and stores by view
     *  @param  hitToPtrMap           This maps our Cluster2D hits back to art Ptr's to reco Hits
     */
    void CollectArtHits(art::Event&            evt,
                        Hit2DVector&           hit2DVector,
                        ViewToHitVectorMap&    viewToHitVector,
                        ViewToWireToHitSetMap& viewToWireToHitSetMap,
                        RecobHitToPtrMap&      hitToPtrMap) const;

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();
    
    /**
     *  @brief Given the results of running DBScan, format the clusters so that they can be 
     *         easily transferred back to the larsoft world
     *
     *  @param hitPairClusterMap      map between view and a list of 3D hits
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     *
     *                                The last two parameters are passed through to the FillClusterParams method
     */
    void BuildClusterInfo(reco::HitPairClusterMap& hitPairClusterMap, ClusterParametersList& clusterParametersList, double rejectFraction = 0.5) const;
    
    /**
     *  @brief A generic routine to actually fill the clusterParams
     *
     *  @param clusterParametersList  a container for our candidate 3D clusters
     *  @param rejectionFraction      Used for determine "hit purity" when rejecting clusters
     */
    void FillClusterParams(ClusterParameters& clusterParams, double rejectFraction = 0.5, double maxLostRatio = 0.75) const;
    
    /**
     *  @brief An interface to the seed finding algorithm
     *
     *  @param evt          the ART event
     *  @param cluster      structure of information representing a single cluster
     *  @param hitToPtrMap  This maps our Cluster2D hits back to art Ptr's to reco Hits
     *  @param seedVec      the output vector of candidate seeds
     *  @param seedHitAssns the associations between the seeds and the 2D hits making them
     */
    void findTrackSeeds(art::Event&                         evt,
                        ClusterParameters&                  cluster,
                        RecobHitToPtrMap&                   hitToPtrMap,
                        std::vector<recob::Seed>&           seedVec,
                        art::Assns<recob::Seed,recob::Hit>& seedHitAssns) const;
    
    /**
     *  @brief Attempt to split clusters by using a minimum spanning tree
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    void splitClustersWithMST(ClusterParameters& clusterParameters, ClusterParametersList& clusterParametersList) const;
    
    /**
     *  @brief Attempt to split clusters using the output of the Hough Filter
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    void splitClustersWithHough(ClusterParameters&       clusterParameters,
                                reco::HitPairClusterMap& hitPairClusterMap,
                                ClusterParametersList&   clusterParametersList) const;
    
    /**
     *  @brief Produces the art output from all the work done in this producer module
     *
     *  @param evt                   the ART event
     *  @param hitPairList           List of all 3D Hits in internal Cluster3D format
     *  @param clusterParametersList Data structure containing the cluster information to output
     *  @param  hitToPtrMap           This maps our Cluster2D hits back to art Ptr's to reco Hits
     */
    void ProduceArtClusters(art::Event&              evt,
                            HitPairList&             hitPairList,
                            reco::HitPairClusterMap& hitPairClusterMap,
                            ClusterParametersList&   clusterParametersList,
                            RecobHitToPtrMap&        hitToPtrMap) const;
    
    /**
     *  @brief There are several places we will want to know if a candidate cluster is a
     *         "parallel hits" type of cluster so encapsulate that here.
     *
     *  @param pca  The Principal Components Analysis parameters for the cluster
     */
    bool aParallelHitsCluster(const reco::PrincipalComponents& pca) const
    {
        return fabs(pca.getEigenVectors()[2][0]) > m_parallelHitsCosAng && 3. * sqrt(pca.getEigenValues()[1]) > m_parallelHitsTransWid;
    }

    /**
     *   Algorithm parameters
     */
    bool                      m_enableMonitoring;      ///< Turn on monitoring of this algorithm
    std::string               m_hitfinderModuleLabel;  ///< Producer of the reco hits
    double                    m_clusHitRejectionFrac;  ///< Cluster hit purity must exceed this to be kept
    double                    m_parallelHitsCosAng;    ///< Cut for PCA 3rd axis angle to X axis
    double                    m_parallelHitsTransWid;  ///< Cut on transverse width of cluster (PCA 2nd eigenvalue)

    /**
     *   Tree variables for output
     */
    TTree*                    m_pRecoTree;             ///<
    int                       m_run;                   ///<
    int                       m_event;                 ///<
    int                       m_hits;                  ///< Keeps track of the number of hits seen
    float                     m_totalTime;             ///< Keeps track of total execution time
    float                     m_artHitsTime;           ///< Keeps track of time to recover hits
    float                     m_makeHitsTime;          ///< Keeps track of time to build 3D hits
    float                     m_buildNeighborhoodTime; ///< Keeps track of time to build epsilon neighborhood
    float                     m_dbscanTime;            ///< Keeps track of time to run DBScan
    float                     m_finishTime;            ///< Keeps track of time to run output module
    
    /** 
     *   Other useful variables
     */
    geo::Geometry*            m_geometry;              ///<  pointer to the Geometry service
  const detinfo::DetectorProperties* m_detector;              ///<  Pointer to the detector properties
    
    DBScanAlg                 m_dbScanAlg;             ///<  Algorithm to cluster hits
    PrincipalComponentsAlg    m_pcaAlg;                ///<  Principal Components algorithm
    SkeletonAlg               m_skeletonAlg;           ///<  Skeleton point finder
    HoughSeedFinderAlg        m_seedFinderAlg;         ///<  Seed finder
    PCASeedFinderAlg          m_pcaSeedFinderAlg;      ///<  Use PCA axis to find seeds
    ParallelHitsSeedFinderAlg m_parallelHitsAlg;       ///<  Deal with parallel hits clusters
};

DEFINE_ART_MODULE(Cluster3D)

} // namespace lar_cluster3d

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

Cluster3D::Cluster3D(fhicl::ParameterSet const &pset) :
    m_dbScanAlg(pset.get<fhicl::ParameterSet>("DBScanAlg")),
    m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg")),
    m_skeletonAlg(pset.get<fhicl::ParameterSet>("SkeletonAlg")),
    m_seedFinderAlg(pset.get<fhicl::ParameterSet>("SeedFinderAlg")),
    m_pcaSeedFinderAlg(pset.get<fhicl::ParameterSet>("PCASeedFinderAlg")),
    m_parallelHitsAlg(pset.get<fhicl::ParameterSet>("ParallelHitsAlg"))
{
    this->reconfigure(pset);

    produces< std::vector<recob::PCAxis> >();
    produces< std::vector<recob::PFParticle> >();
    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::SpacePoint> >();
    produces< std::vector<recob::Seed> >();
    produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
    produces< art::Assns<recob::PFParticle, recob::Seed> >();
    produces< art::Assns<recob::Seed,       recob::Hit> >();
    produces< art::Assns<recob::Cluster,    recob::Hit> >();
    produces< art::Assns<recob::SpacePoint, recob::Hit> >();
}

//------------------------------------------------------------------------------------------------------------------------------------------

Cluster3D::~Cluster3D()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::reconfigure(fhicl::ParameterSet const &pset)
{
    m_hitfinderModuleLabel = pset.get<std::string>("HitFinderModuleLabel", "gaushit");
    m_enableMonitoring     = pset.get<bool>       ("EnableMonitoring",         false);
    m_clusHitRejectionFrac = pset.get<double>     ("ClusterHitRejectionFrac",    0.5);
    m_parallelHitsCosAng   = pset.get<double>     ("ParallelHitsCosAng",       0.999);
    m_parallelHitsTransWid = pset.get<double>     ("ParallelHitsTransWid",      25.0);
    
    m_dbScanAlg.reconfigure(pset.get<fhicl::ParameterSet>("DBScanAlg"));
    m_pcaAlg.reconfigure(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"));
    m_skeletonAlg.reconfigure(pset.get<fhicl::ParameterSet>("SkeletonAlg"));
    m_seedFinderAlg.reconfigure(pset.get<fhicl::ParameterSet>("SeedFinderAlg"));
    m_pcaSeedFinderAlg.reconfigure(pset.get<fhicl::ParameterSet>("PCASeedFinderAlg"));
    m_parallelHitsAlg.reconfigure(pset.get<fhicl::ParameterSet>("ParallelHitsAlg"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::beginJob()
{
    /**
     *  @brief beginJob will be tasked with initializing monitoring, in necessary, but also to init the 
     *         geometry and detector services (and this probably needs to go in a "beginEvent" method?)
     */
    if (m_enableMonitoring)
        this->InitializeMonitoring();
    
    art::ServiceHandle<geo::Geometry>            geometry;
    
    m_geometry = &*geometry;
    m_detector = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::produce(art::Event &evt)
{
    /**
     *  @brief Producer method for reovering the 2D hits and driving the 3D reconstruciton
     */
    mf::LogInfo("Cluster3D") << " *** Cluster3D::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] Starting Now! *** " << std::endl;

    // Set up for monitoring the timing... at some point this should be removed in favor of
    // external profilers
    cet::cpu_timer theClockTotal;
    cet::cpu_timer theClockArtHits;
    cet::cpu_timer theClockFinish;
    
    if (m_enableMonitoring)
    {
        theClockTotal.start();
        theClockArtHits.start();
    }
    
    // This really only does anything if we are monitoring since it clears our tree variables
    this->PrepareEvent(evt);

    // Get instances of the primary data structures needed
    Hit2DVector                  clusterHit2DMasterVec;
    ViewToHitVectorMap           viewToHitVectorMap;
    ViewToWireToHitSetMap        viewToWireToHitSetMap;
    reco::HitPairClusterMap      hitPairClusterMap;
    ClusterParametersList        clusterParametersList;
    RecobHitToPtrMap             clusterHitToArtPtrMap;
    std::unique_ptr< HitPairList > hitPairList(new HitPairList); // Potentially lots of hits, use heap instead of stack
    
    // Recover the 2D hits and then organize them into data structures which will be used in the
    // DBscan algorithm for building the 3D clusters
    this->CollectArtHits(evt, clusterHit2DMasterVec, viewToHitVectorMap, viewToWireToHitSetMap, clusterHitToArtPtrMap);
    
    if (m_enableMonitoring) theClockArtHits.stop();
    
    // If there are no hits in our view/wire data structure then do not proceed with the full analysis
    if (!viewToWireToHitSetMap.empty())
    {
        // Call the main workhorse algorithm for building the local version of candidate 3D clusters
        m_dbScanAlg.ClusterHitsDBScan(viewToHitVectorMap, viewToWireToHitSetMap, *hitPairList, hitPairClusterMap);
        
        // Given the work above, process and build the list of 3D clusters to output
        BuildClusterInfo(hitPairClusterMap, clusterParametersList, m_clusHitRejectionFrac);
    }
    
    if(m_enableMonitoring) theClockFinish.start();

    // Call the module that does the end processing (of which there is quite a bit of work!)
    // This goes here to insure that something is always written to the data store
    ProduceArtClusters(evt, *hitPairList, hitPairClusterMap, clusterParametersList, clusterHitToArtPtrMap);
    
    if (m_enableMonitoring) theClockFinish.stop();
    
    // If monitoring then deal with the fallout
    if (m_enableMonitoring)
    {
        theClockTotal.stop();

        m_run                   = evt.run();
        m_event                 = evt.id().event();
        m_totalTime             = theClockTotal.accumulated_real_time();
        m_artHitsTime           = theClockArtHits.accumulated_real_time();
        m_makeHitsTime          = m_dbScanAlg.getTimeToExecute(DBScanAlg::BUILDTHREEDHITS);
        m_buildNeighborhoodTime = m_dbScanAlg.getTimeToExecute(DBScanAlg::BUILDHITTOHITMAP);
        m_dbscanTime            = m_dbScanAlg.getTimeToExecute(DBScanAlg::RUNDBSCAN);
        m_finishTime            = theClockFinish.accumulated_real_time();
        m_hits                  = static_cast<int>(clusterHit2DMasterVec.size());
        m_pRecoTree->Fill();
        
        mf::LogDebug("Cluster3D") << "*** Cluster3D total time: " << m_totalTime << ", art: " << m_artHitsTime << ", make: " << m_makeHitsTime
        << ", build: " << m_buildNeighborhoodTime << ", dbscan: " << m_dbscanTime << ", finish: " << m_finishTime << std::endl;
    }
    
    // Will we ever get here? ;-)
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run",                  &m_run,                   "run/I");
    m_pRecoTree->Branch("event",                &m_event,                 "event/I");
    m_pRecoTree->Branch("hits",                 &m_hits,                  "hits/I");
    m_pRecoTree->Branch("totalTime",            &m_totalTime,             "time/F");
    m_pRecoTree->Branch("artHitsTime",          &m_artHitsTime,           "time/F");
    m_pRecoTree->Branch("makeHitsTime",         &m_makeHitsTime,          "time/F");
    m_pRecoTree->Branch("buildneigborhoodTime", &m_buildNeighborhoodTime, "time/F");
    m_pRecoTree->Branch("dbscanTime",           &m_dbscanTime,            "time/F");
    m_pRecoTree->Branch("finishTime",           &m_finishTime,            "time/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Cluster3D::PrepareEvent(const art::Event &evt)
{
    m_run                   = evt.run();
    m_event                 = evt.id().event();
    m_hits                  = 0;
    m_totalTime             = 0.f;
    m_artHitsTime           = 0.f;
    m_makeHitsTime          = 0.f;
    m_buildNeighborhoodTime = 0.f;
    m_dbscanTime            = 0.f;
    m_finishTime            = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool SetHitTimeOrder(const reco::ClusterHit2D* left, const reco::ClusterHit2D* right)
{
    // Sort by "modified start time" of pulse
    return left->getHit().PeakTime() < right->getHit().PeakTime();
}
    
bool Hit2DSetCompare::operator() (const reco::ClusterHit2D* left, const reco::ClusterHit2D* right) const
{
    return left->getHit().PeakTime() < right->getHit().PeakTime();
}

//------------------------------------------------------------------------------------------------------------------------------------------
void Cluster3D::CollectArtHits(art::Event&            evt,
                               Hit2DVector&           hitVector,
                               ViewToHitVectorMap&    viewToHitVectorMap,
                               ViewToWireToHitSetMap& viewToWireToHitSetMap,
                               RecobHitToPtrMap&      hitToPtrMap) const
{
    /**
     *  @brief Recover the 2D hits from art and fill out the local data structures for the 3D clustering
     */
    art::Handle< std::vector<recob::Hit> > recobHitHandle;
    evt.getByLabel(m_hitfinderModuleLabel, recobHitHandle);
    
    if (!recobHitHandle.isValid()) return;
    
    // We'll need the offsets for each plane
    std::map<geo::View_t, double> viewOffsetMap;
    
    viewOffsetMap[geo::kU] = m_detector->GetXTicksOffset(geo::kU, 0, 0)-m_detector->TriggerOffset();
    viewOffsetMap[geo::kV] = m_detector->GetXTicksOffset(geo::kV, 0, 0)-m_detector->TriggerOffset();
    viewOffsetMap[geo::kW] = m_detector->GetXTicksOffset(geo::kW, 0, 0)-m_detector->TriggerOffset();
    
    // Reserve memory for the hit vector
    hitVector.reserve(recobHitHandle->size());
    
    // Cycle through the recob hits to build ClusterHit2D objects and insert
    // them into the map
    for (size_t cIdx = 0; cIdx < recobHitHandle->size(); cIdx++)
    {
        art::Ptr<recob::Hit> recobHit(recobHitHandle, cIdx);
        
        const geo::WireID& hitWireID(recobHit->WireID());
        
        double hitPeakTime(recobHit->PeakTime() - viewOffsetMap[recobHit->View()]);
        double xPosition(m_detector->ConvertTicksToX(recobHit->PeakTime(), hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat));
        
        hitVector.emplace_back(reco::ClusterHit2D(0, 0., 0., xPosition, hitPeakTime, *recobHit));

        viewToHitVectorMap[recobHit->View()].push_back(&hitVector.back());
        viewToWireToHitSetMap[recobHit->View()][recobHit->WireID().Wire].insert(&hitVector.back());
        
        const recob::Hit* recobHitPtr = recobHit.get();
        hitToPtrMap[recobHitPtr]      = recobHit;
    }
    
    // Make a loop through to sort the recover hits in time order
    std::sort(viewToHitVectorMap[geo::kU].begin(), viewToHitVectorMap[geo::kU].end(), SetHitTimeOrder);
    std::sort(viewToHitVectorMap[geo::kV].begin(), viewToHitVectorMap[geo::kV].end(), SetHitTimeOrder);
    std::sort(viewToHitVectorMap[geo::kW].begin(), viewToHitVectorMap[geo::kW].end(), SetHitTimeOrder);

    mf::LogDebug("Cluster3D") << ">>>>> Number of ART hits: " << recobHitHandle->size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void RecobClusterParameters::UpdateParameters(const reco::ClusterHit2D* clusterHit)
{
    /**
     *  @brief a utility routine for building 3D clusters to keep basic info up to date
     *         (a candidate for a better way to do this)
     */
    const recob::Hit& hit = clusterHit->getHit();
    
    // Need to keep track of stuff so we can form cluster
    if (hit.WireID().Wire < m_startWire)
    {
        m_startWire      = hit.WireID().Wire;
        m_startTime      = hit.PeakTimeMinusRMS();
        m_sigmaStartTime = hit.SigmaPeakTime();
    }
    
    if (hit.WireID().Wire > m_endWire)
    {
        m_endWire      = hit.WireID().Wire;
        m_endTime      = hit.PeakTimePlusRMS();
        m_sigmaEndTime = hit.SigmaPeakTime();
    }
    
    m_totalCharge += hit.Integral();
    m_view         = hit.View();
    
    m_hitVector.push_back(clusterHit);
    
    return;
}
    
struct HitPairClusterOrder
{
    bool operator()(const reco::HitPairClusterMap::iterator& left, const reco::HitPairClusterMap::iterator& right)
    {
        // Watch out for the case where two clusters can have the same number of hits!
        if (left->second.size() == right->second.size())
            return left->first < right->first;
        
        return left->second.size() > right->second.size();
    }
};
    
void Cluster3D::BuildClusterInfo(reco::HitPairClusterMap& hitPairClusterMap, ClusterParametersList& clusterParametersList, double rejectFraction) const
{
    /** 
     *  @brief Given a list of a list of candidate cluster hits, build these out into the intermediate 
     *         3D cluster objects to pass to the final stage
     *
     *         Note that this routine will also reject unworthy clusters, in particular those that share too
     *         many hits with other clusters. The criteria is that a larger cluster (more hits) will be superior
     *         to a smaller one, if the smaller one shares too many hits with the larger it is zapped. 
     *         *** THIS IS AN AREA FOR CONTINUED STUDY ***
     */
    
    // This is a remote possibility but why not check?
    if (!hitPairClusterMap.empty())
    {
        size_t minHitsPerCluster(3);
        
        // We want to order our clusters on by largest (most number hits) to smallest. So, we'll loop through the clusters,
        // weeding out the unwanted ones and keep track of things in a set of "good" clusters which we'll order
        // by cluster size.
        std::set<reco::HitPairClusterMap::iterator, HitPairClusterOrder> hitPairClusterSet;
        
        // Loop over the "Clusters" in our map where this loop serves a double purpose
        // In the first we are weeding out clusters which fall below what we think is the minimum number of hits
        // More importantly, we are transferring the cluster ownership to a set so it can order the clusters by
        // number of associated hits
        for (reco::HitPairClusterMap::iterator mapItr = hitPairClusterMap.begin(); mapItr != hitPairClusterMap.end(); mapItr++)
        {
            // Weed out the little people
            if (mapItr->second.size() < minHitsPerCluster) continue;
            
            // Add to the set
            hitPairClusterSet.insert(mapItr);
        }
        
        // What remains is an order set of clusters, largest first
        // Now go through and obtain cluster parameters
        for(std::set<reco::HitPairClusterMap::iterator, HitPairClusterOrder>::iterator setItr = hitPairClusterSet.begin(); setItr != hitPairClusterSet.end(); setItr++)
        {
            // Recover original map iterator
            reco::HitPairClusterMap::iterator hitPairClusterMapItr = *setItr;
            
            // Create a new cluster params object in the vector
            clusterParametersList.push_back(ClusterParameters(hitPairClusterMapItr));
            
            // Can we get a reference to what we just created?
            ClusterParameters& clusterParams = clusterParametersList.back();
            
            // Do the actual work of filling the parameters
            FillClusterParams(clusterParams, rejectFraction);
            
            // If this cluster is rejected then the parameters will be empty
            if (clusterParams.m_clusterParams.empty() || !clusterParams.m_fullPCA.getSvdOK())
            {
                clusterParametersList.pop_back();
            }
        }
    }
    
    return;
}
    
void Cluster3D::FillClusterParams(ClusterParameters& clusterParams, double rejectFraction, double maxLostRatio) const
{
    /**
     *  @brief Given a list of hits fill out the remaining parameters for this cluster and evaluate the
     *         candidate's worthiness to achieve stardom in the event display
     */
    
    // Recover the HitPairListPtr from the input clusterParams (which will be the
    // only thing that has been provided)
    reco::HitPairListPtr& hitPairVector = clusterParams.m_hitPairListPtr;
    
    // To be sure, we should clear the other data members
    clusterParams.m_clusterParams.clear();
    clusterParams.m_fullPCA = reco::PrincipalComponents();
    
    // See if we can avoid duplicates by temporarily transferring to a set
    std::set<const reco::ClusterHit2D*> hitSet;
    
    int numTotal(0);
    int numUniqueHits(0);
    int numSharedHits(0);
    int numLostHits(0);
    
    // Create a list to hold 3D hits which are already in use (criteria below)
    reco::HitPairListPtr usedHitPairList;
    
    // First loop through the 3D hits
    // The goal of this loop is to build a set of unique hits from the hit pairs (which may contain many
    // ambiguous duplicate combinations).
    // The secondary goal is to remove 3D hits marked by hit arbitration to be tossed
    for(const auto& hit3D : hitPairVector)
    {
        int nHitsAlreadyUsed(0);
        
        numTotal += hit3D->getHits().size();
        
        // loop over the hits in this 3D Cluster hit
        for(const auto& hit2D : hit3D->getHits())
        {
            if (hit2D->getStatusBits() & 0x1) nHitsAlreadyUsed++;
        }
        
        numSharedHits += nHitsAlreadyUsed;
        
        if (nHitsAlreadyUsed < int(hit3D->getHits().size())-1)   // was 2
        {
            numUniqueHits += int(hit3D->getHits().size()) - nHitsAlreadyUsed;
            for(const auto& hit2D : hit3D->getHits()) hitSet.insert(hit2D);
        }
        else
        {
            numLostHits += hit3D->getHits().size();
//            hit3D->setStatusBit(reco::ClusterHit3D::REJECTED);
            usedHitPairList.emplace_back(hit3D);
        }
    }
    
    // If we have something left then at this point we make one more check
    // This check is intended to weed out clusters made from isolated groups of ambiguous hits which
    // really belong to a larger cluster
    if (numUniqueHits > 3)
    {
        // Look at reject to accept ratio
        //double rejectToAccept = double(numRejected) / double(numAccepted);
        double acceptRatio = double(numUniqueHits) / double(numTotal);
        double lostRatio   = double(numLostHits)   / double(numTotal);
        
        // Arbitrary rejection criteria... need to understand
        // Anyway, if we get past this we're making a cluster
        //if (rejectToAccept < rejectFraction)
        if (acceptRatio > rejectFraction && lostRatio < maxLostRatio)  // lostRatio cut was 1. - off
        {
            // Why do we need to explicitly define this?
            unsigned int usedBit(0x1);
            
            // Add the "good" hits to our cluster parameters
            for(const auto& hit2D : hitSet)
            {
                hit2D->setStatusBit(usedBit);
                clusterParams.UpdateParameters(hit2D);
            }
            
            // Final selection cut, need at least 3 hits each view
            if (!(clusterParams.m_clusterParams[geo::kU].m_hitVector.size() < 2 ||
                  clusterParams.m_clusterParams[geo::kV].m_hitVector.size() < 2 ||
                  clusterParams.m_clusterParams[geo::kW].m_hitVector.size() < 2)  )
            {
                // First task is to remove the hits already in use
                if (!usedHitPairList.empty())
                {
                    hitPairVector.sort();
                    usedHitPairList.sort();
                    
                    reco::HitPairListPtr::iterator newListEnd =
                    std::set_difference(hitPairVector.begin(),   hitPairVector.end(),
                                        usedHitPairList.begin(), usedHitPairList.end(),
                                        hitPairVector.begin() );
                    
                    hitPairVector.erase(newListEnd, hitPairVector.end());
                }
                
                // First stage of feature extraction runs here
                m_pcaAlg.PCAAnalysis_3D(clusterParams.m_hitPairListPtr, clusterParams.m_fullPCA);
                
                // Must have a valid pca
                if (clusterParams.m_fullPCA.getSvdOK())
                {
                    // If any hits were thrown away, see if we can rescue them
                    if (!usedHitPairList.empty())
                    {
                        double maxDoca = 2. * sqrt(clusterParams.m_fullPCA.getEigenValues()[1]);
                    
                        if (maxDoca < 5.)
                        {
                            size_t curHitVectorSize = hitPairVector.size();
                        
                            m_pcaAlg.PCAAnalysis_calc3DDocas(usedHitPairList, clusterParams.m_fullPCA);
                            
                            for(const auto& hit3D : usedHitPairList)
                                if (hit3D->getDocaToAxis() < maxDoca) hitPairVector.push_back(hit3D);
                        
                            if (hitPairVector.size() > curHitVectorSize)
                                m_pcaAlg.PCAAnalysis_3D(clusterParams.m_hitPairListPtr, clusterParams.m_fullPCA);
                        }
                    }
                
                    // Set the skeleton PCA to make sure it has some value
                    clusterParams.m_skeletonPCA = clusterParams.m_fullPCA;
                }
            }
        }
    }

    return;
}
    
void Cluster3D::findTrackSeeds(art::Event&                         evt,
                               ClusterParameters&                  cluster,
                               RecobHitToPtrMap&                   hitToPtrMap,
                               std::vector<recob::Seed>&           seedVec,
                               art::Assns<recob::Seed,recob::Hit>& seedHitAssns) const
{
    /**
     *  @brief This method provides an interface to various algorithms for finding candiate 
     *         recob::Seed objects and, as well, their candidate related seed hits
     */
    
    // Make sure we are using the right pca
    reco::PrincipalComponents& fullPCA        = cluster.m_fullPCA;
    reco::PrincipalComponents& skeletonPCA    = cluster.m_skeletonPCA;
    reco::HitPairListPtr&      hitPairListPtr = cluster.m_hitPairListPtr;
    reco::HitPairListPtr       skeletonListPtr;

    // We want to work with the "skeleton" hits so first step is to call the algorithm to
    // recover only these hits from the entire input collection
    m_skeletonAlg.GetSkeletonHits(hitPairListPtr, skeletonListPtr);
    
    // Skeleton hits are nice but we can do better if we then make a pass through to "average"
    // the skeleton hits position in the Y-Z plane
    m_skeletonAlg.AverageSkeletonPositions(skeletonListPtr);
    
    SeedHitPairListPairVec seedHitPairVec;
    
    // Some combination of the elements below will be used to determine which seed finding algorithm
    // to pursue below
    double eigenVal0 = 3. * sqrt(skeletonPCA.getEigenValues()[0]);
    double eigenVal1 = 3. * sqrt(skeletonPCA.getEigenValues()[1]);
    double eigenVal2 = 3. * sqrt(skeletonPCA.getEigenValues()[2]);
    double transRMS  = sqrt(std::pow(eigenVal1,2) + std::pow(eigenVal2,2));
    
    bool   foundGoodSeed(false);

    // Choose a method for finding the seeds based on the PCA that was run...
    // Currently we have an ad hoc if-else block which I hope will be improved soon!
    if (aParallelHitsCluster(fullPCA))
    {
        // In this case we have a track moving relatively parallel to the wire plane with lots of
        // ambiguous 3D hits. Your best bet here is to use the "parallel hits" algorithm to get the
        // best axis and seeds
        // This algorithm does not fail (foundGoodSeed will always return true)
        foundGoodSeed = m_parallelHitsAlg.findTrackSeeds(hitPairListPtr, skeletonPCA, seedHitPairVec);
    }
    else if (eigenVal0 > 40. && transRMS < 5.)
    {
        // If the input cluster is relatively "straight" then chances are it is a single straight track,
        // probably a CR muon, and we can simply use the PCA to determine the seed
        // This algorithm will check both "ends" of the input hits and if the angles become inconsistent
        // then it will "fail"
        foundGoodSeed = m_pcaSeedFinderAlg.findTrackSeeds(skeletonListPtr, skeletonPCA, seedHitPairVec);
    }
    
    // In the event the above two methods failed then we hit it with the real seed finder
    if (!foundGoodSeed)
    {
        // If here then we have a complicated 3D cluster and we'll use the hough transform algorithm to
        // return a list of candidate seeds and seed hits
        m_seedFinderAlg.findTrackSeeds(skeletonListPtr, skeletonPCA, seedHitPairVec);
    }

    // Go through the returned lists and build out the art friendly seeds and hits
    for(const auto& seedHitPair : seedHitPairVec)
    {
        seedVec.push_back(seedHitPair.first);
        
        // We use a set here because our 3D hits can share 2D hits
        // The set will make sure we get unique combinations of 2D hits
        std::set<art::Ptr<recob::Hit> > seedHitSet;
        
        for(const auto& hit3D : seedHitPair.second)
        {
            for(const auto& hit2D : hit3D->getHits())
            {
                const recob::Hit* recobHit = &hit2D->getHit();
                
                seedHitSet.insert(hitToPtrMap[recobHit]);
            }
        }
        
        RecobHitVector seedHitVec;
        
        for(const auto& hit2D : seedHitSet) seedHitVec.push_back(hit2D);
        
        util::CreateAssn(*this, evt, seedVec, seedHitVec, seedHitAssns);
    }
    
    return;
}
    
struct Hit3DDistanceOrder
{
    bool operator()(const std::pair<double, const reco::ClusterHit3D*>& left, const std::pair<double, const reco::ClusterHit3D*>& right)
    {
        return left.first < right.first;
    }
};
    
void Cluster3D::splitClustersWithMST(ClusterParameters& clusterParameters, ClusterParametersList& clusterParametersList) const
{
    // This is being left in place for future development. Essentially, it was an attempt to implement
    // a Minimum Spanning Tree as a way to split a particular cluster topology, one where two straight
    // tracks cross closely enought to appear as one cluster. As of Feb 2, 2015 I think the idea is still
    // worth merit so am leaving this module in place for now.
    //
    // If this routine is called then we believe we have a cluster which needs splitting.
    // The way we will do this is to use a Minimum Spanning Tree algorithm to associate all
    // hits together by their distance apart. In theory, we should be able to split the cluster
    // by finding the largest distance and splitting at that point.
    //
    // Typedef some data structures that we will use.
    // Start with the adjacency map
    typedef std::pair<double, const reco::ClusterHit3D*>                DistanceHit3DPair;
    typedef std::list<DistanceHit3DPair >                               DistanceHit3DPairList;
    typedef std::map<const reco::ClusterHit3D*, DistanceHit3DPairList > Hit3DToDistanceMap;
    
    // Now typedef the lists we'll keep
    typedef std::list<const reco::ClusterHit3D*>                        Hit3DList;
    typedef std::pair<Hit3DList::iterator, Hit3DList::iterator>         Hit3DEdgePair;
    typedef std::pair<double, Hit3DEdgePair >                           DistanceEdgePair;
    typedef std::list<DistanceEdgePair >                                DistanceEdgePairList;
    
    struct DistanceEdgePairOrder
    {
        bool operator()(const DistanceEdgePair& left, const DistanceEdgePair& right) const
        {
            return left.first > right.first;
        }
    };
    
    // Recover the hits we'll work on.
    // Note that we use on the skeleton hits so will need to recover them
    reco::HitPairListPtr& hitPairListPtr = clusterParameters.m_hitPairListPtr;
    reco::HitPairListPtr  skeletonListPtr;
    
    // We want to work with the "skeleton" hits so first step is to call the algorithm to
    // recover only these hits from the entire input collection
    m_skeletonAlg.GetSkeletonHits(hitPairListPtr, skeletonListPtr);
    
    // Skeleton hits are nice but we can do better if we then make a pass through to "average"
    // the skeleton hits position in the Y-Z plane
    m_skeletonAlg.AverageSkeletonPositions(skeletonListPtr);
    
    // First task is to define and build the adjacency map
    Hit3DToDistanceMap hit3DToDistanceMap;
    
    for(reco::HitPairListPtr::const_iterator hit3DOuterItr = skeletonListPtr.begin(); hit3DOuterItr != skeletonListPtr.end(); )
    {
        const reco::ClusterHit3D* hit3DOuter   = *hit3DOuterItr++;
        DistanceHit3DPairList&    outerHitList = hit3DToDistanceMap[hit3DOuter];
        TVector3                  outerPos(hit3DOuter->getPosition()[0], hit3DOuter->getPosition()[1], hit3DOuter->getPosition()[2]);
        
        for(reco::HitPairListPtr::const_iterator hit3DInnerItr = hit3DOuterItr; hit3DInnerItr != skeletonListPtr.end(); hit3DInnerItr++)
        {
            const reco::ClusterHit3D* hit3DInner = *hit3DInnerItr;
            TVector3                  innerPos(hit3DInner->getPosition()[0], hit3DInner->getPosition()[1], hit3DInner->getPosition()[2]);
            TVector3                  deltaPos = innerPos - outerPos;
            double                    hitDistance(deltaPos.Mag());
            
            if (hitDistance > 20.) continue;
            
            hit3DToDistanceMap[hit3DInner].emplace_back(DistanceHit3DPair(hitDistance,hit3DOuter));
            outerHitList.emplace_back(DistanceHit3DPair(hitDistance,hit3DInner));
        }
        
        // Make sure our membership bit is clear
        hit3DOuter->clearStatusBits(reco::ClusterHit3D::SELECTEDBYMST);
    }
    
    // Make pass through again to order each of the lists
    for(auto& mapPair : hit3DToDistanceMap)
    {
        mapPair.second.sort(Hit3DDistanceOrder());
    }
    
    // Get the containers for the MST to operate on/with
    Hit3DList            hit3DList;
    DistanceEdgePairList distanceEdgePairList;
    
    // Initialize with first element
    hit3DList.emplace_back(skeletonListPtr.front());
    distanceEdgePairList.emplace_back(DistanceEdgePair(0.,Hit3DEdgePair(hit3DList.begin(),hit3DList.begin())));
    
    skeletonListPtr.front()->setStatusBit(reco::ClusterHit3D::SELECTEDBYMST);
    
    double largestDistance(0.);
    double averageDistance(0.);
    
    // Now run the MST
    // Basically, we loop until the MST list is the same size as the input list
    while(hit3DList.size() < skeletonListPtr.size())
    {
        Hit3DList::iterator bestHit3DIter = hit3DList.begin();
        double              bestDist      = 10000000.;
        
        // Loop through all hits currently in the list and look for closest hit not in the list
        for(Hit3DList::iterator hit3DIter = hit3DList.begin(); hit3DIter != hit3DList.end(); hit3DIter++)
        {
            const reco::ClusterHit3D* hit3D = *hit3DIter;
            
            // For the given 3D hit, find the closest to it that is not already in the list
            DistanceHit3DPairList& nearestList = hit3DToDistanceMap[hit3D];
            
            while(!nearestList.empty())
            {
                const reco::ClusterHit3D* hit3DToCheck = nearestList.front().second;
                
                if (!hit3DToCheck->bitsAreSet(reco::ClusterHit3D::SELECTEDBYMST))
                {
                    if (nearestList.front().first < bestDist)
                    {
                        bestHit3DIter = hit3DIter;
                        bestDist      = nearestList.front().first;
                    }
                    
                    break;
                }
                else nearestList.pop_front();
            }
        }
        
        if (bestDist > largestDistance) largestDistance = bestDist;
        
        averageDistance += bestDist;
        
        // Now we add the best hit not in the list to our list, keep track of the distance
        // to the object it was closest to
        const reco::ClusterHit3D* bestHit3D = *bestHit3DIter;                               // "best" hit already in the list
        const reco::ClusterHit3D* nextHit3D = hit3DToDistanceMap[bestHit3D].front().second; // "next" hit we are adding to the list
        
        Hit3DList::iterator nextHit3DIter = hit3DList.insert(hit3DList.end(),nextHit3D);
        
        distanceEdgePairList.emplace_back(DistanceEdgePair(bestDist,Hit3DEdgePair(bestHit3DIter,nextHit3DIter)));
        
        nextHit3D->setStatusBit(reco::ClusterHit3D::SELECTEDBYMST);
    }
    
    averageDistance /= double(hit3DList.size());
    
    double thirdDist = 2.*sqrt(clusterParameters.m_skeletonPCA.getEigenValues()[2]);
    
    // Ok, find the largest distance in the iterator map
    distanceEdgePairList.sort(DistanceEdgePairOrder());
    
    DistanceEdgePairList::iterator largestDistIter = distanceEdgePairList.begin();
    
    for(DistanceEdgePairList::iterator edgeIter = distanceEdgePairList.begin(); edgeIter != distanceEdgePairList.end(); edgeIter++)
    {
        if (edgeIter->first < thirdDist) break;
        
        largestDistIter = edgeIter;
    }
    
    reco::HitPairListPtr::iterator breakIter = largestDistIter->second.second;
    reco::HitPairListPtr bestList;
    
    bestList.resize(std::distance(hit3DList.begin(), breakIter));
    
    std::copy(hit3DList.begin(), breakIter, bestList.begin());
    
    // Remove from the grand hit list and see what happens...
    // The pieces below are incomplete and were really for testing only.
    hitPairListPtr.sort();
    bestList.sort();
    
    reco::HitPairListPtr::iterator newListEnd =
    std::set_difference(hitPairListPtr.begin(), hitPairListPtr.end(),
                        bestList.begin(),       bestList.end(),
                        hitPairListPtr.begin() );
    
    hitPairListPtr.erase(newListEnd, hitPairListPtr.end());
    
    return;
}

class CopyIfInRange
{
public:
    CopyIfInRange(double maxRange) : m_maxRange(maxRange) {}
    
    bool operator()(const reco::ClusterHit3D* hit3D)
    {
        return hit3D->getDocaToAxis() < m_maxRange;
    }
private:
    double m_maxRange;
};

void Cluster3D::splitClustersWithHough(ClusterParameters&       clusterParameters,
                                       reco::HitPairClusterMap& hitPairClusterMap,
                                       ClusterParametersList&   clusterParametersList) const
{
    // @brief A method for splitted "crossed tracks" clusters into separate clusters
    //
    // If this routine is called then we believe we have a cluster which needs splitting.
    // The specific topology we are looking for is two long straight tracks which cross at some
    // point in close proximity so their hits were joined into a single 3D cluster. The method
    // to split this topology is to let the hough transform algorithm find the two leading candidates
    // and then to see if we use those to build two clusters instead of one.
    
    // Recover the hits we'll work on.
    // Note that we use on the skeleton hits so will need to recover them
    reco::HitPairListPtr& hitPairListPtr = clusterParameters.m_hitPairListPtr;
    reco::HitPairListPtr  skeletonListPtr;
    
    // We want to work with the "skeleton" hits so first step is to call the algorithm to
    // recover only these hits from the entire input collection
    m_skeletonAlg.GetSkeletonHits(hitPairListPtr, skeletonListPtr);
    
    // Skeleton hits are nice but we can do better if we then make a pass through to "average"
    // the skeleton hits position in the Y-Z plane
    m_skeletonAlg.AverageSkeletonPositions(skeletonListPtr);
    
    // Define the container for our lists of hits
    reco::HitPairListPtrList hitPairListPtrList;
    
    // Now feed this to the Hough Transform to find candidate straight lines
    m_seedFinderAlg.findTrackHits(skeletonListPtr, clusterParameters.m_skeletonPCA, hitPairListPtrList);
    
    // We need at least two lists or else there is nothing to do
    if (hitPairListPtrList.size() < 2) return;
    
    // The game plan will be the following:
    // 1) Take the first list of hits and run the PCA on this to get an axis
    //    - Then calculate the 3d doca for ALL hits in the cluster to this axis
    //    - Move all hits within "3 sigam" of the axis to a new list
    // 2) run the PCA on the second list of hits to get that axis
    //    - Then calculate the 3d doca for all hits in our first list
    //    - Copy hits in the first list which are within 3 sigma of the new axis
    //      back into our original cluster - these are shared hits
    reco::HitPairListPtrList::iterator hitPairListIter = hitPairListPtrList.begin();
    reco::HitPairListPtr&              firstHitList    = *hitPairListIter++;
    reco::PrincipalComponents          firstHitListPCA;
    
    m_pcaAlg.PCAAnalysis_3D(firstHitList, firstHitListPCA);
    
    // Make sure we have a successful calculation.
    if (firstHitListPCA.getSvdOK())
    {
        // The fill routines below will expect to see unused 2D hits so we need to clear the
        // status bits... and I am not sure of a better way...
        for(const auto& hit3D : hitPairListPtr)
        {
            for(const auto& hit2D : hit3D->getHits()) hit2D->clearStatusBits(0x1);
        }
        
        // Calculate the 3D doca's for the hits which were used to make this PCA
        m_pcaAlg.PCAAnalysis_calc3DDocas(firstHitList, firstHitListPCA);
        
        // Divine from the ether some maximum allowed range for transfering hits
        double allowedHitRange = 6. * firstHitListPCA.getAveHitDoca();
        
        // Now go through and calculate the 3D doca's for ALL the hits in the original cluster
        m_pcaAlg.PCAAnalysis_calc3DDocas(hitPairListPtr, firstHitListPCA);
        
        // Get a new container for these hits, assume worst case that everything transfers
        // We'll make this new cluster the last in the list
        int newClusterKey = hitPairClusterMap.rbegin()->first + 1;
        
        reco::HitPairListPtr& newClusterHitList = hitPairClusterMap[newClusterKey];
        
        newClusterHitList.resize(hitPairListPtr.size());
        
        // Do the actual copy of the hits we want
        reco::HitPairListPtr::iterator newListEnd =
            std::copy_if(hitPairListPtr.begin(), hitPairListPtr.end(), newClusterHitList.begin(), CopyIfInRange(allowedHitRange));
        
        // Shrink to fit
        newClusterHitList.resize(std::distance(newClusterHitList.begin(), newListEnd));
        
        // And now remove these hits from the original cluster
        hitPairListPtr.remove_if(CopyIfInRange(allowedHitRange));
        
        // Let's make a new cluster from this set of hits
        clusterParametersList.push_back(ClusterParameters(newClusterHitList));
        
        // Can we get a reference to what we just created?
        ClusterParameters& newClusterParams = clusterParametersList.back();

        // Now "fill" the cluster parameters but turn off the hit rejection
        FillClusterParams(newClusterParams, 0., 1.);

        // Set the skeleton pca to the value calculated above on input
        clusterParameters.m_skeletonPCA = firstHitListPCA;

        // We are done with splitting out one track. Because the two tracks cross in
        // close proximity, this is the one case where we might consider sharing 3D hits
        // So let's make a little detour here to try to copy some of those hits back into
        // the main hit list
        reco::HitPairListPtr&     secondHitList = *hitPairListIter;
        reco::PrincipalComponents secondHitListPCA;
        
        m_pcaAlg.PCAAnalysis_3D(secondHitList, secondHitListPCA);
        
        // Make sure we have a successful calculation.
        if (secondHitListPCA.getSvdOK())
        {
            // Calculate the 3D doca's for the hits which were used to make this PCA
            m_pcaAlg.PCAAnalysis_calc3DDocas(secondHitList, secondHitListPCA);
            
            // Since this is the "other" cluster, we'll be a bit more generous in adding back hits
            double newAllowedHitRange = 6. * secondHitListPCA.getAveHitDoca();
            
            // Go through and calculate the 3D doca's for the hits in our new candidate cluster
            m_pcaAlg.PCAAnalysis_calc3DDocas(newClusterHitList, secondHitListPCA);
            
            // Create a temporary list to fill with the hits we might want to save
            reco::HitPairListPtr tempHitList(newClusterHitList.size());
            
            // Do the actual copy of the hits we want...
            reco::HitPairListPtr::iterator tempListEnd =
                std::copy_if(newClusterHitList.begin(), newClusterHitList.end(), tempHitList.begin(), CopyIfInRange(newAllowedHitRange));

            hitPairListPtr.insert(hitPairListPtr.end(), tempHitList.begin(), tempListEnd);
        }
 
        // Of course, now we need to modify the original cluster parameters
        ClusterParameters originalParams(hitPairListPtr);
        
        // Now "fill" the cluster parameters but turn off the hit rejection
        FillClusterParams(originalParams, 0., 1.);

        // Overwrite original cluster parameters with our new values
        clusterParameters.m_clusterParams = originalParams.m_clusterParams;
        clusterParameters.m_fullPCA       = originalParams.m_fullPCA;
        clusterParameters.m_skeletonPCA   = secondHitListPCA;
    }

    return;
}
    
void Cluster3D::ProduceArtClusters(art::Event&              evt,
                                   HitPairList&             hitPairVector,
                                   reco::HitPairClusterMap& hitPairClusterMap,
                                   ClusterParametersList&   clusterParametersList,
                                   RecobHitToPtrMap&        hitToPtrMap) const
{
    /**
     *  @brief The workhorse to take the candidate 3D clusters and produce all of the necessary art output
     */
    
    mf::LogDebug("Cluster3D") << " *** Cluster3D::ProduceArtClusters() *** " << std::endl;
    
    std::unique_ptr< std::vector<recob::PCAxis> >     artPCAxisVector( new std::vector<recob::PCAxis>         );
    std::unique_ptr< std::vector<recob::PFParticle> > artPFParticleVector( new std::vector<recob::PFParticle> );
    std::unique_ptr< std::vector<recob::Cluster> >    artClusterVector( new std::vector<recob::Cluster>       );
    std::unique_ptr< std::vector<recob::SpacePoint> > artSpacePointVector( new std::vector<recob::SpacePoint> );
    std::unique_ptr< std::vector<recob::Seed> >       artSeedVector( new std::vector<recob::Seed>             );
    
    std::unique_ptr< art::Assns<recob::Cluster,    recob::Hit> >         artClusterAssociations(    new art::Assns<recob::Cluster,    recob::Hit>          );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::PCAxis> >      artPFPartAxisAssociations( new art::Assns<recob::PFParticle, recob::PCAxis>       );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> >     artPFPartClusAssociations( new art::Assns<recob::PFParticle, recob::Cluster>      );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> >  artPFPartSPAssociations(   new art::Assns<recob::PFParticle, recob::SpacePoint>   );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Seed> >        artPFPartSeedAssociations( new art::Assns<recob::PFParticle, recob::Seed>         );
    std::unique_ptr< art::Assns<recob::Seed,       recob::Hit> >         artSeedHitAssociations(    new art::Assns<recob::Seed,       recob::Hit>         );
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> >         artSPHitAssociations(      new art::Assns<recob::SpacePoint, recob::Hit>          );
    
    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here, except that we override selected items
    // (so, thanks to metaprogramming, we finally have wrappers of wrappers);
    // configuration would happen here, but we are using the default
    // configuration for that algorithm
    using OverriddenClusterParamsAlg_t
      = cluster::OverriddenClusterParamsAlg<cluster::StandardClusterParamsAlg>;
    cluster::ClusterParamsImportWrapper<OverriddenClusterParamsAlg_t>
      ClusterParamAlgo;
    
    // Create id for space points
    int    spacePointID(0);
    int    pcaAxisID(0);
    
    // Make sure there is something to do here!
    if (!clusterParametersList.empty())
    {
        // Grab an instance of the Geometry as it is needed for the association of hits

        // indices for the clusters created
        int    clusterIdx(0);
        size_t pfParticleIdx(0);
        
        // Keep track of the position of the latest added 2D clusters
        size_t clusterStart(0);
        size_t clusterEnd(0);
        
        // Create id for space points
        int    spacePointID(0);
        
        // This is the loop over candidate 3D clusters
        // Note that it might be that the list of candidate clusters is modified by splitting
        // So we use the following construct to make sure we get all of them
        ClusterParametersList::iterator clusterParametersListItr = clusterParametersList.begin();
        
        while(clusterParametersListItr != clusterParametersList.end())
        {
            // Dereference to get the cluster paramters
            ClusterParameters& clusterParameters = *clusterParametersListItr;
            
            // It should be straightforward at this point to transfer information from our vector of clusters
            // to the larsoft objects... of course we still have some work to do first, in particular to
            // find the candidate seeds and their seed hits
            
            // We keep track of 2 PCA axes, the first is the "full" PCA run over all the 3D hits in the
            // candidate cluster. The second will be that derived from just using the "skeleton" hits.
            // Make a copy of the full PCA to keep that, then get a reference for the skeleton PCA
            reco::PrincipalComponents& fullPCA     = clusterParameters.m_fullPCA;
            reco::PrincipalComponents& skeletonPCA = clusterParameters.m_skeletonPCA;
            
            // The chances of getting here and this condition not being true are probably zero... but check anyway
            if (!fullPCA.getSvdOK())
            {
                mf::LogDebug("Cluster3D") << "--> no feature extraction done on this cluster!!" << std::endl;
                clusterParametersListItr++;
                continue;
            }
            
            // As tracks become more parallel to the wire plane the number of "ambiguous" 3D hits can increase
            // rapidly. Now that we have more information we can go back through these hits and do a better job
            // selecting "the right ones". Here we call the "medial skeleton" algorithm which uses a modification
            // of a standard medial skeleton procedure to get the 3D hits we want
            // But note that even this is hopeless in the worst case and, in fact, it can be a time waster
            // So bypass when you recognize that condition
            if (!aParallelHitsCluster(fullPCA))
            {
                int nSkeletonPoints = m_skeletonAlg.FindMedialSkeleton(clusterParameters.m_hitPairListPtr);
            
                // If enough skeleton points then rerun pca with only those
                if (nSkeletonPoints > 10)
                {
                    // Now rerun the principal components axis on just those points
                    m_pcaAlg.PCAAnalysis_3D(clusterParameters.m_hitPairListPtr, skeletonPCA, true);
                
                    // If there was a failure (can that happen?) then restore the full PCA
                    if (!skeletonPCA.getSvdOK()) skeletonPCA = fullPCA;
                }
            
                // Here we can try to handle a specific case. It can happen that two tracks (think CR muons here) pass so
                // close together at some point to get merged into one cluster. Now that we have skeletonized the hits and
                // have run the PCA on the skeleton points we can try to divide these two tracks. The signature will be that
                // their are a large number of total hits, that the PCA will have a large spread in two dimensions. The
                // spread in the third dimension will be an indicator of the actual separation between the two tracks
                // which we might try to exploit in the actual algorithm.
                // hardwire for now to see what is going on...
                if (skeletonPCA.getNumHitsUsed() > 1000 && skeletonPCA.getEigenValues()[1] > 100. && fabs(skeletonPCA.getEigenVectors()[2][0]) < m_parallelHitsCosAng)
                {
                    mf::LogDebug("Cluster3D") << "--> Detected crossed axes!! Total # hits: " << fullPCA.getNumHitsUsed() <<
                        "\n    Skeleton PCA # hits: " << skeletonPCA.getNumHitsUsed() << ", eigenValues: " <<
                        skeletonPCA.getEigenValues()[0] << ", " <<skeletonPCA.getEigenValues()[1] << ", " <<skeletonPCA.getEigenValues()[2] << std::endl;
                    
                    splitClustersWithHough(clusterParameters, hitPairClusterMap, clusterParametersList);
                }
            }
        
            // Start loop over views to build out the hit lists and the 2D cluster objects
            for(ViewToClusterParamsMap::const_iterator viewItr = clusterParameters.m_clusterParams.begin(); viewItr != clusterParameters.m_clusterParams.end(); viewItr++)
            {
                const RecobClusterParameters& clusParams = viewItr->second;
                
                // We love looping. In this case, our list of hits is comprised of "ClusterHits" and we need to get a RecobHitVector instead...
                RecobHitVector recobHits;
                
                for(HitVectorConst::const_iterator hitItr = clusParams.m_hitVector.begin(); hitItr != clusParams.m_hitVector.end(); hitItr++)
                {
                    art::Ptr<recob::Hit> hitPtr = hitToPtrMap[&(*hitItr)->getHit()];
                    recobHits.push_back(hitPtr);
                }
                
                // And sorting! Sorting is good for the mind, soul and body
                // ooopsss... don't do this else event display will look funky
//                std::sort(recobHits.begin(), recobHits.end());
                
                // Get the tdc/wire slope... from the unit vector...
                double dTdW(0.);
                double startWire(clusParams.m_startWire);
                double endWire(clusParams.m_endWire);
                double startTime(clusParams.m_startTime);
                double endTime(clusParams.m_endTime);
                double wirePitch(m_geometry->WirePitch(clusParams.m_view));
                double midWire(0);

                // Get wire number corresponding to the current position
                try
                {
                    midWire = 1.*m_geometry->NearestWire(skeletonPCA.getAvePosition(), clusParams.m_view);
                } catch (cet::exception& e)
                {
                    mf::LogWarning("Cluster3D") << "Exception caught finding nearest wire, position - " << e.what() << std::endl;
                    midWire  = skeletonPCA.getAvePosition()[2];
                    midWire /= wirePitch;
                }
                
                // Now sort out the slope in this view
                // This follows the code in the display pacakage
                double thetaWire     = m_geometry->Plane(clusParams.m_view).Wire(0).ThetaZ();
                double driftVelocity = m_detector->DriftVelocity();
                double timeTick      = m_detector->SamplingRate()*1.e-3;
                
                //rotate coord system CCW around x-axis by pi-thetawire
                //   new yprime direction is perpendicular to the wire direction
                //   in the same plane as the wires and in the direction of
                //   increasing wire number
                //use yprime-component of dir cos in rotated coord sys to get
                //   dTdW (number of time ticks per unit of wire pitch)
                double rotAng = 3.1416-thetaWire;
                double yPrime = std::cos(rotAng)*skeletonPCA.getEigenVectors()[0][1]+std::sin(rotAng)*skeletonPCA.getEigenVectors()[0][2];
                
                if (fabs(yPrime) < 0.0000001) yPrime = 0.0000001;
                
                dTdW = skeletonPCA.getEigenVectors()[0][0]*wirePitch/(driftVelocity*timeTick*yPrime);
                
                // Finally, use the time corresponding to this position
                double midTime = m_detector->ConvertXToTicks(skeletonPCA.getAvePosition()[0], clusParams.m_view, 0, 0);
                
                // plane ID is not a part of clusParams... get the one from the first hit
                geo::PlaneID plane; // invalid by default
                if (!recobHits.empty())
                  plane = recobHits.front()->WireID().planeID();
                
                // Ok, now adjust everything to draw a nice long line through our cluster
                double deltaWires = 0.5 * (endWire - startWire);
                
                startWire = midWire - deltaWires;
                endWire   = midWire + deltaWires;
                startTime = midTime - dTdW * deltaWires;
                endTime   = midTime + dTdW * deltaWires;
                
                // feed the algorithm with all the cluster hits
                ClusterParamAlgo.ImportHits(recobHits);
                
                // override the end angles: instead of using the standard
                // algorithm, we use a precomputed value
                double tan_dTdW(std::tan(dTdW));
                
                ClusterParamAlgo.OverrideParameter
                  (OverriddenClusterParamsAlg_t::cpStartAngle, tan_dTdW);
                ClusterParamAlgo.OverrideParameter
                  (OverriddenClusterParamsAlg_t::cpEndAngle, tan_dTdW);
                
                // create the recob::Cluster directly in the vector
                cluster::ClusterCreator artCluster(
                  ClusterParamAlgo,                     // algo
                  startWire,                            // start_wire
                  0.,                                   // sigma_start_wire
                  startTime,                            // start_tick
                  clusParams.m_sigmaStartTime,          // sigma_start_tick
                  endWire,                              // end_wire
                  0.,                                   // sigma_end_wire,
                  endTime,                              // end_tick
                  clusParams.m_sigmaEndTime,            // sigma_end_tick
                  clusterIdx++,                         // ID
                  clusParams.m_view,                    // view
                  plane,                                // plane
                  recob::Cluster::Sentry                // sentry
                  );
                
                artClusterVector->emplace_back(artCluster.move());
                             
                util::CreateAssn(*this, evt, *artClusterVector, recobHits, *artClusterAssociations);
                clusterEnd++;
            }
            
            // Deal with converting the Hit Pairs to art
            // Recover the hit pairs and start looping! Love to loop!
            reco::HitPairListPtr& clusHitPairVector = clusterParameters.m_hitPairListPtr;
            
            // Last, let's try to get seeds for tracking..
            // Keep track of how many we have so far
            size_t numSeedsStart = artSeedVector->size();
            
            // Call the magical algorith to do the dirty work
            findTrackSeeds(evt, clusterParameters, hitToPtrMap, *artSeedVector, *artSeedHitAssociations);
            
            // Right now error matrix is uniform...
            double spError[] = {1., 0., 1., 0., 0., 1.};
            
            // Keep track of current start for space points
            int spacePointStart(spacePointID);
            
            // Copy these hits to the vector to be stored with the event
            for (auto& hitPair : clusHitPairVector)
            {
                // Don't make space point if this hit was "rejected"
                if (hitPair->bitsAreSet(reco::ClusterHit3D::REJECTEDHIT)) continue;
                
                double chisq = 1.;    // secret handshake...
                
                if      ( hitPair->bitsAreSet(reco::ClusterHit3D::SKELETONHIT) && !hitPair->bitsAreSet(reco::ClusterHit3D::EDGEHIT)) chisq = -1.;  // pure skeleton point
                else if (!hitPair->bitsAreSet(reco::ClusterHit3D::SKELETONHIT) &&  hitPair->bitsAreSet(reco::ClusterHit3D::EDGEHIT)) chisq = -2.;  // pure edge point
                else if ( hitPair->bitsAreSet(reco::ClusterHit3D::SKELETONHIT) &&  hitPair->bitsAreSet(reco::ClusterHit3D::EDGEHIT)) chisq = -3.;  // skeleton and edge point

                if      (hitPair->bitsAreSet(reco::ClusterHit3D::SEEDHIT)                                                          ) chisq = -4.;  // Seed point
                
                if (hitPair->getHits().size() < 3) chisq = -10.;
                
                // Mark this hit pair as in use
                hitPair->setStatusBit(reco::ClusterHit3D::MADESPACEPOINT);
                double spacePointPos[] = {hitPair->getPosition()[0],hitPair->getPosition()[1],hitPair->getPosition()[2]};
                artSpacePointVector->push_back(recob::SpacePoint(spacePointPos, spError, chisq, spacePointID++));
            }
            
            // Empty daughter vector for now
            std::vector<size_t> nullVector;
            
            // Create the PFParticle to tie the pieces together
            size_t parentID(recob::PFParticle::kPFParticlePrimary);
            
            recob::PFParticle pfParticle(13, pfParticleIdx++, parentID, nullVector);
            artPFParticleVector->push_back(pfParticle);
            
            // Look at making the PCAxis and associations - for both the skeleton (the first) and the full
            recob::PCAxis skelPcAxis(skeletonPCA.getSvdOK(),
                                     skeletonPCA.getNumHitsUsed(),
                                     skeletonPCA.getEigenValues(),
                                     skeletonPCA.getEigenVectors(),
                                     skeletonPCA.getAvePosition(),
                                     skeletonPCA.getAveHitDoca(),
                                     pcaAxisID++);
            
            artPCAxisVector->push_back(skelPcAxis);
            
            recob::PCAxis fullPcAxis(fullPCA.getSvdOK(),
                                     fullPCA.getNumHitsUsed(),
                                     fullPCA.getEigenValues(),
                                     fullPCA.getEigenVectors(),
                                     fullPCA.getAvePosition(),
                                     fullPCA.getAveHitDoca(),
                                     pcaAxisID++);
            
            artPCAxisVector->push_back(fullPcAxis);
            
            util::CreateAssn(*this, evt, *artPFParticleVector, *artPCAxisVector, *artPFPartAxisAssociations, artPCAxisVector->size()-2, artPCAxisVector->size());
            
            // Create associations to the PFParticle
            util::CreateAssn(*this, evt, *artPFParticleVector, *artSeedVector, *artPFPartSeedAssociations, numSeedsStart, artSeedVector->size());
            
            // Make associations to the 2D cluster objects
            util::CreateAssn(*this, evt, *artPFParticleVector, *artClusterVector, *artPFPartClusAssociations, clusterStart, clusterEnd);
            
            // Make associations to the SpacePoints
            util::CreateAssn(*this, evt, *artPFParticleVector, *artSpacePointVector, *artPFPartSPAssociations, spacePointStart, spacePointID);
            
            // Update the start/end indices
            clusterStart = clusterEnd;
            
            // Go to next cluster parameters object
            clusterParametersListItr++;
        }
    }
    
    // Right now error matrix is uniform...
    double spError[] = {1., 0., 1., 0., 0., 1.};
    
    // Run through the HitPairVector and add any unused hit pairs to the list
    for(auto& hitPair : hitPairVector)
    {
        if (!hitPair->bitsAreSet(reco::ClusterHit3D::MADESPACEPOINT))
        {
            double spacePointPos[] = {hitPair->getPosition()[0],hitPair->getPosition()[1],hitPair->getPosition()[2]};
            artSpacePointVector->push_back(recob::SpacePoint(spacePointPos, spError, 1., spacePointID++));
        }
    }
    
    // Finaly done, now output everything to art
    evt.put(std::move(artPCAxisVector));
    evt.put(std::move(artPFParticleVector));
    evt.put(std::move(artClusterVector));
    evt.put(std::move(artSpacePointVector));
    evt.put(std::move(artSeedVector));
    evt.put(std::move(artPFPartAxisAssociations));
    evt.put(std::move(artPFPartClusAssociations));
    evt.put(std::move(artClusterAssociations));
    evt.put(std::move(artPFPartSPAssociations));
    evt.put(std::move(artPFPartSeedAssociations));
    evt.put(std::move(artSeedHitAssociations));
    evt.put(std::move(artSPHitAssociations));
    
    return;
}

} // namespace lar_cluster3d
