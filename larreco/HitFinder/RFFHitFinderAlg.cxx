/*!
 * Title:   RFFHitFinderAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Class that runs the RFF HitFinder. Implements an RFFHitFitter, and takes 
 * the result and stores it in recob::Hit objects. 
 *
 * Input:  recob::Wire
 * Output: recob::Hit
*/

#include "RFFHitFinderAlg.h"
#include <numeric>
#include <functional>
#include "/usr/local/include/zmq.h"

hit::RFFHitFinderAlg::RFFHitFinderAlg(fhicl::ParameterSet const& p)
{
  fMatchThresholdVec = p.get< std::vector<float> >("MeanMatchThreshold");
  fMergeMultiplicityVec = p.get< std::vector<unsigned int> >("MinMergeMultiplicity");
  fAmpThresholdVec = p.get< std::vector<float> >("AmplitudeThreshold",std::vector<float>(1,0.0));

  fNWorkers = p.get<size_t>("NWorkers",1);

  fContext = zmq_ctx_new();

  //establish the sockets...

  fReceiver_Hits = zmq_socket(fContext, ZMQ_PULL);
  zmq_bind (fReceiver_Hits, "inproc://hits");

  fReceiver_fin = zmq_socket(fContext, ZMQ_PULL);
  zmq_bind (fReceiver_fin, "inproc://finished");
  
  fController = zmq_socket (fContext, ZMQ_PUB);
  zmq_bind (fController, "inproc://control");

  fPairNROI =  zmq_socket(fContext, ZMQ_PAIR);
  int linger_time=0;
  zmq_setsockopt(fPairNROI,ZMQ_LINGER,&linger_time,sizeof(int));
  zmq_bind (fPairNROI, "inproc://nrois");
}

hit::RFFHitFinderAlg::~RFFHitFinderAlg()
{
  
  zmq_close(fPairNROI);
  zmq_close(fReceiver_Hits);
  zmq_close(fReceiver_fin);
  zmq_close(fController);
  zmq_ctx_destroy(fContext);
}

void hit::RFFHitFinderAlg::SetFitterParamsVectors(geo::Geometry const& geo)
{
  const unsigned int n_planes = geo.Nplanes();

  //If size zero, throw. If size one, assume same for all planes.
  //If size > 1 but < n_planes, throw. If size = n_plane, good.

  if(fMatchThresholdVec.size()==0 || 
     fMergeMultiplicityVec.size()==0 ||
     fAmpThresholdVec.size()==0)
    throw std::runtime_error("Error in RFFHitFinderAlg: Configured with zero planes.");

  if( (fMatchThresholdVec.size()>1 && fMatchThresholdVec.size()<n_planes) ||
      (fMergeMultiplicityVec.size()>1 && fMergeMultiplicityVec.size()<n_planes) ||
      (fAmpThresholdVec.size()>1 && fAmpThresholdVec.size()<n_planes) )
    throw std::runtime_error("Error in RFFHitFinderAlg: Configured with incorrect n_planes.");
    
  if(fMatchThresholdVec.size()==1)
    fMatchThresholdVec.resize(n_planes,fMatchThresholdVec[0]);

  if(fMergeMultiplicityVec.size()==1)
    fMergeMultiplicityVec.resize(n_planes,fMergeMultiplicityVec[0]);

  if(fAmpThresholdVec.size()==1)
    fAmpThresholdVec.resize(n_planes,fAmpThresholdVec[0]);
}
/*
void hit::RFFHitFinderAlg::SetFitterParams(unsigned int p)
{
  fFitter.SetFitterParams(fMatchThresholdVec[p],fMergeMultiplicityVec[p],fAmpThresholdVec[p]);
}
*/

void hit::RFFHitFinderAlg::SendDataToWorkers(std::vector<recob::Wire> const& wireVector,
					     geo::Geometry const& geo,
					     void* context,
					     std::vector<float> const& vec_match_thresh,
					     std::vector<unsigned int> const& vec_merge_multi,
					     std::vector<float> const& vec_amp_thresh)
{

  void *sender_ROI = zmq_socket(context, ZMQ_PUSH);
  zmq_bind (sender_ROI, "inproc://workers");

  void* n_ROIs_sock = zmq_socket(context, ZMQ_PAIR);
  zmq_connect (n_ROIs_sock, "inproc://nrois");

  size_t total_rois=0;
  
  for(auto const& wire : wireVector){
    
    //printf("\tRunning channel %u\n",wire.Channel());
    
    geo::SigType_t const& sigtype = geo.SignalType(wire.Channel());
    geo::WireID const& wireID = geo.ChannelToWire(wire.Channel()).at(0);
    for(auto const& roi : wire.SignalROI().get_ranges()){
      
      FitterInput_t f_input;
      f_input.sigtype = sigtype;//geo.SignalType(wire.Channel());
      f_input.wireid  = wireID;//geo.ChannelToWire(wire.Channel()).at(0);
      f_input.view    = wire.View();
      f_input.channel = wire.Channel();
      f_input.roi_start = roi.begin_index();
      f_input.roi_end = roi.begin_index()+roi.size();
      f_input.roi_dataptr = &(roi.data());
      f_input.match_thresh = vec_match_thresh[f_input.view];
      f_input.merge_multi = vec_merge_multi[f_input.view];
      f_input.amp_thresh = vec_amp_thresh[f_input.view];
      
      zmq_msg_t msg_roi;
      zmq_msg_init_size(&msg_roi,sizeof(FitterInput_t));
      memcpy(zmq_msg_data(&msg_roi),&f_input,sizeof(FitterInput_t));
      //zmq_msg_init_data(&msg_roi,&f_input,sizeof(FitterInput_t),NULL,NULL);
      zmq_sendmsg(sender_ROI,&msg_roi,0);
      ++total_rois;
      ////printf("\t\tSent roi ...");
      
    }//end loop over ROIs on wire
    
    //printf("\tSent %lu rois ...\n",total_rois);
    
  }//end loop over wires

  zmq_msg_t msg_nroi;
  zmq_msg_init_size(&msg_nroi,sizeof(size_t));
  memcpy(zmq_msg_data(&msg_nroi),&total_rois,sizeof(size_t));
  zmq_sendmsg(n_ROIs_sock,&msg_nroi,0);
  
  zmq_close(sender_ROI);
  zmq_close(n_ROIs_sock);


  //printf("DataSender done.\n");
}

void hit::RFFHitFinderAlg::Run(std::vector<recob::Wire> const& wireVector,
			       std::vector<recob::Hit>& hitVector,
			       geo::Geometry const& geo)
{

  //printf("Running RFF Hit Finder Alg\n");

  std::thread data_sender(SendDataToWorkers,std::ref(wireVector),std::ref(geo),fContext,
			  std::ref(fMatchThresholdVec),std::ref(fMergeMultiplicityVec),std::ref(fAmpThresholdVec));
  
  std::vector<std::thread> workers;
  for(size_t i_w=0; i_w<fNWorkers; ++i_w)
    workers.emplace_back(FitterWorker,fContext);

  size_t total_rois=999999;

  hitVector.reserve(5*wireVector.size());
  
  //now we need to collect the hits. Exit only when we verify that our vector contains all hits we expect
  size_t rois_received=0;
  size_t expected_hits=0;
  while (rois_received<total_rois || hitVector.size() < expected_hits) {
    zmq_pollitem_t items [] = {
      { fReceiver_Hits, 0, ZMQ_POLLIN, 0 },
      { fReceiver_fin, 0, ZMQ_POLLIN, 0 },
      { fPairNROI, 0, ZMQ_POLLIN, 0 }      
    };
    zmq_poll (items, 3, -1);
    if (items [0].revents & ZMQ_POLLIN) {

      zmq_msg_t msg_recv;
      zmq_msg_init(&msg_recv);
      zmq_recvmsg(fReceiver_Hits,&msg_recv,0);

      hitVector.emplace_back(*(recob::Hit*)(zmq_msg_data(&msg_recv)));

      //printf("\t\t\tHit received! Channel %u\n",hitVector.back().Channel());
      
    }
    if (items [1].revents & ZMQ_POLLIN){

      zmq_msg_t msg_recv;
      zmq_msg_init(&msg_recv);
      zmq_recvmsg(fReceiver_fin,&msg_recv,0);

      expected_hits += *(size_t*)(zmq_msg_data(&msg_recv));
      ++rois_received;
      //printf("\tReceived last hit for channel! %lu / %lu rois recevied, %lu / %lu hits received\n",rois_received,total_rois,hitVector.size(),expected_hits); 
    }
    if(items[2].revents & ZMQ_POLLIN){
      zmq_msg_t msg_nroi;
      zmq_msg_init(&msg_nroi);
      zmq_recvmsg(fPairNROI,&msg_nroi,0);
      
      memcpy(&total_rois,zmq_msg_data(&msg_nroi),sizeof(size_t));
    }
  }

  //printf("FINISHED! Gonna kill our workers now.\n");
  
  //upon exit, we need to kill the workers
  zmq_send(fController, "KILL",4,0);

  data_sender.detach();
  for(size_t i_w=0; i_w<fNWorkers; ++i_w)
    workers[i_w].join();
  
}

void hit::RFFHitFinderAlg::FitterWorker(void* context){

  void *receiver = zmq_socket(context, ZMQ_PULL);
  zmq_connect (receiver, "inproc://workers");

  void *sender_hits = zmq_socket(context, ZMQ_PUSH);
  zmq_connect (sender_hits, "inproc://hits");

  void *sender_fin = zmq_socket(context, ZMQ_PUSH);
  zmq_connect (sender_fin, "inproc://finished");
  
  void *controller = zmq_socket (context, ZMQ_SUB);
  zmq_connect (controller, "inproc://control");
  zmq_setsockopt (controller, ZMQ_SUBSCRIBE, "", 0);

  //printf("Worker started\n");

  
  while (1) {
    zmq_pollitem_t items [] = {
      { receiver, 0, ZMQ_POLLIN, 0 },
      { controller, 0, ZMQ_POLLIN, 0 }
    };
    zmq_poll (items, 2, -1);
    if (items [0].revents & ZMQ_POLLIN) {

      //printf("Worker received data!\n");

      zmq_msg_t msg_recv;
      zmq_msg_init(&msg_recv);
      zmq_recvmsg(receiver,&msg_recv,0);

      FitterInput_t f_input;
      memcpy(&f_input,zmq_msg_data(&msg_recv),zmq_msg_size(&msg_recv));

      zmq_msg_close(&msg_recv);

      //printf("Received channe %u, roi_start %u\n",f_input.channel,f_input.roi_start);

      
      RFFHitFitter worker_fitter(f_input.match_thresh,
				 f_input.merge_multi,
				 f_input.amp_thresh);
      worker_fitter.RunFitter(*(f_input.roi_dataptr));

      std::vector<recob::Hit> hitVector;      
      size_t n_hits = EmplaceHit(hitVector,
				 worker_fitter,
				 std::accumulate(f_input.roi_dataptr->begin(),
						 f_input.roi_dataptr->end(),
						 0.0),
				 f_input.roi_start,
				 f_input.roi_end,
				 f_input.sigtype,
				 f_input.wireid,
				 f_input.channel,
				 f_input.view);

      //printf("\tFitter results, channel %u: %lu hits\n",f_input.channel,n_hits);
      
      for(size_t i_hit=0; i_hit<hitVector.size(); ++i_hit){
	zmq_msg_t msg_send;
	zmq_msg_init_size(&msg_send,sizeof(recob::Hit));
	memcpy(zmq_msg_data(&msg_send),&(hitVector[i_hit]),sizeof(recob::Hit));
	zmq_sendmsg(sender_hits,&msg_send,0);
      }
      zmq_msg_t msg_send;
      zmq_msg_init_size(&msg_send,sizeof(size_t));
      memcpy(zmq_msg_data(&msg_send),&n_hits,sizeof(size_t));
      zmq_sendmsg(sender_fin,&msg_send,0);
      //printf("\tFinal fitter results sent! channel %u: %lu hits\n",f_input.channel,n_hits);
      
    }
    //  Any waiting controller command acts as 'KILL'
    if (items [1].revents & ZMQ_POLLIN){
      //printf("Thread got KILL message.\n");
      break;                      //  Exit loop
    }
  }


  
  zmq_close(receiver);
  zmq_close(sender_hits);
  zmq_close(sender_fin);
  zmq_close(controller);

  //printf("WorkerThread done.\n");
}

size_t hit::RFFHitFinderAlg::EmplaceHit(std::vector<recob::Hit>& hitVector,
					RFFHitFitter & worker_fitter,
					float const& summedADCTotal,
					raw::TDCtick_t const& startTick, raw::TDCtick_t const& endTick,
					geo::SigType_t const& sigtype, geo::WireID const& wireID,
					raw::ChannelID_t const& channel, geo::View_t const& view)
{

  float totalArea = 0.0;
  std::vector<float> areaVector(worker_fitter.NHits());
  std::vector<float> areaErrorVector(worker_fitter.NHits());
  std::vector<float> areaFracVector(worker_fitter.NHits());

  for(size_t ihit=0; ihit < worker_fitter.NHits(); ihit++){
    areaVector[ihit] = worker_fitter.AmplitudeVector()[ihit]*worker_fitter.SigmaVector()[ihit]*SQRT_TWO_PI;
    areaErrorVector[ihit] =
      SQRT_TWO_PI*std::sqrt(worker_fitter.AmplitudeVector()[ihit]*worker_fitter.SigmaErrorVector()[ihit]*worker_fitter.AmplitudeVector()[ihit]*worker_fitter.SigmaErrorVector()[ihit] +
			    worker_fitter.AmplitudeErrorVector()[ihit]*worker_fitter.SigmaVector()[ihit]*worker_fitter.AmplitudeErrorVector()[ihit]*worker_fitter.SigmaVector()[ihit]);
    totalArea += areaVector[ihit];
  }

  for(size_t ihit=0; ihit < worker_fitter.NHits(); ihit++){
    areaFracVector[ihit] = areaVector[ihit]/totalArea;
    
    hitVector.emplace_back(channel,
			   startTick,
			   endTick,
			   worker_fitter.MeanVector()[ihit]+(float)startTick,
			   worker_fitter.MeanErrorVector()[ihit],
			   worker_fitter.SigmaVector()[ihit],
			   worker_fitter.AmplitudeVector()[ihit],
			   worker_fitter.AmplitudeErrorVector()[ihit],
			   summedADCTotal*areaFracVector[ihit],
			   areaVector[ihit],
			   areaErrorVector[ihit],
			   worker_fitter.NHits(),
			   ihit,
			   -999.,
			   -999,
			   view,
			   sigtype,
			   wireID);
  }

  return worker_fitter.NHits();
}
