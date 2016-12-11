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
#include <zmq.hpp>

hit::RFFHitFinderAlg::RFFHitFinderAlg(fhicl::ParameterSet const& p)
{
  fMatchThresholdVec = p.get< std::vector<float> >("MeanMatchThreshold");
  fMergeMultiplicityVec = p.get< std::vector<unsigned int> >("MinMergeMultiplicity");
  fAmpThresholdVec = p.get< std::vector<float> >("AmplitudeThreshold",std::vector<float>(1,0.0));

  fNWorkers = p.get<size_t>("NWorkers",1);
  
  fContext = zmq_new_ctx();

  //establish the sockets...

  fSender_ROI = zmq_socket(context, ZMQ_PUSH);
  zmq_bind (fSender_ROI, "inproc://workers");

  fReceiver_Hits = zmq_socket(context, ZMQ_PULL);
  zmq_bind (fReceiver_Hits, "inproc://hits");

  fReceiver_fin = zmq_socket(context, ZMQ_PULL);
  zmq_bind (fReceiver_fin, "inproc://finished");
  
  fController = zmq_socket (context, ZMQ_PUB);
  zmq_bind (fController, "inproc://control");

}

hit::RFFHitFinderAlg::~RFFHitFinderAlg()
{
  zmq_close(fSender_ROI);
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

void hit::RFFHitFinderAlg::SetFitterParams(unsigned int p)
{
  fFitter.SetFitterParams(fMatchThresholdVec[p],fMergeMultiplicityVec[p],fAmpThresholdVec[p]);
}

void hit::RFFHitFinderAlg::Run(std::vector<recob::Wire> const& wireVector,
			       std::vector<recob::Hit>& hitVector,
			       geo::Geometry const& geo)
{
  for(size_t i_w=0; i_w<fNWorkers; ++i_w)
    fWorkers.emplace_back(FitterWorker);

  size_t total_rois = 0;
  
  hitVector.reserve(wireVector.size());
  for(auto const& wire : wireVector){

    geo::SigType_t const& sigtype = geo.SignalType(wire.Channel());
    geo::WireID const& wireID = geo.ChannelToWire(wire.Channel()).at(0);

    /*
    zmq_msg_t* refmsg_sigtype,refmsg_wireID,refmsg_view;
    
    zmq_msg_init_size(msg_view,sizeof(unsigned int));
    zmq_msg_init_size(msg_sigtype,sizeof(geo::SigType_t));
    zmq_msg_init_size(msg_wireID,sizeof(geo::WireID));

    memcpy(msg_view->data(),wire.View(),sizeof(unsigned int));
    memcpy(msg_sigtype->data(),geo.SignalType(wire.Channel()),sizeof(geo::SigType_t));
    memcpy(msg_wireID->data(),geo.ChannelToWire(wire.Channel()).at(0),sizeof(geo::WireID));
    */
        
    //SetFitterParams(wire.View());
    //for(size_t i_roi=0; i_roi<wire.SignalROI().n_ranges(); ++i_roi){
    for(auto const& roi : wire.SignalROI().get_ranges()){

      FitterInput_t f_input;
      f_input.sigtype = sigtype;//geo.SignalType(wire.Channel());
      f_input.wireid  = wireID;//geo.ChannelToWire(wire.Channel()).at(0);
      f_input.view    = wire.View();
      f_input.channel = wire.Channel();
      f_input.roi_start = roi.begin_index();
      f_input.roi_end = roi.begin_index()+roi.size();
      f_input.roi_dataptr = &(roi.data());

      zmq_msg_t* msg_roi;
      zmq_msg_init_data(msg_roi,&f_input,sizeof(FitterInput_t),NULL,NULL);

      zmq_msg_send(fSender_ROI,msg_roi,0);
      ++total_rois;
      /*
      fFitter.RunFitter(roi.data());

      const float summedADCTotal = std::accumulate(roi.data().begin(),roi.data().end(),0.0);
      const raw::TDCtick_t startTick = roi.begin_index();
      const raw::TDCtick_t endTick = roi.begin_index()+roi.size();
      
      EmplaceHit(hitVector,summedADCTotal,startTick,endTick,sigtype,wireID,wire.Channel(),wire.View());
      */
      
    }//end loop over ROIs on wire

  }//end loop over wires

  //now we need to collect the hits. Exit only when we verify that our vector contains all hits we expect
  size_t rois_received=0;
  size_t expected_hits=0;
  while (rois_received<total_rois || hitVector.size() < expected_hits) {
    zmq_pollitem_t items [] = {
      { fReceiver_Hits, 0, ZMQ_POLLIN, 0 },
      { fReceiver_fin, 0, ZMQ_POLLIN, 0 }
    };
    zmq_poll (items, 2, -1);
    if (items [0].revents & ZMQ_POLLIN) {

      zmq_msg_t *msg_recv;
      zmq_msg_init(msg_recv);
      zmq_msg_recv(fReceiver_Hits,msg_recv);

      hitVector.emplace_back(*(recob::Hit*)(zmq_msg_data(msg_recv)));
      
    }
    if (items [1].revents & ZMQ_POLLIN){

      zmq_msg_t *msg_recv;
      zmq_msg_init(msg_recv);
      zmq_msg_recv(fReceiver_fin,msg_recv);

      expected_hits += *(size_t*)(zmq_msg_data(msg_recv));
      ++rois_received;
    }
  }

  //upon exit, we need to kill the workers
  zmq_send(fController, "KILL",4,0);

}

void hit::RFFHitFinderAlg::FitterWorker(){

  void *receiver = zmq_socket(context, ZMQ_PULL);
  zmq_connect (receiver, "inproc://workers");

  void *sender_hits = zmq_socket(context, ZMQ_PUSH);
  zmq_connect (sender_hits, "inproc://hits");

  void *sender_fin = zmq_socket(context, ZMQ_PUSH);
  zmq_connect (sender_fin, "inproc://finished");
  
  void *controller = zmq_socket (context, ZMQ_SUB);
  zmq_connect (controller, "inproc://control");
  zmq_setsockopt (controller, ZMQ_SUBSCRIBE, "", 0);

  while (1) {
    zmq_pollitem_t items [] = {
      { receiver, 0, ZMQ_POLLIN, 0 },
      { controller, 0, ZMQ_POLLIN, 0 }
    };
    zmq_poll (items, 2, -1);
    if (items [0].revents & ZMQ_POLLIN) {

      zmq_msg_t *msg_recv;
      zmq_msg_init(msg_recv);
      zmq_msg_recv(receiver,msg_recv);

      FitterInput_t f_input;
      memcpy(&f_input,zmq_msg_data(msg_recv),zmq_msg_size(msg_recv));

      zmq_msg_close(msg_recv);
      
      RFFHitFitter worker_fitter;
      worker_fitter.SetFitterParams(fMatchThresholdVec[f_input.view],
				    fMergeMultiplicityVec[f_input.view],
				    fAmpThresholdVec[f_input.view]);
      fFitter.RunFitter(*(f_input.roi_dataptr));

      std::vector<recob::Hit> hitVector;      
      EmplaceHit(hitVector,
		 std::accumulate(f_input.roi_dataptr->begin(),
				 f_input.roi_dataptr->end(),
				 0.0),
		 f_input.roi_start,
		 f_input.roi_end,
		 f_input.roi_sigtype,
		 f_input.roi_wireid,
		 f_input.channel,
		 f_input.view);

      for(size_t i_hit=0; i_hit<hitVector.size(); ++i_hit){
	zmq_msg_t *msg_send;
	zmq_msg_init_size(msg_send,sizeof(recob::Hit));
	memcpy(zmq_msg_data(msg_send),hitVector[i_hit],sizeof(recob::Hit));
	zmq_msg_send(sender_hits,msg_send,0);
      }
      zmq_msg_t *msg_send;
      zmq_msg_init_size(msg_send,sizeof(size_t));
      memcpy(zmq_msg_data(msg_send),hitVector.size(),sizeof(size_t));
      zmq_msg_send(sender_fin,msg_send,0);
      
    }
    //  Any waiting controller command acts as 'KILL'
    if (items [1].revents & ZMQ_POLLIN)
      break;                      //  Exit loop
  }
}

void hit::RFFHitFinderAlg::EmplaceHit(std::vector<recob::Hit>& hitVector,
				      float const& summedADCTotal,
				      raw::TDCtick_t const& startTick, raw::TDCtick_t const& endTick,
				      geo::SigType_t const& sigtype, geo::WireID const& wireID,
				      raw::ChannelID_t const& channel, geo::View_t const& view)
{

  float totalArea = 0.0;
  std::vector<float> areaVector(fFitter.NHits());
  std::vector<float> areaErrorVector(fFitter.NHits());
  std::vector<float> areaFracVector(fFitter.NHits());

  for(size_t ihit=0; ihit < fFitter.NHits(); ihit++){
    areaVector[ihit] = fFitter.AmplitudeVector()[ihit]*fFitter.SigmaVector()[ihit]*SQRT_TWO_PI;
    areaErrorVector[ihit] =
      SQRT_TWO_PI*std::sqrt(fFitter.AmplitudeVector()[ihit]*fFitter.SigmaErrorVector()[ihit]*fFitter.AmplitudeVector()[ihit]*fFitter.SigmaErrorVector()[ihit] +
			    fFitter.AmplitudeErrorVector()[ihit]*fFitter.SigmaVector()[ihit]*fFitter.AmplitudeErrorVector()[ihit]*fFitter.SigmaVector()[ihit]);
    totalArea += areaVector[ihit];
  }

  for(size_t ihit=0; ihit < fFitter.NHits(); ihit++){
    areaFracVector[ihit] = areaVector[ihit]/totalArea;
    
    hitVector.emplace_back(channel,
			   startTick,
			   endTick,
			   fFitter.MeanVector()[ihit]+(float)startTick,
			   fFitter.MeanErrorVector()[ihit],
			   fFitter.SigmaVector()[ihit],
			   fFitter.AmplitudeVector()[ihit],
			   fFitter.AmplitudeErrorVector()[ihit],
			   summedADCTotal*areaFracVector[ihit],
			   areaVector[ihit],
			   areaErrorVector[ihit],
			   fFitter.NHits(),
			   ihit,
			   -999.,
			   -999,
			   view,
			   sigtype,
			   wireID);
  }
			 
}
