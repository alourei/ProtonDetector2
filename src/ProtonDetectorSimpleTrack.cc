/*
 * ProtonDetectorSimpleTrack.cc
 *
 *  Created on: Dec 9, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorSimpleTrack.hh>

ClassImp(ProtonDetectorSimpleTrack)

ProtonDetectorSimpleTrack::ProtonDetectorSimpleTrack() {
	 //
	 // Constructor initializing to defaults
	 //

	  xPre = 0.;
	  yPre = 0.;
	  zPre = 0.;
	  xPost = 0.;
	  yPost = 0.;
	  zPost = 0.;
	  energyStride = 0.;
	  strideLength = 0.;
	  timePre = 0.;
	  timePost = 0.;
	  numberSteps = 0;
	  strideOrdinal = 0;
	  parentTrackID = -999;
	  trackID = -999;
	  eventID = -999;
	  runID = -999;
}

ProtonDetectorSimpleTrack::~ProtonDetectorSimpleTrack() {
	 //
	 // Destructor. Makes nothing
	 //
}

void ProtonDetectorSimpleTrack::Reset(){

	  //
	  // Clearing to defaults
	  //
	  xPre = 0.;
	  yPre = 0.;
	  zPre = 0.;
	  xPost = 0.;
	  yPost= 0.;
	  zPost = 0.;
	  energyStride = 0.;
	  particleCharge=0.;
	  particleMass=0.;
	  particleID=0;
	  strideLength = 0.;
	  timePre = 0.;
	  timePost = 0.;
	  numberSteps = 0;
	  strideOrdinal = 0;
	  parentTrackID = -999;
	  trackID = -999;
	  eventID = -999;
	  runID = -999;
}


ProtonDetectorSimpleTrack & ProtonDetectorSimpleTrack::operator=(const ProtonDetectorSimpleTrack &right){
  //
  // overloading the copy operator...
  //

  if (this != &right){
  xPre = right.xPre;
  yPre = right.yPre;
  zPre = right.zPre;
  xPost = right.xPost;
  yPost = right.yPost;
  zPost = right.zPost;
  energyStride = right.energyStride;
  particleCharge = right.particleCharge;
  particleMass = right.particleMass;
  particleID = right.particleID;
  strideLength = right.strideLength;
  timePre = right.timePre;
  timePost = right.timePost;
  numberSteps = right.numberSteps;
  strideOrdinal = right.strideOrdinal;
  parentTrackID = right.parentTrackID;
  trackID = right.trackID;
  eventID = right.eventID;
  runID = right.runID;
  }
  return *this;
}



