/*
 * ProtonDetectorTrack.cc
 *
 *  Created on: Dec 10, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorTrack.hh>

ClassImp(ProtonDetectorTrack)

ProtonDetectorTrack::ProtonDetectorTrack() {

	xCoord = 0.;
	  yCoord = 0.;
	  zCoord = 0.;
	  xPreCoord = 0.;
	  yPreCoord = 0.;
	  zPreCoord = 0.;
	  energyStep = 0.;
	  parentTrackID = 0;
	  trackID = 0;
	  eventID = 0;
	  runID = 0;


}

ProtonDetectorTrack::~ProtonDetectorTrack() {
	// TODO Auto-generated destructor stub
}

