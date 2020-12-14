/*
 * ProtonDetectorTrack.hh
 *
 *  Created on: Dec 10, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORTRACK_HH_
#define PROTONDETECTORTRACK_HH_

#include "TROOT.h"  //for including Rtypes.h

class ProtonDetectorTrack {
private:
	 Double_t xCoord;
	  Double_t yCoord;
	  Double_t zCoord;
	  Double_t xPreCoord;
	  Double_t yPreCoord;
	  Double_t zPreCoord;
	  Double_t energyStep;
	  Int_t parentTrackID;
	  Int_t trackID;
	  Int_t eventID;
	  Int_t runID;

public:
	ProtonDetectorTrack();
	~ProtonDetectorTrack();

	 Double_t GetXCoord(){return xCoord;}
	  Double_t GetYCoord(){return yCoord;}
	  Double_t GetZCoord(){return zCoord;}
	  Double_t GetXPreCoord(){return xPreCoord;}
	  Double_t GetYPreCoord(){return yPreCoord;}
	  Double_t GetZPreCoord(){return zPreCoord;}
	  Double_t GetEnergyStep(){return energyStep;}
	  Int_t GetTrackID(){return trackID;}
	  Int_t GetParentTrackID(){return parentTrackID;}
	  Int_t GetEventID(){return eventID;}
	  Int_t GetRunID(){return runID;}

	  void SetXCoord(Double_t xc){xCoord = xc;}
	  void SetYCoord(Double_t yc){yCoord = yc;}
	  void SetZCoord(Double_t zc){zCoord = zc;}
	  void SetXPreCoord(Double_t xc){xPreCoord = xc;}
	  void SetYPreCoord(Double_t yc){yPreCoord = yc;}
	  void SetZPreCoord(Double_t zc){zPreCoord = zc;}
	  void SetEnergyStep(Double_t energy){energyStep = energy;}
	  void SetTrackID(Int_t track){trackID = track;}
	  void SetParentTrackID(Int_t pt){parentTrackID = pt;}
	  void SetEventID(Int_t ev){eventID = ev;}
	  void SetRunID(Int_t ev){runID = ev;}

	  ClassDef(ProtonDetectorTrack,1) //ROOT CINT

};

#endif /* PROTONDETECTORTRACK_HH_ */
