/*
 * ProtonDedtectorROOTAnalysis.hh
 *
 *  Created on: Dec 10, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORROOTANALYSIS_HH_
#define PROTONDETECTORROOTANALYSIS_HH_

//G4includes
#include "G4ClassificationOfNewTrack.hh"
#include "G4TrackStatus.hh"
#include "G4Types.hh"
#include "G4PrimaryParticle.hh"

//ROOT includes
#include "TROOT.h"


//Geant4 classes
class G4VPhysicalVolume;
class G4Event;
class G4Run;
class G4Track;
class G4Step;

//ROOT classes
class TTree;
class TBranch;
class TFile;
class TClonesArray;

//ProtonDetector
#include <ProtonDetectorGeantHit.hh>

//ProtonDetector classes
class ProtonDetectorAnalysisMessenger;


class ProtonDetectorConstruction;
class ProtonDetectorPrimaryGeneratorAction;

//class ProtonDetectorData;
class ProtonDetectorTrack;
class ProtonDetectorSimpleTrack;
class ProtonDetectorBeamInfo;
class ProtonDetectorData;

class ProtonDetectorROOTAnalysis;
extern ProtonDetectorROOTAnalysis *gProtonDetectorROOTAnalysis; // global


class ProtonDetectorROOTAnalysis {

private:
	  time_t LastDoItTime; // used in OnceAWhileDoIt method

	  TFile* simFile;     // The ROOT File
	  char* newDirName;

	  TTree* eventTree; //Tree
	  TTree* tracksTree; //Tree

	  G4PrimaryParticle* primary;               //Storing the primary for accesing during UserStep
	  //ProtonDetectorAnalysisMessenger* analMessenger; // pointer to messenger

	  ProtonDetectorBeamInfo *theBeamInfo;
	  ProtonDetectorData *theData;

	  TClonesArray*  simpleTrackCA;
	  ProtonDetectorSimpleTrack** simpleTrack; //the two simple data track

	  TBranch* simpleTrackBranch; //Local branch
	  TBranch* trackBranch; //Local branch
	  TBranch* beamInfoBrach;

	  ProtonDetectorTrack *theTracks;

	  //ProtonDetectorPrimaryInfo** thePrimaryInfo; //Primary particles data
	  //TClonesArray* primaryInfoCA;

	  G4int theRunID; //To keep some numbers on the Tree
	  G4int theEventID; //To keep some numbers on the Tree

	  //Flags for control of gas tracks...
	  G4String  storeTracksFlag;
	  G4String  storeTrackHistosFlag;
	  G4String  storeEventsFlag;
	  G4String  storeSimpleTracksFlag;
	  G4String  beamInteractionFlag;  //flag to turn "on"/"off" the beam interaction analysis

	  G4double minStrideLength;

public:
	ProtonDetectorROOTAnalysis();
	~ProtonDetectorROOTAnalysis();

	void InitAnalysis();

	  TTree* GetEventTree(){return eventTree;}
	  void SetEventTree(TTree* tree) {eventTree = tree;}
	  TTree* GetTracksTree(){return tracksTree;}
	  void SetTracksTree(TTree* tree) {tracksTree = tree;}

	  G4int GetTheEventID(){return theEventID;}
	  void SetTheEventID(G4int id){theEventID = id;}
	  G4int GetTheRunID(){return theRunID;}
	  void SetTheRunID(G4int id){theRunID = id;}

	  TBranch* GetTrackBranch(){return trackBranch;}
	  void SetTrackBranch(TBranch* aBranch) {trackBranch= aBranch;}
	  TBranch* GetSimpleTrackBranch(){return simpleTrackBranch;}
	  void SetSimpleTrackBranch(TBranch* aBranch) {simpleTrackBranch= aBranch;}

	  TClonesArray* getSimpleTrackCA(void){return simpleTrackCA;}
	  void SetSimpleTrackCA(TClonesArray* CA) {simpleTrackCA = CA;}



	  //Messenger actions
	  void SetStoreTracksFlag(G4String val) {storeTracksFlag = val;};
	  void SetStoreTrackHistosFlag(G4String val) {storeTrackHistosFlag = val;};
	  void SetStoreEventsFlag(G4String val) {storeEventsFlag = val;};
	  void SetStoreSimpleTracksFlag(G4String val) {storeSimpleTracksFlag=val;};

	  void InitAnalysisForExistingDetectors();

	  //DPL 29NOV2012
	  void SetMinStrideLength(Double_t value);

	 // G4VUserDetectorConstruction
	  void Construct(const G4VPhysicalVolume*);

	  // G4VUserPhysicsList
	  void ConstructParticle();
	  void ConstructProcess();
	  void SetCuts();

	  // G4VUserPrimaryGeneratorAction
	  //TODO->Solve this assymetry!
	  //void GeneratePrimaries(const G4Event*);
	  //void GeneratePrimaries(const G4Event*,G4double,G4double,G4double,G4double);
	  void GenerateBeam(const G4Event*);

	  // G4UserRunAction
	  void BeginOfRunAction(const G4Run*);
	  void EndOfRunAction(const G4Run*);

	  // G4UserEventAction
	  void BeginOfEventAction(const G4Event*);
	  void EndOfEventAction(const G4Event*);

	  // G4UserStackingAction
	  void ClassifyNewTrack(const G4Track*, G4ClassificationOfNewTrack*);
	  void NewStage();
	  void PrepareNewEvent();

	  // G4UserTrackingAction
	  void PreUserTrackingAction(const G4Track*);
	  void PostUserTrackingAction(const G4Track*, G4TrackStatus*);

	  // G4UserSteppingAction
	  void UserSteppingAction(const G4Step*); // original

	  // once a while do "something"
	  void OnceAWhileDoIt(const G4bool DoItNow = false);

};

#endif /* PROTONDETECTORROOTANALYSIS_HH_ */
