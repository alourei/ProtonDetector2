//==========================================================================
// IonChamberSD.hh
// Based on DemonScintSD.hh and ExN04TrackerSD.hh
//
// Written/Modified by: Brian Roeder, LPC Caen 02/14/07
//                      email - roeder@lpccaen.in2p3.fr
//
// Purpose: Defines IonChamberSD class for simulation data readout
//          Sets and accesses Data in SiliconDetHit.hh methods
//
//==========================================================================
//
// - See UM Hits Presentation and JLab Scoring 2 Talk for more info.
//
// 8/8/2007 - Redefined Readout in ProcessHits to act like UserSteppingAction
// and thus records all hits for neutrons and not just ones that deposit 
// energy.
//


#ifndef IonChamberSD_h
#define IonChamberSD_h 1

#include "G4VSensitiveDetector.hh"
#include "PerfectDetHit.hh"
#include "MargotDataRecordTree.hh"

class G4Step;

class G4HCofThisEvent;

class G4TouchableHistory;

class IonChamberSD : public G4VSensitiveDetector
{
public:

  IonChamberSD(G4String DetName);
  ~IonChamberSD();
  void Initialize(G4HCofThisEvent *HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);  // Will add TH Later
  void EndOfEvent(G4HCofThisEvent* );

private:

  PerfectDetHitsCollection* IonChamberHitsCollection;
  MargotDataRecordTree* MargotDataOutCT;
 
       
};

#endif
