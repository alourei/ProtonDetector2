//==========================================================================
// IonChamberSD.cc
// Based on DemonScintSD.cc and ExN04TrackerSD.cc
//
// Written/Modified by: Brian Roeder, TAMU 11/15/08
//                      email - broeder@comp.tamu.edu
//
// Purpose: Defines IonChamberSD member functions for simulation data readout
//          Sets and accesses Data in SiliconDetHit.hh methods
//
//==========================================================================
//
//
//


#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "PerfectDetHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"

#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4VProcess.hh"

#include "IonChamberSD.hh"


IonChamberSD::IonChamberSD(G4String name) : G4VSensitiveDetector(name)
{
  /* Constructor of Detector Method */
  // Note: Constructor runs before "run" begins
  // G4cout << "You are Using a IonChamber SD! " << G4endl;
  G4String HCname;
  collectionName.insert(HCname="IonChamberHitsCollection"); // Set name of Event Hits Collection
  MargotDataOutCT = MargotDataRecordTree::MargotPointer;
}

IonChamberSD::~IonChamberSD()
{/* Destructor */ }

void IonChamberSD::Initialize(G4HCofThisEvent* HCE)    // Member Function of IonChamberSD
{
  //- This method is run at the beginning of each event
  //- Sets up HitsCollection, etc.
  //- Note that "GetCollectionID()" is a slow op. - Following Demon Setup!

  IonChamberHitsCollection = new PerfectDetHitsCollection(SensitiveDetectorName,collectionName[0]);
  static int DetCollectionID = -1;

  if(DetCollectionID<0)
    { DetCollectionID = GetCollectionID(0); }  // Sets Hit CollectionID

  HCE->AddHitsCollection( DetCollectionID, IonChamberHitsCollection ); // Sets Events Hits Collection

  //The below is for debugging purposes
  /*
   G4cout << "DetCollection ID = " << DetCollectionID << G4endl;
   G4cout << "The Detector was initialized!" << G4endl;
   G4cout << "This is the address of the Silicon Hits Collection : " << SiliconHitsCollection << G4endl;
   G4cout << "This is the address of the  ->HitsCollection (HCE) : " << HCE << G4endl;
  */
}

G4bool IonChamberSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // Creates a Hit if EnergyDep > 0!

  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep == 0.)
    {return false;}

  PerfectDetHit* newHit = new PerfectDetHit();

  // newHit->SetKinEng( aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy() );
  newHit->SetKEDep ( edep );
  //newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() ); // Note - collects reaction point! 
  newHit->SetPos( aStep->GetPostStepPoint()->GetPosition() ); // Note - collects reaction point! 
  
  newHit->SetPDGCharge ( aStep->GetTrack()->GetDefinition()->GetPDGCharge() );
  newHit->SetPDGMass( aStep->GetTrack()->GetDefinition()->GetPDGMass() );

  G4double theTime = aStep->GetTrack()->GetGlobalTime();
  newHit->SetTOF( theTime );

  newHit->SetParticleName( aStep->GetTrack()->GetDefinition()->GetParticleName() );

  newHit->SetDetCopyNumber(aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0));
  newHit->SetDetName(aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(0)->GetName());

  const G4VProcess* theProcess = aStep->GetTrack()->GetCreatorProcess();
  if(theProcess != 0)
    {
      const G4String theProcessName = theProcess->GetProcessName();
      newHit->SetParticleProcess(theProcessName);
    }
  else
    { newHit->SetParticleProcess("NoReaction"); }

  IonChamberHitsCollection->insert( newHit );

  return true;
}

void IonChamberSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  // Note "EndOfEvent()" is name, NOT EndofEvent !

 G4SDManager* SDMan = G4SDManager::GetSDMpointer();
 G4String detName = "IonChamberSD";
 G4int collectionID = SDMan->GetCollectionID("IonChamberSD/IonChamberHitsCollection");

 PerfectDetHitsCollection* myCollection = (PerfectDetHitsCollection*)(HCE->GetHC(collectionID));

 G4String theParticleName = "none";
 G4String theDetName = "none";
 G4double KinEngDep = 0.;

 G4double MicroBulk_DeltaE[5]; 
 for(G4int i=0;i<5;i++)
   {MicroBulk_DeltaE[i] = 0.;}
 G4double MicroBulk_DeltaE_Total = 0.;
 G4double BulkFront_DeltaE = 0.;
 G4double BulkBack_DeltaE = 0.;

 G4ThreeVector thePos(0.0,0.0,0.0);
 G4ThreeVector theMicroBulkPos(0.0,0.0,0.0);
 G4ThreeVector theBulkFrontPos(0.0,0.0,0.0);
 G4ThreeVector theBulkBackPos(0.0,0.0,0.0);

 if(myCollection)
  {
    int n_hit = myCollection->entries();
 
    // Event readout loop!
    
    for(int i=0;i<n_hit;i++)
    {
       
      // (*myCollection)[i] is the pointer to the i-th hit of the event.
      // All data analysis output is read out here and then sent to the DataTree!

      PerfectDetHit* theCurrentHit = (*myCollection)[i];
     
      theParticleName = theCurrentHit->GetParticleName();
      KinEngDep = theCurrentHit->GetKEDep();
     
      //G4String theCurrentProcess = theCurrentHit->GetParticleProcess(); 
      G4double theCharge = theCurrentHit->GetPDGCharge();

      theParticleName = theCurrentHit->GetParticleName();
      theDetName = theCurrentHit->GetDetName(); 
      thePos = theCurrentHit->GetPos();

 
      if(theCharge > 2 && theDetName == "MicroBulk1")
	{
	  MicroBulk_DeltaE[0] += KinEngDep;
	  MicroBulk_DeltaE_Total += KinEngDep;
	  theMicroBulkPos = thePos;
	}

      if(theCharge > 2 && theDetName == "MicroBulk2")
	{
	  MicroBulk_DeltaE[1] += KinEngDep;
	  MicroBulk_DeltaE_Total += KinEngDep;
	  theMicroBulkPos = thePos;
	}

      if(theCharge > 2 && theDetName == "MicroBulk3")
	{
	  MicroBulk_DeltaE[2] += KinEngDep;
	  MicroBulk_DeltaE_Total += KinEngDep;
	  theMicroBulkPos = thePos;
	}

      if(theCharge > 2 && theDetName == "MicroBulk4")
	{
	  MicroBulk_DeltaE[3] += KinEngDep;
	  MicroBulk_DeltaE_Total += KinEngDep;
	  theMicroBulkPos = thePos;
	}

      if(theCharge > 2 && theDetName == "MicroBulk5")
	{
	  MicroBulk_DeltaE[4] += KinEngDep;
	  MicroBulk_DeltaE_Total += KinEngDep;
	  theMicroBulkPos = thePos;
	}


      
      if(theCharge > 2 && theDetName == "BulkFront")
	{
	  BulkFront_DeltaE += KinEngDep;
	  theBulkFrontPos = thePos;
	}

      if(theCharge > 2 && theDetName == "BulkBack")
	{
	  BulkBack_DeltaE += KinEngDep;
	  theBulkBackPos = thePos;
	}
    }
  
  }
 else
   {
     /* No hit collection! (although, always a HitCollection even if no Hits ! */
     G4cout << "Warning! No Hits Collection! " << G4endl;
   }

 MargotDataOutCT->sendDeltaEData(MicroBulk_DeltaE,theMicroBulkPos,MicroBulk_DeltaE_Total,BulkFront_DeltaE,theBulkFrontPos,BulkBack_DeltaE,theBulkBackPos);

 //G4cout <<<< ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>End of Event Analysis!" << G4endl;

} // End of EndOfEvent() Method
