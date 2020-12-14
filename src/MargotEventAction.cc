//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: MargotEventAction.cc,v 1.2 2007/02/14 17:45:14 roeder Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Modified by Brian Roeder, LPC Caen, monitor of events for Margot
// email - roeder@lpccaen.in2p3.fr
//
// 2/15/07 - Now reads events generated in MargotSD from MargotDetHit
 
#include "MargotEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"


#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

// Added following classes for Hit Map Printing

#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"

// These classes are for scorers (being phased out)
/*
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
*/

// 19 Feb 07 - BTR - Root Files now made in MargotDataRecordTree 
//                   Data Analysis/Readout done in MargotSD.cc - EndOfEvent()
//

MargotEventAction::MargotEventAction() :
  numberOfEvent(-1)
{
  /* Constructor - Initializes Counters and Analysis Pointer */
  MargotGetDataEV = MargotDataRecordTree::MargotPointer;
}

//
MargotEventAction::~MargotEventAction()
{
/* Destructor */
}

//
void MargotEventAction::BeginOfEventAction(const G4Event*)
{}

//
void MargotEventAction::EndOfEventAction(const G4Event* evt)
{
  
  // Periodically Print out Event Results!
  // Analysis of Event done in MargotSD.cc EndOfEvent() Method.

 numberOfEvent++;
 // G4cout << "************************ End of Event # " << numberOfEvent << G4endl;

  //  Periodic Printout of event  //

 //G4RunManager* runManager = G4RunManager::GetRunManager();
 // G4SDManager* SDMan = G4SDManager::GetSDMpointer();

 // G4String detName = "MargotSD";
 // G4int collectionID = SDMan->GetCollectionID("MargotSD/MargotHitsCollection");
 //const G4Event* currentEvent = runManager->GetCurrentEvent();

 // G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
 // MargotDetHitsCollection* myCollection = (MargotDetHitsCollection*)(HCE->GetHC(collectionID));

      G4int event_show = 2000;   
      if((numberOfEvent) % event_show == 0) 
      {
      G4cout << "====================================================" << G4endl;
      G4cout << ">>> Event " << evt->GetEventID() << G4endl;

      // G4int n_hits = myCollection->entries();
      // G4cout << n_hits << " Hits were registered in this event! " << G4endl; 
      MargotGetDataEV->ShowDataFromEvent();
      }
  
}  // End of "End of Event Action"
