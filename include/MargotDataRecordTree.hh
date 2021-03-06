//-----------------------------------------------------------------
// MargotDataRecordTree.hh
//
// - written by: Brian Roeder, LPC Caen, 18 Jan 07 (in Serra)
// - email - roeder@lpccaen.in2p3.fr
//
// - modified: 14/02/07 for use with Margot Tracking
//
// - Usage: A C++ class for GEANT4 for creating a ROOT tree for 
// -        the Margot neutron detector simulation event by event.
//
///////////////////////////////////////////////////////////////////
//
//

#ifndef DATARECORD_H
#define DATARECORD_H

#include <iostream>
using namespace std;

// C++ formatting and file io headers

#include <iomanip>
#include <fstream>
#include <cmath>

// Root Analysis header files for C++

#include "Riostream.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "globals.hh"
#include "G4UnitsTable.hh"

#include "G4ThreeVector.hh"

class MargotDataRecordTree
{

private:

  // Initialized in class constructor in MargotDataRecordTree.cc 
  
  TFile* DataFile;
  TTree* IonChDistTree;

  TH1F* h1_ZStopPos;

  G4double eng_int;
  G4double theta_int;
  G4double phi_int;

  G4double eng_MicroBulk[5];
  G4double eng_MicroBulk_Tot;
  G4double X_MicroBulk;
  G4double Y_MicroBulk;
  G4double Z_MicroBulk;
  G4double PosMag_MicroBulk;

  G4double eng_BulkFront;
  G4double X_BulkFront;
  G4double Y_BulkFront;
  G4double Z_BulkFront;
  G4double PosMag_BulkFront;

  G4double eng_BulkBack;
  G4double X_BulkBack;
  G4double Y_BulkBack;
  G4double Z_BulkBack;
  G4double PosMag_BulkBack;

  G4double engDepEvent;
  G4double PosMagEvent;
 
  // Particle Counters

  int event_counter;
  int number_detected;
  int number_MicroBulk;
  int number_BulkFront;
  int number_BulkBack;
  int number_PunchThrough;
    
 public:
    MargotDataRecordTree();
   ~MargotDataRecordTree();
  
 // Creates Pointer to Data Analysis in main() program (Margot.cc).    
  static MargotDataRecordTree* MargotPointer;
 
  // Definition of member functions - access and distribute data
  // See "MargotDataRecordTree.cc for the function listings.

  void senddataPG(double value1, G4double theta, G4double phi);
  void ShowDataFromEvent();
  void GetParticleTotals();
  void sendDeltaEData(G4double DE1[], G4ThreeVector pos1,G4double DE1_Tot,G4double DE2, G4ThreeVector pos2, G4double DE3, G4ThreeVector pos3);
  
}; 
#endif
