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
// $Id: MargotDetectorConstruction.hh,v 1.6 2008/11/15 17:47:13 roeder Exp $
// GEANT4 tag $Name: geant4-09-01p02 $
//
// Modified by Brian Roeder, TAMU on 11/15/2008
// email - broeder@comp.tamu.edu
// 
// --------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other examples.
// --------------------------------------------------------------
//
//

#ifndef MargotDetectorConstruction_H
#define MargotDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class G4Box;
class G4Tubs;
class G4Sphere;
class G4Material;

class G4UnionSolid;
class G4SubtractionSolid;
class G4IntersectionSolid;


#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class MargotDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

  MargotDetectorConstruction(G4String DetType, G4double Thick, G4double DetRot);
  ~MargotDetectorConstruction();

  G4double CalculateDetTheta(G4double X, G4double Y, G4double Z);
  G4double CalculateDetPhi(G4double Theta, G4double X, G4double Y);
  G4double VectModulus(G4double X, G4double Y, G4double Z);

  G4RotationMatrix* Create_Rotation(G4double phi, G4double theta, G4double delta);
   
  G4VPhysicalVolume* Construct();
 
  private:

  // G4 - Elements and Materials
    G4Material* Al;
    G4Material* Fe;
    G4Material* Air;
    G4Material* Vacuum;

  G4Material* Ar;
  //G4Material* Isobutane;
  G4Material* Methane;

  G4Material* Kapton;
  //G4Material* ArIsobutane;
  G4Material* ArMethane;

//World
    G4Box*             solidWorld;    // pointer to the solid envelope 
    G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
    G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope

//Aluminum Degrader
    G4Box*             solidDegrader;    // pointer to the solid envelope 
    G4LogicalVolume*   logicDegrader;    // pointer to the logical envelope
    G4VPhysicalVolume* physiDegrader;    // pointer to the physical envelope

//Ion Chamber
    G4VSolid*            solidIonChamberP1;    // pointer to the solid envelope 
    G4VSolid*            solidIonChamberP2;    // pointer to the solid envelope 
    G4VSolid*      solidIonChamber;    // pointer to the solid envelope 
    G4LogicalVolume*   logicIonChamber;    // pointer to the logical envelope
    G4VPhysicalVolume* physiIonChamber;    // pointer to the physical envelope

//Ion Chamber Gas
    G4VSolid*          solidIonChamberGasP1;    // pointer to the solid envelope 
    G4VSolid*          solidIonChamberGasP2;    // pointer to the solid envelope 
    G4VSolid*          solidIonChamberGas;    // pointer to the solid envelope 
    G4LogicalVolume*   logicIonChamberGas;    // pointer to the logical envelope
    G4VPhysicalVolume* physiIonChamberGas;    // pointer to the physical envelope

//Ion Chamber Window
    G4Tubs*            solidWindow;    // pointer to the solid envelope 
    G4LogicalVolume*   logicWindow;    // pointer to the logical envelope
    G4VPhysicalVolume* physiWindowFront;    // pointer to the physical envelope
    G4VPhysicalVolume* physiWindowBack;    // pointer to the physical envelope

//Ion Chamber Gas - Inner Dets., MicroBulk
    G4VSolid*          solidMicroBulk;    // pointer to the solid envelope 
    G4LogicalVolume*   logicMicroBulk;    // pointer to the logical envelope
    G4VPhysicalVolume* physiMicroBulk1;    // pointer to the physical envelope
    G4VPhysicalVolume* physiMicroBulk2;    // pointer to the physical envelope
    G4VPhysicalVolume* physiMicroBulk3;    // pointer to the physical envelope

//Ion Chamber Gas - Outer Dets., Bulk Top/Bottom
    G4VSolid*          solidBulkTB;    // pointer to the solid envelope 
    G4LogicalVolume*   logicBulkTB;    // pointer to the logical envelope
    G4VPhysicalVolume* physiBulkTop1;    // pointer to the physical envelope
    G4VPhysicalVolume* physiBulkTop2;    // pointer to the physical envelope
    G4VPhysicalVolume* physiBulkBottom1; // pointer to the physical envelope
    G4VPhysicalVolume* physiBulkBottom2; // pointer to the physical envelope

//Ion Chamber Gas - Outer Dets., Bulk Front/Back
    G4VSolid*          solidBulkFB;    // pointer to the solid envelope 
    G4LogicalVolume*   logicBulkFB;    // pointer to the logical envelope
    G4VPhysicalVolume* physiBulkFront;    // pointer to the physical envelope
    G4VPhysicalVolume* physiBulkBack; // pointer to the physical envelope

//Ion Chamber Gas - Beam Counter, Bulk
    G4VSolid*          solidBeamCounter;    // pointer to the solid envelope 
    G4LogicalVolume*   logicBeamCounter;    // pointer to the logical envelope
    G4VPhysicalVolume* physiBeamCounterFront;    // pointer to the physical envelope
    G4VPhysicalVolume* physiBeamCounterBack; // pointer to the physical envelope

//Si Det
    G4Box*            solidSiDet;    // pointer to the solid envelope 
    G4LogicalVolume*   logicSiDet;    // pointer to the logical envelope
    G4VPhysicalVolume* physiSiDet;    // pointer to the physical envelope

private:

  G4double WorldSize;
  G4String DetectorType;
  G4double DegThickness;
  G4double DegRotAngle;

public:

  // These functions allow other parts of the program to access 
  // detector construction information

  G4double GetWorldLength()
  {return WorldSize;}

};

#endif

