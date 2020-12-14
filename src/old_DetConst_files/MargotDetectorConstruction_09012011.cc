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
// $Id: MargotDetectorConstruction.cc,v 1.9 2010/09/29 17:47:19 roeder Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// Modified by Brian Roeder, TAMU on 09/29/2010
// email - broeder@comp.tamu.edu
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
//


#include <iomanip>
#include <fstream>
#include <cmath>

#include "MargotDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"
#include "G4ios.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"

// The methods below implement a user defined sensitive detector

#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "IonChamberSD.hh" 

MargotDetectorConstruction::MargotDetectorConstruction(G4String DetType, G4double Thick, G4double DegRot)
 :solidWorld(0),  logicWorld(0),  physiWorld(0),
  solidDegrader(0), logicDegrader(0), physiDegrader(0),
  solidIonChamberP1(0), solidIonChamberP2(0),
  solidIonChamber(0), logicIonChamber(0), physiIonChamber(0),
  solidIonChamberGasP1(0), solidIonChamberGasP2(0),
  solidIonChamberGas(0), logicIonChamberGas(0), physiIonChamberGas(0),
  solidWindow(0), logicWindow(0), physiWindowFront(0), physiWindowBack(0),
  solidMicroBulk(0), logicMicroBulk(0), physiMicroBulk1(0), physiMicroBulk2(0), physiMicroBulk3(0),
  solidBulkTB(0), logicBulkTB(0), physiBulkTop1(0), physiBulkTop2(0), physiBulkBottom1(0), physiBulkBottom2(0),
  solidBulkFB(0), logicBulkFB(0), physiBulkFront(0), physiBulkBack(0),
  solidBeamCounter(0), logicBeamCounter(0), physiBeamCounterFront(0), physiBeamCounterBack(0),
  solidSiDet(0), logicSiDet(0), physiSiDet(0),
  DetectorType(DetType), DegThickness(Thick), DegRotAngle(DegRot)
{;}

MargotDetectorConstruction::~MargotDetectorConstruction()
{;}

G4double MargotDetectorConstruction::VectModulus(G4double X, G4double Y, G4double Z)
{
  //Gives Modulus (normal) of a 3-vector 
  G4double V_Modulus = sqrt(pow(X,2)+pow(Y,2)+pow(Z,2));
  return V_Modulus;
} 

G4double MargotDetectorConstruction::CalculateDetTheta(G4double X, G4double Y, G4double Z)
{
  // Calculates Rotation angle of demon module based on position of det.
  // w.r.t. to target.
  G4double R = VectModulus(X,Y,Z);
  G4double theta = acos(Z/R);

  return theta;  // Returns theta in Radians
}

G4double MargotDetectorConstruction::CalculateDetPhi(G4double Theta, G4double X, G4double Y)
{
  const G4double Pi = CLHEP::pi;
  G4double Phi = 0.;
  if(Theta != 0.)
    {
      G4double X_sq = pow(X,2);
      G4double Y_sq = pow(Y,2);
      G4double Proj_Mag = sqrt(X_sq+Y_sq);

      //Phi = (asin(Y/Proj_Mag));
      
      if(X >= 0. && Y >= 0.)  // quad 1
	{ Phi = (asin(Y/Proj_Mag)); }
      else if(X < 0. && Y > 0.)  // quad 2
	{ Phi = (Pi - asin(Y/Proj_Mag)); } 
      else if(X <= 0. && Y <= 0.)  // quad 3
	{ Phi = (Pi - asin(Y/Proj_Mag)); } 
      else if(X > 0. && Y < 0.)  // quad 4
	{ Phi = (2*Pi + asin(Y/Proj_Mag)); }
     
    }
 else
   {Phi = 0.;}
  return Phi;    // Returns Phi in Radians
}

G4RotationMatrix* MargotDetectorConstruction::Create_Rotation(G4double phi, G4double theta, G4double delta)
{
  CLHEP::HepRotation r1,r2,r3;
  r1.rotateZ(-phi);
  r2.rotateY(-theta);
  r3.rotateZ(-delta);

  CLHEP::HepRotation rot = r3*r2*r1;

  G4ThreeVector col1 = rot.colX();
  G4ThreeVector col2 = rot.colY();
  G4ThreeVector col3 = rot.colZ(); 
  
  G4RotationMatrix* rot1 = new G4RotationMatrix(col1,col2,col3);
  return rot1;
}


G4VPhysicalVolume* MargotDetectorConstruction::Construct()
{
  // Elements and  materials used ---------------------------------------

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density, temperature, pressure;
  G4int nel; // number of elements in compound

  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
  //G4Element* H  = new G4Element("Hydrogen","H" , z= 1., a= 1.01*g/mole);
  //G4Element* C  = new G4Element("Carbon"  ,"C" , z= 6., a= 12.011*g/mole); 
  
  G4NistManager* NistMan = G4NistManager::Instance();
 
  Al = NistMan->FindOrBuildMaterial("G4_Al");
  Kapton = NistMan->FindOrBuildMaterial("G4_KAPTON");
  Fe =  NistMan->FindOrBuildMaterial("G4_Fe");

  Ar = NistMan->FindOrBuildMaterial("G4_Ar");
  Methane = NistMan->FindOrBuildMaterial("G4_METHANE");

  //Isobutane = NistMan->FindOrBuildMaterial("G4_BUTANE");
  //Ar = NistMan->ConstructNewGasMaterial("Argon","G4_Ar",temperature=273.15*kelvin,pressure=13332.*pascal,0);

  // LISE gas densities
  // Density Ar at 100 torr, 293.15K is 0.00021852 g/cm3
  // Density Butane at 100 torr, 293.15K is 0.00031795 g/cm3
  // 100 torr - Assume density of 95% Ar, 5% Butane mix = 0.0002234915 g/cm3
 
  // --- 1.0 atm ---
  // --- Construct others with ratio.
  // Density Ar at 760 torr, 293.15K is 0.00166078 g/cm3
  // Density Butane at 760 torr, 293.15K is 0.0024164 g/cm3
  // 760 torr - Assume density of 95% Ar, 5% Butane mix = 0.001698561 g/cm3
  // 760 torr = 1.0*atm
  //density = 0.001698561*g/cm3;
  //pressure = 1.0*atmosphere;

  // We used this pressure for Astrobox 1! 
  // Density Ar at 800 torr, 293.15K is 0.00166078 g/cm3
  // Density Butane at 800 torr, 293.15K is 0.0024164 g/cm3
  // 800 torr - Assume density of 95% Ar, 5% Butane mix = 0.002038 g/cm3
  // 800 torr = 1.0526*atm
  //density = (1.0526*0.001698561)*g/cm3;
  //pressure = 1.05026*atmosphere;

  //ArIsobutane = new G4Material("ArIsobutane", density, nel=2, kStateGas, 
  //			       temperature= 273.15*kelvin, pressure);
  //ArIsobutane->AddMaterial(Ar, 0.95);
  //ArIsobutane->AddMaterial(Isobutane, 0.05);

  //---------------------------------------------------------------------
  // Used in the experiments Oct. 2011
  // P5 at 800 torr
  // Density Ar at 760 torr, 293.15K is 0.00166078 g/cm3 
  // Density Methane at 760 torr, 293.15K is 0.000667 g/cm3
  // --- P5 density at 760 torr  -> 0.001611091 g/cm3

  G4double PressureFactor = 1.05263;  // 800 torr
  density = PressureFactor*0.001611091*g/cm3;
  pressure = PressureFactor*atmosphere; 
  ArMethane = new G4Material("ArMethane", density, nel=2, kStateGas, 
			       temperature= 273.15*kelvin, pressure);
  ArMethane->AddMaterial(Ar, 0.95);
  ArMethane->AddMaterial(Methane, 0.05);  //P5



  //Silicon (for Silicon Detector)
  density = 2.33*g/cm3;
  G4Material* Si = new G4Material("Si", z=14., a=28.0855*g/mole, density);

  // Air
    density = 1.205*mg/cm3;
    Air = new G4Material("Air", density, nel=2);
    Air->AddElement(N, .7);
    Air->AddElement(O, .3);

 //  Vacuum

    density = 0.00000001*mg/cm3;
    Vacuum = new G4Material("Vacuum", density, nel=2, kStateGas, 
					temperature= 273.15*kelvin, pressure=0.0001*pascal);
    Vacuum->AddElement(N, .7);
    Vacuum->AddElement(O, .3);

  // Print all the materials defined.
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //World Volume - A cube of 1 m^3 
  //Creates cubic world, note pointer "solidWorld" predefined in header

  WorldSize = 0.5*m;

  solidWorld = new G4Box("World",WorldSize/2,WorldSize/2,WorldSize/2); 

  //World logic volume - filled with vacuum 
  logicWorld = new G4LogicalVolume(solidWorld,         //solid of logical World
				   Vacuum,             //material filling logic volume
				   "World",0,0,0);     //Name, other properties 0.

  //World Physical volume (placement)
  physiWorld = new G4PVPlacement(0,                // no rotation
				 G4ThreeVector(),  // position at 0,0,0
				 logicWorld,       // its logical volume
				 "World",          // its name
				 0,                // its mother volume;(0) for world
				 false,            // no boolean op.
				 0);               // not a copy

  //World Visualization (see Demon)
  G4VisAttributes* WorldVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  WorldVisAtt->SetForceWireframe(true);
  //logicWorld->SetVisAttributes(WorldVisAtt);
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  // Size Variables for the Beta Detector Setup ----------------------------

  // Degrader -------------------------
 
   G4double RotatedDegThickness = DegThickness/(cos(DegRotAngle));
   G4cout << "The Degrader z-axis thickness = " 
	  << G4BestUnit(RotatedDegThickness,"Length") << G4endl;

 // Aluminum energy degrader -
  
    G4double DegraderSize = 5.*cm;

    solidDegrader = new G4Box("AlDegrader",DegraderSize/2.,DegraderSize/2.,RotatedDegThickness/2.);
 
    logicDegrader = new G4LogicalVolume(solidDegrader,    // Solid
				       Al,                // Material 
				       "AlDegrader",0,0,0); // Name, other properties 0.

    G4ThreeVector positionDegrader = G4ThreeVector(0.0,0.0,-5.0*cm);

    G4RotationMatrix* DegRotation = new G4RotationMatrix();
    DegRotation->rotateY(DegRotAngle);

    physiDegrader = new G4PVPlacement(DegRotation,       // no rotation
				      positionDegrader, // position defined above
				      logicDegrader,    // its logical volume
                                      "AlDegrader",       // its name
				      logicWorld,       // its mother  volume
				      false,            // no boolean operations
				      0);               // copy number

 // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* DegraderVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));  //Red 
  DegraderVisAtt->SetForceSolid(true);
  logicDegrader->SetVisAttributes(DegraderVisAtt);


 // Ion Chamber Mother Volume - Vacuum
 // Updated 8/16/2011 to be a "Union solid", changed to be
 // retangular for AstroBox2!
  
  // Solid for the length of the chamber - a box for AstroBox2
  // Right now only considers "field cage"
  // box is 160mm long (along beam axis)
  //        90mm  wide (along detection, horizontal axis)
  //        90mm  high (along vertical axis)

  // Assume same circular inlets as from AstroBox1.

  G4double IonChamberP1x = 95.0*mm;
  G4double IonChamberP1y = 95.0*mm;
  G4double IonChamberP1z = 165.00*mm;   

  solidIonChamberP1 = new G4Box("IonChamberP1",IonChamberP1x/2.,IonChamberP1y/2.,IonChamberP1z/2.);
 
    G4double StartAngle = 0.*deg;
    G4double SpanAngle = 360.*deg;
    G4double ICRminP2 = 0.*mm;
    G4double ICRmaxP2 = (25.0/2.)*mm;
    G4double halfHeightICP2 = (190.0/2.)*mm; // center is at 95.0mm w.r.t. to front window

    solidIonChamberP2 = new G4Tubs("IonChamberP2",ICRminP2,ICRmaxP2,halfHeightICP2,StartAngle,SpanAngle);

    G4ThreeVector ICUnionTrans = G4ThreeVector(0.,0.,0.);

    G4RotationMatrix* ICUnionRot = new G4RotationMatrix();
    ICUnionRot->rotateY(0.*deg);

    solidIonChamber = new G4UnionSolid("IonChamber",solidIonChamberP1,solidIonChamberP2,ICUnionRot,ICUnionTrans);
 
    logicIonChamber = new G4LogicalVolume(solidIonChamber,    // Solid
				       Fe,                    // Material 
				       "IonChamber",0,0,0);   // Name, other properties 0.

    G4ThreeVector positionIonChamber = G4ThreeVector(0.0,0.0,halfHeightICP2);

    G4RotationMatrix* ICRot = new G4RotationMatrix();
    ICRot->rotateY(0.*deg);

    physiIonChamber = new G4PVPlacement(ICRot,                // rotation
				      positionIonChamber, // position defined above
				      logicIonChamber,    // its logical volume
                                      "IonChamber",       // its name
				      logicWorld,         // its mother  volume
				      false,              // no boolean operations
				      0);                 // copy number

 // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* IonChamberVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.2,0.0));  //Green
  IonChamberVisAtt->SetForceWireframe(true);
  //logicIonChamber->SetVisAttributes(G4VisAttributes::Invisible);
  logicIonChamber->SetVisAttributes(IonChamberVisAtt);

  // The chamber interior with gas
  
    G4double ICGasP1x = 90.0*mm;
    G4double ICGasP1y = 90.0*mm;
    G4double ICGasP1z = 160.0*mm;

    solidIonChamberGasP1 = new G4Box("IonChamberGasP1",ICGasP1x/2.,ICGasP1y/2.,ICGasP1z/2.);

    G4double ICGasRminP2 = 0.*mm;
    G4double ICGasRmaxP2 = (20.00/2.)*mm;
    G4double halfHeightICGasP2 = ((190.0/2.))*mm;

    solidIonChamberGasP2 = new G4Tubs("IonChamberGasP2",ICGasRminP2,ICGasRmaxP2,halfHeightICGasP2,StartAngle,SpanAngle);

    G4ThreeVector ICGasUnionTrans = G4ThreeVector(0.,0.,0.);

    G4RotationMatrix* ICGasUnionRot = new G4RotationMatrix();
    ICGasUnionRot->rotateY(0.*deg);

    solidIonChamberGas = new G4UnionSolid("IonChamberGas",solidIonChamberGasP1,solidIonChamberGasP2,ICGasUnionRot,ICGasUnionTrans);
 
    logicIonChamberGas = new G4LogicalVolume(solidIonChamberGas,    // Solid
				             ArMethane,                // Material 
				            "IonChamberGas",0,0,0);   // Name, other properties 0.

    G4ThreeVector positionIonChamberGas = G4ThreeVector(0.0,0.0,0.0);

    physiIonChamberGas = new G4PVPlacement(0,                // rotation
				      positionIonChamberGas, // position defined above
				      logicIonChamberGas,    // its logical volume
                                      "IonChamberGas",       // its name
				      logicIonChamber,         // its mother  volume
				      false,              // no boolean operations
				      0);                 // copy number

 // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* IonChamberGasVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.4,0.0));  //Green 
  IonChamberGasVisAtt->SetForceSolid(true);
  logicIonChamberGas->SetVisAttributes(IonChamberVisAtt);
  //logicIonChamberGas->SetVisAttributes(G4VisAttributes::Invisible);
  
// Ion Chamber Window
  
    G4double ICRminWin = 0.*mm;
    G4double ICRmaxWin = (20.0/2.)*mm;  //diameter/2.
    G4double halfHeightWin = (0.0508/2.)*mm; // 2 mils thick
    StartAngle = 0.*deg;
    SpanAngle = 360.*deg;

    solidWindow = new G4Tubs("KaptonWindow",ICRminWin,ICRmaxWin,halfHeightWin,StartAngle,SpanAngle);
 
    logicWindow = new G4LogicalVolume(solidWindow,        // Solid
				       Kapton,            // Material 
				       "KaptonWindow",0,0,0);   // Name, other properties 0.

    G4ThreeVector positionWindowFront = G4ThreeVector(0.0,0.0,-(halfHeightICGasP2-halfHeightWin)); // w.r.t. Ion Chamber

    G4RotationMatrix* WindowRot = new G4RotationMatrix();
    WindowRot->rotateY(0.*deg);

    physiWindowFront     = new G4PVPlacement(WindowRot,                // no rotation
					 positionWindowFront, // position defined above
					 logicWindow,    // its logical volume
					 "KaptonWindowFront",       // its name
					 logicIonChamberGas,      // its mother  volume
					 false,              // no boolean operations
					 0);                 // copy number

    G4ThreeVector positionWindowBack = G4ThreeVector(0.0,0.0,(halfHeightICGasP2-halfHeightWin)); // w.r.t. Ion Chamber

     physiWindowBack     = new G4PVPlacement(WindowRot,                // no rotation
					 positionWindowBack, // position defined above
					 logicWindow,    // its logical volume
					 "KaptonWindowBack",       // its name
					 logicIonChamberGas,      // its mother  volume
					 false,              // no boolean operations
					 1);                 // copy number


 // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* WindowVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));  //Blue 
  WindowVisAtt->SetForceSolid(true);
  logicWindow->SetVisAttributes(WindowVisAtt);

 // Ion Chamber Gas - Micro-Bulk Volumes --------------------------------------
 // Now simulated as "Box" Volumes

    G4double MicroBulkX = 90.0*mm;
    G4double MicroBulkY = 30.0*mm;
    G4double MicroBulkZ = 30.0*mm;
   
    solidMicroBulk = new G4Box("MicroBulk",MicroBulkX/2.,MicroBulkY/2.,MicroBulkZ/2.);
 
    logicMicroBulk = new G4LogicalVolume(solidMicroBulk,        // Solid
				       ArMethane,               // Material 
				       "MicroBulk",0,0,0);   // Name, other properties 0.

    G4ThreeVector positionMicroBulk1 = G4ThreeVector(0.0,0.0,-(MicroBulkZ+1.0*mm)); // w.r.t. Ion Chamber Gas

    physiMicroBulk1 = new G4PVPlacement(0,                // no rotation
				      positionMicroBulk1, // position defined above
				      logicMicroBulk,    // its logical volume
                                      "MicroBulk1",       // its name
				      logicIonChamberGas,  // its mother  volume
				      false,              // no boolean operations
				      1);                 // copy number

    G4ThreeVector positionMicroBulk2 = G4ThreeVector(0.0,0.0,0.0); // w.r.t. Ion Chamber Gas

    physiMicroBulk2 = new G4PVPlacement(0,                // no rotation
				      positionMicroBulk2, // position defined above
				      logicMicroBulk,    // its logical volume
                                      "MicroBulk2",       // its name
				      logicIonChamberGas,  // its mother  volume
				      false,              // no boolean operations
				      2);                 // copy number

    G4ThreeVector positionMicroBulk3 = G4ThreeVector(0.0,0.0,(MicroBulkZ+1.0*mm)); // w.r.t. Ion Chamber Gas

    physiMicroBulk3 = new G4PVPlacement(0,                // no rotation
				      positionMicroBulk3, // position defined above
				      logicMicroBulk,    // its logical volume
                                      "MicroBulk3",       // its name
				      logicIonChamberGas,  // its mother  volume
				      false,              // no boolean operations
				      3);                 // copy number




 // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* MicroBulkVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));  //Red 
  MicroBulkVisAtt->SetForceSolid(true);
  logicMicroBulk->SetVisAttributes(MicroBulkVisAtt);
 
  // Ion Chamber Gas - Outer Dets. Made of Bulk Detectors --------------------------------------
  
  G4double BulkTBx = 90.0*mm;
  G4double BulkTBy = 20.0*mm;  // Total active Y length = 70mm, so 70mm-30mm(MicroBulk) = 40.mm/2 (top/bottom)
  G4double BulkTBz = 45.0*mm;
 
    solidBulkTB = new G4Box("BulkTB",BulkTBx/2.,BulkTBy/2.,BulkTBz/2.);
 
    logicBulkTB = new G4LogicalVolume(solidBulkTB,        // Solid
				       ArMethane,               // Material 
				       "BulkTB",0,0,0);   // Name, other properties 0.

    G4ThreeVector positionBulkTop1 = G4ThreeVector(0.0,(25.0*mm+1.0*mm),-23.5*mm); // w.r.t. Ion Chamber

    physiBulkTop1 = new G4PVPlacement(0,                // no rotation
				      positionBulkTop1, // position defined above
				      logicBulkTB,    // its logical volume
                                      "BulkTop1",       // its name
				      logicIonChamberGas,    // its mother  volume
				      false,              // no boolean operations
				      1);                 // copy number

    G4ThreeVector positionBulkTop2 = G4ThreeVector(0.0,(25.0*mm+1.0*mm),23.5*mm); // w.r.t. Ion Chamber

    physiBulkTop2 = new G4PVPlacement(0,                // no rotation
				      positionBulkTop2, // position defined above
				      logicBulkTB,    // its logical volume
                                      "BulkTop2",       // its name
				      logicIonChamberGas,    // its mother  volume
				      false,              // no boolean operations
				      2);                 // copy number

    G4ThreeVector positionBulkBottom1 = G4ThreeVector(0.0,-(25.0*mm+1.0*mm),-23.5*mm); // w.r.t. Ion Chamber

    physiBulkBottom1 = new G4PVPlacement(0,                // no rotation
				      positionBulkBottom1, // position defined above
				      logicBulkTB,    // its logical volume
                                      "BulkBottom1",       // its name
				      logicIonChamberGas,    // its mother  volume
				      false,              // no boolean operations
				      1);                 // copy number

    G4ThreeVector positionBulkBottom2 = G4ThreeVector(0.0,-(25.0*mm+1.0*mm),23.5*mm); // w.r.t. Ion Chamber

    physiBulkBottom2 = new G4PVPlacement(0,                // no rotation
				      positionBulkBottom2, // position defined above
				      logicBulkTB,    // its logical volume
                                      "BulkBottom2",       // its name
				      logicIonChamberGas,    // its mother  volume
				      false,              // no boolean operations
				      2);                 // copy number

  // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* BulkTBVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.5,0.6));  //Blue 
  BulkTBVisAtt->SetForceSolid(true);
  logicBulkTB->SetVisAttributes(BulkTBVisAtt);

  G4double BulkFBx = 90.0*mm;
  G4double BulkFBy = 70.0*mm;  // Total active Y length = 70mm, so 70mm-30mm(MicroBulk) = 40.mm/2 (top/bottom)
  G4double BulkFBz = 20.0*mm;
 
    solidBulkFB = new G4Box("BulkFB",BulkFBx/2.,BulkFBy/2.,BulkFBz/2.);
 
    logicBulkFB = new G4LogicalVolume(solidBulkFB,        // Solid
				       ArMethane,               // Material 
				       "BulkFB",0,0,0);   // Name, other properties 0.

    G4ThreeVector positionBulkFront = G4ThreeVector(0.0,0.0,-57.0*mm); // w.r.t. Ion Chamber

    physiBulkFront = new G4PVPlacement(0,                // no rotation
				      positionBulkFront, // position defined above
				      logicBulkFB,    // its logical volume
                                      "BulkFront",       // its name
				      logicIonChamberGas,    // its mother  volume
				      false,              // no boolean operations
				      1);                 // copy number

    G4ThreeVector positionBulkBack = G4ThreeVector(0.0,0.0,57.0*mm); // w.r.t. Ion Chamber

    physiBulkBack = new G4PVPlacement(0,                // no rotation
				      positionBulkBack, // position defined above
				      logicBulkFB,    // its logical volume
                                      "BulkBack",       // its name
				      logicIonChamberGas,    // its mother  volume
				      false,              // no boolean operations
				      2);                 // copy number

 // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* BulkFBVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.5,0.6));  //Blue 
  BulkFBVisAtt->SetForceSolid(true);
  logicBulkFB->SetVisAttributes(BulkFBVisAtt);

  G4double BeamCounterX = 90.0*mm;
  G4double BeamCounterY = 10.0*mm;
  G4double BeamCounterZ = 10.0*mm;
 
    solidBeamCounter = new G4Box("BeamCounter",BeamCounterX/2.,BeamCounterY/2.,BeamCounterZ/2.);
 
    logicBeamCounter = new G4LogicalVolume(solidBeamCounter,        // Solid
				       ArMethane,               // Material 
				       "BeamCounter",0,0,0);   // Name, other properties 0.

    G4ThreeVector positionBeamCounterFront = G4ThreeVector(0.0,0.0,-73.5*mm); // w.r.t. Ion Chamber

    physiBeamCounterFront = new G4PVPlacement(0,                // no rotation
				      positionBeamCounterFront, // position defined above
				      logicBeamCounter,    // its logical volume
                                      "BeamCounterFront",       // its name
				      logicIonChamberGas,    // its mother  volume
				      false,              // no boolean operations
				      1);                 // copy number

    G4ThreeVector positionBeamCounterBack = G4ThreeVector(0.0,0.0,73.5*mm); // w.r.t. Ion Chamber

    physiBeamCounterBack = new G4PVPlacement(0,                // no rotation
				      positionBeamCounterBack, // position defined above
				      logicBeamCounter,    // its logical volume
                                      "BeamCounterBack",       // its name
				      logicIonChamberGas,    // its mother  volume
				      false,              // no boolean operations
				      2);                 // copy number

 // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* BeamCounterVisAtt
    = new G4VisAttributes(G4Colour(0.0,1.0,0.0));  //Green 
  BeamCounterVisAtt->SetForceSolid(true);
  logicBeamCounter->SetVisAttributes(BeamCounterVisAtt);

  
 // Back Si Detector for beam tuning --------------------------------------
  G4double SiDetSize = 3.5*cm;
  G4double SiDetThickness = 0.3*mm;

    solidSiDet = new G4Box("SiDet",SiDetSize/2.,SiDetSize/2.,SiDetThickness/2.);
 
    logicSiDet = new G4LogicalVolume(solidSiDet,        // Solid
				       Si,               // Material 
				       "SiDet",0,0,0);   // Name, other properties 0.

    G4ThreeVector positionSiDet = G4ThreeVector(0.0,0.0,2.*(halfHeightICGasP2+halfHeightWin+SiDetThickness)); // w.r.t. World, 0 at front Window

    physiSiDet = new G4PVPlacement(0,                // no rotation
				      positionSiDet, // position defined above
				      logicSiDet,    // its logical volume
                                      "SiDet",       // its name
				      logicWorld,    // its mother  volume
				      false,              // no boolean operations
				      0);                 // copy number

 // Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* SiDetVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.5,0.0));  
  SiDetVisAtt->SetForceSolid(true);
  logicSiDet->SetVisAttributes(SiDetVisAtt);

  // Sensitive Detector Setup ------------------------------------------------------- 
  
    // Make Proton Detector Sensitive to Heavy Ions
  
    G4SDManager* SDMan = G4SDManager::GetSDMpointer();

    G4String IonChamberSDname = "IonChamberSD";   // To access use "IonChamberSD/IonChamberHitsCollection"
    IonChamberSD* theIonChamberSD = new IonChamberSD(IonChamberSDname);
    SDMan->AddNewDetector( theIonChamberSD );
  
    logicMicroBulk->SetSensitiveDetector(theIonChamberSD);
    logicBulkTB->SetSensitiveDetector(theIonChamberSD);
    logicBulkFB->SetSensitiveDetector(theIonChamberSD);
    logicBeamCounter->SetSensitiveDetector(theIonChamberSD);
    logicSiDet->SetSensitiveDetector(theIonChamberSD);
  

  return physiWorld;    // Returns Detector Construction to RunManager
}


