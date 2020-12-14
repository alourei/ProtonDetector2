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
// $Id: ProtonConstruction.cc,v 0.1 2013/11/27 17:47:13 D. Perez-Loureiro Exp $
// GEANT4 tag $Name: geant4-09-06p02 $
//
//
// --------------------------------------------------------------
// Based on     GEANT 4 - AstroBox2 and ActarSim2013
// --------------------------------------------------------------
//
//

#include <ProtonDetectorConstruction.hh>

#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4UnionSolid.hh>
#include <G4RunManager.hh>

#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include "G4SDManager.hh"

#include "IonChamberSD.hh"
#include "ProtonDetectorSD.hh"
#include "DegraderSD.hh"
#include "VetoSD.hh"
#include <ProtonDetectorConstructionMessenger.hh>


ProtonDetectorConstruction::ProtonDetectorConstruction() {


	 //detectorGeometry="box";
	 detectorGeometry="tube";
	 degraderIncludedFlag = "off";
	 //Default dimensions of gas box
	 XboxLength = 100*mm;
	 YboxLength = 120*mm;
	 ZboxLength = 180*mm;

	 radiusGasTube = 60*mm;
	 lengthGasTube = 320*mm;

	 //innerPadRadius = 0*mm;
	 innerPadRadius = 20*mm;
	 outerPadRadius = 45*mm;
	 //outerPadRadius = 200*mm;

	 degraderThickness = 2.8*mm;
	 degraderAngle = 0*deg;
	 degraderPosition = G4ThreeVector(0,0,-25*cm);



	 theMessenger = new ProtonDetectorConstructionMessenger(this);

	 //Make Sensitive Detectors

	 G4SDManager* SDMan = G4SDManager::GetSDMpointer();

	 G4String IonChamberSDname = "IonChamberSD";   // To access use "IonChamberSD/IonChamberHitsCollection"
	 //theIonChamberSD = new IonChamberSD(IonChamberSDname);
	 //SDMan->AddNewDetector(theIonChamberSD);


	 G4String ProtonDetectorSDName = "ProtonDetectorSD";
	 theProtonDetectorSD = new ProtonDetectorSD(ProtonDetectorSDName);
	 SDMan->AddNewDetector(theProtonDetectorSD);

	 G4String VetoSDName = "VetoDetectorSD";
	 theVetoSD = new VetoSD(VetoSDName);
	 SDMan->AddNewDetector(theVetoSD);

	 G4String DegraderSDName = "DegraderSD";
	 theDegraderSD = new DegraderSD(DegraderSDName);
	 SDMan->AddNewDetector(theDegraderSD);
	 //define materials and set medium material
	 DefineMaterials();


	 SetGasMaterial("ArIso-2.5%_1atm");



}

ProtonDetectorConstruction::~ProtonDetectorConstruction() {

	delete theMessenger;

}

void ProtonDetectorConstruction::DefineMaterials(){

	  G4double a;  // atomic mass
	  G4double z;  // atomic number
	  G4double density, temperature, pressure;
	  G4int nel;   // number of elements in compound


	  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
	  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
	  G4Element* Xe = new G4Element("Xenon"  , "Xe", z=54., a= 131.293*g/mole);
	  G4Element* H  = new G4Element("Hydrogen","H" , z= 1., a= 1.01*g/mole);
	  G4Element* C  = new G4Element("Carbon"  ,"C" , z= 6., a= 12.011*g/mole);

	  G4NistManager* NistMan = G4NistManager::Instance();

	  NistMan->FindOrBuildMaterial("G4_Al");
	  NistMan->FindOrBuildMaterial("G4_KAPTON");
	  NistMan->FindOrBuildMaterial("G4_Fe");

	  G4Material *Ar = NistMan->FindOrBuildMaterial("G4_Ar");
	  G4Material *Ne = NistMan->FindOrBuildMaterial("G4_Ne");
	  G4Material *Methane = NistMan->FindOrBuildMaterial("G4_METHANE");
	  G4Material *isobutane = NistMan->FindOrBuildMaterial("G4_BUTANE");
	  G4Material *CO2 = NistMan->FindOrBuildMaterial("G4_CARBON_DIOXIDE");

	  G4Material *mylar = NistMan->FindOrBuildMaterial("G4_MYLAR");
	  G4Material *Al = NistMan->FindOrBuildMaterial("G4_Al");

	  //Aluminized mylar
	  G4Material *Almylar = new G4Material("Almylar",density=0.999999*g/cm3,nel=2);
	  Almylar->AddMaterial(Al,0.2);
	  Almylar->AddMaterial(mylar,0.8);
	  //Vacuum
	  G4Material* Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
			  	  	  	  	  	  	  	  kStateGas, 3.e-18*pascal, 2.73*kelvin);
	  Vacuum->SetChemicalFormula("NOMATTER");

	  //Air
	   density = 1.205*mg/cm3;
	   G4Material *Air = new G4Material("Air", density, nel=2);
	   Air->AddElement(N, .7);
	   Air->AddElement(O, .3);
	  

	   //Xenon 

	   G4double PressureFactor = 1.;  // 760 torr
	   density = PressureFactor*0.0054586*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *Xegas_1atm = new G4Material("Xenon_1atm", density, nel=1, kStateGas,
					      temperature= 293.15*kelvin, pressure);
	   Xegas_1atm->AddElement(Xe,1);

	   PressureFactor = 1.05263;  // 800 torr
	   density = PressureFactor*0.0054586*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *Xegas_800torr = new G4Material("Xenon_800torr", density, nel=1, kStateGas,
					      temperature= 293.15*kelvin, pressure);
	   Xegas_800torr->AddElement(Xe,1);

	   PressureFactor = 2;  //2 atm
	   density = PressureFactor*0.0054586*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *Xegas_2atm = new G4Material("Xenon_2atm", density, nel=1, kStateGas,
					      temperature= 293.15*kelvin, pressure);
	   Xegas_2atm->AddElement(Xe,1);

	   PressureFactor = 5.0000;  // 5.0 atm
	   density = PressureFactor*0.0054586*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *Xegas_5atm = new G4Material("Xenon_5atm", density, nel=1, kStateGas,
					      temperature= 293.15*kelvin, pressure);

	   Xegas_5atm->AddElement(Xe,1);

	   //P5 Gas 95%Ar 5%Methane

	   PressureFactor = 0.5;  // 0.5 atm
	   density = PressureFactor*0.001611091*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArMethane_05atm = new G4Material("P5_0.5atm", density, nel=2, kStateGas,
						  temperature= 293.15*kelvin, pressure);
	   ArMethane_05atm->AddMaterial(Ar, 0.95);
	   ArMethane_05atm->AddMaterial(Methane, 0.05);

	   PressureFactor = 1.;  // 0.5 atm
	   density = PressureFactor*0.001611091*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArMethane_1atm = new G4Material("P5_1atm", density, nel=2, kStateGas,
						  temperature= 293.15*kelvin, pressure);
	   ArMethane_1atm->AddMaterial(Ar, 0.95);
	   ArMethane_1atm->AddMaterial(Methane, 0.05);

	   PressureFactor = 1.05263;  // 800 torr
	   density = PressureFactor*0.001611091*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArMethane_800torr = new G4Material("P5_800torr", density, nel=2, kStateGas,
							  temperature= 293.15*kelvin, pressure);
	   ArMethane_800torr->AddMaterial(Ar, 0.95);
	   ArMethane_800torr->AddMaterial(Methane, 0.05);

	   PressureFactor = 1.5;  // 1.5 atm
	   density = PressureFactor*0.001611091*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArMethane_15atm = new G4Material("P5_1.5atm", density, nel=2, kStateGas,
							  temperature= 293.15*kelvin, pressure);
	   ArMethane_15atm->AddMaterial(Ar, 0.95);
	   ArMethane_15atm->AddMaterial(Methane, 0.05);



	   //P10 Gas 90%Ar 10%Methane

	   PressureFactor = 0.5;  // 0.5 atm
	   density = PressureFactor*0.0015613*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *P10_05atm = new G4Material("P10_0.5atm", density, nel=2, kStateGas,
						  temperature= 293.15*kelvin, pressure);
	   P10_05atm->AddMaterial(Ar, 0.95);
	   P10_05atm->AddMaterial(Methane, 0.05);

	   PressureFactor = 1.;  // 0.5 atm
	   density = PressureFactor*0.0015613*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *P10_1atm = new G4Material("P10_1atm", density, nel=2, kStateGas,
						  temperature= 293.15*kelvin, pressure);
	   P10_1atm->AddMaterial(Ar, 0.95);
	   P10_1atm->AddMaterial(Methane, 0.05);

	   PressureFactor = 1.05263;  // 800 torr
	   density = PressureFactor*0.0015613*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *P10_800torr = new G4Material("P10_800torr", density, nel=2, kStateGas,
							  temperature= 293.15*kelvin, pressure);
	   P10_800torr->AddMaterial(Ar, 0.95);
	   P10_800torr->AddMaterial(Methane, 0.05);

	   PressureFactor = 1.5;  // 1.5 atm
	   density = PressureFactor*0.0015613*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *P10_15atm = new G4Material("P10_1.5atm", density, nel=2, kStateGas,
							  temperature= 293.15*kelvin, pressure);
	   P10_15atm->AddMaterial(Ar, 0.95);
	   P10_15atm->AddMaterial(Methane, 0.05);



	   //Gas 97.5%Ar 2.5%isoButane
	   PressureFactor = 1;  // 760 torr
	   density = PressureFactor*0.0016807*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_25_1atm = new G4Material("ArIso-2.5%_1atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_25_1atm->AddMaterial(Ar, 0.975);
	   ArisoButane_25_1atm->AddMaterial(isobutane, 0.025);


	   PressureFactor = 1.05263;  // 800 torr
	   density = PressureFactor*0.0016807*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_25_800torr = new G4Material("ArIso-2.5%_800torr", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_25_800torr->AddMaterial(Ar, 0.975);
	   ArisoButane_25_800torr->AddMaterial(isobutane, 0.025);

	   PressureFactor = 1.5;  // 1.5 atm
	   density = PressureFactor*0.0016807*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_25_15atm = new G4Material("ArIso-2.5%_1.5atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_25_15atm->AddMaterial(Ar, 0.975);
	   ArisoButane_25_15atm->AddMaterial(isobutane, 0.025);


	   PressureFactor = 0.5;  // 1.5 atm
	   density = PressureFactor*0.0016807*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_25_05atm = new G4Material("ArIso-2.5%_0.5atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_25_05atm->AddMaterial(Ar, 0.975);
	   ArisoButane_25_05atm->AddMaterial(isobutane, 0.025);


	   //Gas 95%Ar 5%isoButane
	   PressureFactor = 1;  // 760 torr
	   density = PressureFactor*0.0016985*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_5_1atm = new G4Material("ArIso-5%_1atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_5_1atm->AddMaterial(Ar, 0.95);
	   ArisoButane_5_1atm->AddMaterial(isobutane, 0.05);


	   PressureFactor = 1.05263;  // 800 torr
	   density = PressureFactor*0.0016985*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_5_800torr = new G4Material("ArIso-5%_800torr", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_5_800torr->AddMaterial(Ar, 0.95);
	   ArisoButane_5_800torr->AddMaterial(isobutane, 0.05);

	   PressureFactor = 1.5;  // 1.5 atm
	   density = PressureFactor*0.0016985*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_5_15atm = new G4Material("ArIso-5%_1.5atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_5_15atm->AddMaterial(Ar, 0.95);
	   ArisoButane_5_15atm->AddMaterial(isobutane, 0.05);


	   PressureFactor = 0.5;  // 0.5 atm
	   density = PressureFactor*0.0016985*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_5_05atm = new G4Material("ArIso-5%_0.5atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_5_05atm->AddMaterial(Ar, 0.95);
	   ArisoButane_5_05atm->AddMaterial(isobutane, 0.05);

	   PressureFactor = 2;  // 2 atm
	   density = PressureFactor*0.0016985*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_5_2atm = new G4Material("ArIso-5%_2atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_5_2atm->AddMaterial(Ar, 0.95);
	   ArisoButane_5_2atm->AddMaterial(isobutane, 0.05);


	   //Gas 90%Ar 10%isoButane
	   PressureFactor = 1;  // 760 torr
	   density = PressureFactor*0.0017362*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_10_1atm = new G4Material("ArIso-10%_1atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_10_1atm->AddMaterial(Ar, 0.9);
	   ArisoButane_10_1atm->AddMaterial(isobutane, 0.1);


	   PressureFactor = 1.05263;  // 800 torr
	   density = PressureFactor*0.0017362*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_10_800torr = new G4Material("ArIso-10%_800torr", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_10_800torr->AddMaterial(Ar, 0.9);
	   ArisoButane_10_800torr->AddMaterial(isobutane, 0.1);

	   PressureFactor = 1.5;  // 1.5 atm
	   density = PressureFactor*0.0017362*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_10_15atm = new G4Material("ArIso-10%_1.5atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_10_15atm->AddMaterial(Ar, 0.9);
	   ArisoButane_10_15atm->AddMaterial(isobutane, 0.1);


	   PressureFactor = 0.5;  // 1.5 atm
	   density = PressureFactor*0.0017362*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArisoButane_10_05atm = new G4Material("ArIso-10%_0.5atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArisoButane_10_05atm->AddMaterial(Ar, 0.9);
	   ArisoButane_10_05atm->AddMaterial(isobutane, 0.1);



	   //Gas 90%Ar 10%CO2

	   PressureFactor = 1;  // 760 torr
	   density = PressureFactor*0.0015945*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArCO2_10_1atm = new G4Material("ArCO2-10%_1atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArCO2_10_1atm->AddMaterial(Ar, 0.9);
	   ArCO2_10_1atm->AddMaterial(CO2, 0.1);


	   PressureFactor = 1.05263;  // 800 torr
	   density = PressureFactor*0.0015945*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArCO2_10_800torr = new G4Material("ArCO2-10%_800torr", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArCO2_10_800torr->AddMaterial(Ar, 0.9);
	   ArCO2_10_800torr->AddMaterial(CO2, 0.1);

	   PressureFactor = 1.5;  // 1.5 atm
	   density = PressureFactor*0.0015945*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArCO2_10_15atm = new G4Material("ArCO2-10%_1.5atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArCO2_10_15atm->AddMaterial(Ar, 0.9);
	   ArCO2_10_15atm->AddMaterial(CO2, 0.1);


	   PressureFactor = 0.5;  // 0.5 atm
	   density = PressureFactor*0.0015945*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *ArCO2_10_05atm = new G4Material("ArCO2-10%_0.5atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   ArCO2_10_05atm->AddMaterial(Ar, 0.9);
	   ArCO2_10_05atm->AddMaterial(CO2, 0.1);

	   //Gas 95%Ne 5%isoButane
	   PressureFactor = 1;  // 760 torr
	   density = PressureFactor*9.184e-4*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *NeisoButane_5_1atm = new G4Material("NeIso-5%_1atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   NeisoButane_5_1atm->AddMaterial(Ne, 0.95);
	   NeisoButane_5_1atm->AddMaterial(isobutane, 0.05);

	   PressureFactor = 2;  // 2 atm
	   density = PressureFactor*9.184e-4*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *NeisoButane_5_2atm = new G4Material("NeIso-5%_2atm", density, nel=2, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   NeisoButane_5_2atm->AddMaterial(Ne, 0.95);
	   NeisoButane_5_2atm->AddMaterial(isobutane, 0.05);

	   //Pure Ne

	   PressureFactor = 1.;  // 1 atm
	   density = PressureFactor*0.0008385*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *Ne_1atm = new G4Material("Ne_1atm", density, nel=1, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   Ne_1atm->AddMaterial(Ne,1);

	   PressureFactor = 1.;  // 1 atm
	   density = PressureFactor*0.001662*g/cm3;
	   pressure = PressureFactor*atmosphere;
	   G4Material *Ar_1atm = new G4Material("Ar_1atm", density, nel=1, kStateGas,
						    temperature= 293.15*kelvin, pressure);
	   Ar_1atm->AddMaterial(Ar,1);


	   //Silicon (for Silicon Detector)
	   density = 2.33*g/cm3;
	   G4Material* Si = new G4Material("Si", z=14., a=28.0855*g/mole, density);

	   if (Si){;}

	   // Print all the materials defined.
	   G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	   G4cout << *(G4Material::GetMaterialTable()) << G4endl;



}

G4VPhysicalVolume* ProtonDetectorConstruction::Construct(){

return ConstructDetector();

}

G4VPhysicalVolume* ProtonDetectorConstruction::ConstructDetector(){

	G4double WorldSize_X=5*m;
	G4double WorldSize_Y=5*m;
	G4double WorldSize_Z=5*m;

	worldVolume = new G4Box("World",WorldSize_X/2.,WorldSize_Y/2.,WorldSize_Z/2.);

	G4Material *pttoMaterial = G4Material::GetMaterial("Galactic");


	world_log = new G4LogicalVolume(worldVolume,pttoMaterial,"World_log");

	world_log->SetVisAttributes (G4VisAttributes::Invisible);

	world_phys = new G4PVPlacement(0,G4ThreeVector(),world_log,"World",0,false,0);




	if(detectorGeometry == "box"){

		G4cout<< "##################################################################" <<G4endl
				<<"######  ProtonDetectorConstruction::ConstructDetector()  ########" <<G4endl
				<< " Box-like gas geometry" << G4endl;
		G4cout<< "##################################################################"<<G4endl;

		//Constructing the gas chamber

		G4double gasChamberPart1_X = XboxLength/2.;
		G4double gasChamberPart1_Y = YboxLength/2.;
		G4double gasChamberPart1_Z = ZboxLength/2.;


		G4VSolid *GasChamberPart1 =new G4Box("Chamber1",gasChamberPart1_X,gasChamberPart1_Y,gasChamberPart1_Z);

		G4double gasChamberPart2_Radius = (85.0/2.)*mm;
		G4double gasChamberPart2_Length = (ZboxLength+10)/2.*mm;
	    G4double StartAngle = 0.*deg;
	    G4double SpanAngle = 360.*deg;

	    G4VSolid *GasChamberPart2 = new G4Tubs("Chamber2",0.*mm,gasChamberPart2_Radius,gasChamberPart2_Length,StartAngle,SpanAngle);


	    G4VSolid *GasChamber = new G4UnionSolid("Chamber",GasChamberPart1,GasChamberPart2);

	    pttoMaterial= G4Material::GetMaterial("G4_Fe");

	    G4VisAttributes *gasCVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

	    gasChamber_log = new G4LogicalVolume(GasChamber,pttoMaterial,"chamber_log");

	    gasChamber_log->SetVisAttributes(gasCVisAtt);

	    gasChamber_phys= new G4PVPlacement(0,G4ThreeVector(),gasChamber_log,"Chamber",world_log,false,0);

	    //Construct the Gas volume

		G4Box *gasBox = new G4Box("GasBox",(XboxLength-5)/2.,(YboxLength-5)/2.,(ZboxLength-5)/2.);

		G4double gasTube_Radius = (80.0/2.)*mm;
		//G4double gasTube_Length = (ZboxLength+50)/2.*mm;
		G4double gasTube_Length = (ZboxLength+10)/2.*mm;
	    StartAngle = 0.*deg;
	    SpanAngle = 360.*deg;

		G4Tubs *gasTube = new G4Tubs("GasTube",0*mm,gasTube_Radius,gasTube_Length,StartAngle,SpanAngle);


		G4VSolid * gasVolume = new G4UnionSolid("GasVolume",gasBox,gasTube);
		//G4VSolid * gasVolume = gasBox;

		//pttoMaterial= G4Material::GetMaterial("ArMethane");
		//pttoMaterial= G4Material::GetMaterial("Xegas");

		G4VisAttributes *gasBoxVisAtt = new G4VisAttributes(G4Colour(1.,1.0,0.0));

		gasVolume_log = new G4LogicalVolume(gasVolume,gasMaterial,"GasBox_log");

		gasVolume_log->SetVisAttributes(gasBoxVisAtt);

		G4double XboxPos = 0*mm;
		G4double YboxPos = 0*mm;
		G4double ZboxPos = 0*mm;

		//gasVolume_phys = new G4PVPlacement(0,G4ThreeVector(XboxPos,YboxPos,ZboxPos),gasVolume_log,"Gas",gasChamber_log,false,0);
		gasVolume_phys = new G4PVPlacement(0,G4ThreeVector(XboxPos,YboxPos,ZboxPos),gasVolume_log,"Gas",gasChamber_log,false,0);

		//Construct the windows

		G4double ICRminWin = 0.*mm;
		G4double ICRmaxWin = (80.0/2.)*mm;  //diameter/2.
		G4double halfHeightWin = (0.0508/2.)*mm; // 2 mils thick
		StartAngle = 0.*deg;
		SpanAngle = 360.*deg;

		G4Material *Kapton = G4Material:: GetMaterial("G4_KAPTON");

		G4Tubs* solidWindow = new G4Tubs("KaptonWindow",ICRminWin,ICRmaxWin,halfHeightWin,StartAngle,SpanAngle);

		G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow,        // Solid
						       Kapton,            // Material
						       "KaptonWindow",0,0,0);   // Name, other properties 0.

		G4ThreeVector positionWindowFront = G4ThreeVector(0.0,0.0,-(gasTube_Length-halfHeightWin)); // w.r.t. Ion Chamber

		G4RotationMatrix* WindowRot = new G4RotationMatrix();
		WindowRot->rotateY(0.*deg);

		G4VPhysicalVolume* physiWindowFront     = new G4PVPlacement(WindowRot,                // no rotation
							 positionWindowFront, // position defined above
							 logicWindow,    // its logical volume
							 "KaptonWindowFront",       // its name
							 gasVolume_log,      // its mother  volume
							 false,              // no boolean operations
							 0);                 // copy number

		    G4ThreeVector positionWindowBack = G4ThreeVector(0.0,0.0,(gasTube_Length-halfHeightWin)); // w.r.t. Ion Chamber

		G4VPhysicalVolume* physiWindowBack     = new G4PVPlacement(WindowRot,                // no rotation
							 positionWindowBack, // position defined above
							 logicWindow,    // its logical volume
							 "KaptonWindowBack",       // its name
							 gasVolume_log,      // its mother  volume
							 false,              // no boolean operations
							 1);                 // copy number

		  G4VisAttributes* WindowVisAtt
		    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));  //Blue
		  WindowVisAtt->SetForceSolid(true);
		  logicWindow->SetVisAttributes(WindowVisAtt);


	}

	else if(detectorGeometry == "tube"){

	  G4cout<< "##################################################################" <<G4endl
		<<"######  ProtonDetectorConstruction::ConstructDetector()  ########" <<G4endl
		<< " Tube-like gas geometry" << G4endl;
	  G4cout<< "##################################################################"<<G4endl;
	  
	  //Constructing the gas chamber
	  
	  G4double gasChamberPart1_Radius = radiusGasTube;
	  G4double gasChamberPart1_Length = lengthGasTube/2;
	  G4double StartAngle = 0.*deg;
	  G4double SpanAngle = 360.*deg;
	  

	  
	  G4VSolid *GasChamberPart1 =new G4Tubs("Chamber1",0.*mm,gasChamberPart1_Radius,gasChamberPart1_Length,StartAngle,SpanAngle);
	  
	  G4double gasChamberPart2_Radius0 = (radiusGasTube/2-3.)*mm;
	  G4double gasChamberPart2_Radius = (radiusGasTube/2.)*mm;
	  G4double gasChamberPart2_Length = (lengthGasTube)/2.*mm;
	  StartAngle = 0.*deg;
	  SpanAngle = 360.*deg;
	  
	  G4VSolid *GasChamberPart2 = new G4Tubs("Chamber2",0.*mm,gasChamberPart2_Radius,gasChamberPart2_Length,StartAngle,SpanAngle);
	  
	  G4ThreeVector pos(0,0,-.5*cm);
	  //G4ThreeVector pos(0,0,0*cm);
	  
	  G4VSolid *GasChamber = new G4UnionSolid("Chamber",GasChamberPart1,GasChamberPart2,0,pos);
	  
	  pttoMaterial= G4Material::GetMaterial("G4_Al");
	  
	  G4VisAttributes *gasCVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  
	  gasChamber_log = new G4LogicalVolume(GasChamber,pttoMaterial,"chamber_log");
	  
	  gasChamber_log->SetVisAttributes(gasCVisAtt);
	  
	  //gasChamber_phys= new G4PVPlacement(0,G4ThreeVector(),gasChamber_log,"Chamber",world_log,false,0);
	  
	  //DPL Apr2016 new dimensions for the gas Volume

	  G4double XboxPos = 0*mm;
	  G4double YboxPos = 0*mm;
	  G4double ZboxPos = 0*mm;
	  G4double gasTube_Length = (lengthGasTube-2)/2.*mm;
	  if( innerPadRadius>0){
	  
	  G4double gasTube_Radius = innerPadRadius;
	  gasTube_Length = (lengthGasTube-2)/2.*mm;
	  StartAngle = 0.*deg;
	  SpanAngle = 360.*deg;

	  G4Tubs *gasTube1 = new G4Tubs("InnerPad",0,gasTube_Radius,gasTube_Length,StartAngle,SpanAngle);

	  InnerPad_log = new G4LogicalVolume(gasTube1,gasMaterial,"InnerPad_log");

	  XboxPos = 0*mm;
	  YboxPos = 0*mm;
	  ZboxPos = 0*mm;

	  InnerPad_phys= new G4PVPlacement(0,G4ThreeVector(XboxPos,YboxPos,ZboxPos),InnerPad_log,"InnerPad",world_log,false,0);
	  }

	  G4double gasTube_Radius2 = outerPadRadius;
	  gasTube_Length = (lengthGasTube-2)/2.*mm;
	  //StartAngle = 0.*deg;
	  //SpanAngle = 360.*deg;
	  StartAngle = 0.*deg;
	  Int_t ndivisions=4;
	  SpanAngle = 360./ndivisions *deg;

	  G4Tubs *gasTube2 = new G4Tubs("OuterPad",innerPadRadius,gasTube_Radius2,gasTube_Length,StartAngle,SpanAngle);

	  OuterPad_log = new G4LogicalVolume(gasTube2,gasMaterial,"OuterPad_log");
	  
	  //Divide sectors
	  //StartAngle = 0.*deg;
	  //SpanAngle = 360./ndivisions *deg;

	  //G4Tubs *Sector=new G4Tubs("SectorPad",innerPadRadius,gasTube_Radius2,gasTube_Length,StartAngle,SpanAngle);

	  //Sector_log=new G4LogicalVolume(Sector,gasMaterial,"Sector_log");

	  XboxPos = 0*mm;
	  YboxPos = 0*mm;
	  ZboxPos = 0*mm;
	  //OuterPad_phys= new G4PVPlacement(0,G4ThreeVector(XboxPos,YboxPos,ZboxPos),OuterPad_log,"OuterPad",world_log,false,0);
	  
	  for(Int_t k=0;k<ndivisions;k++){
	  
	    G4RotationMatrix rotm  = G4RotationMatrix();
	    rotm.rotateZ((G4double)k/ndivisions*360*deg);
	    G4ThreeVector position(XboxPos,YboxPos,ZboxPos) ;
	    G4Transform3D transform = G4Transform3D(rotm,position);
	    
	    OuterPad_phys= new G4PVPlacement(transform,OuterPad_log,"OuterPad",world_log,false,k);
	  //Pad_phys=new G4PVReplica("PhiSlices",Sector_log,OuterPad_log,kPhi,ndivisions,SpanAngle, 0);
	  }

	  StartAngle = 0.*deg;
	  ndivisions=8;
	  SpanAngle = 360./ndivisions*deg;

	  Double_t Veto_radius=outerPadRadius+9*mm;
	  	  
	  
	  G4Tubs *gasTube3 = new G4Tubs("Veto",outerPadRadius,Veto_radius,gasTube_Length,StartAngle,SpanAngle);
	  
	  //G4ThreeVector pos(0,0,-2*cm);
	  
	  //G4VSolid * gasVolume = new G4UnionSolid("GasVolume",gasTube1,gasTube2,0,pos);
	  
	  G4cout<<"Selected Gas "<<gasMaterial->GetName()<<G4endl;
	  
	  G4VisAttributes *gasBoxVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  G4VisAttributes *vetoBoxVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  
	  
	  Veto_log = new G4LogicalVolume(gasTube3,gasMaterial,"Veto_log");
	  
	  OuterPad_log->SetVisAttributes(gasBoxVisAtt);
	  Veto_log->SetVisAttributes(vetoBoxVisAtt);
	  
	   XboxPos = 0*mm;
	   YboxPos = 0*mm;
	   ZboxPos = 0*mm;

	  //gasVolume_phys = new G4PVPlacement(0,G4ThreeVector(XboxPos,YboxPos,ZboxPos),gasVolume_log,"Gas",gasChamber_log,false,0);
	  //gasVolume_phys= new G4PVPlacement(0,G4ThreeVector(XboxPos,YboxPos,ZboxPos),gasVolume_log,"Gas",world_log,false,0);
	  
	  //gasVolume_phys=new G4PVReplica("PhiSlices",gasVolume_log,world_log,kPhi,4, M_PI*0.5*rad, 0);
	  

	  for(Int_t k=0;k<ndivisions;k++){
	  
	    G4RotationMatrix rotm  = G4RotationMatrix();
	    rotm.rotateZ((G4double)k/ndivisions*360*deg);
	    G4ThreeVector position(XboxPos,YboxPos,ZboxPos) ;
	    G4Transform3D transform = G4Transform3D(rotm,position);

	    Veto_phys = new G4PVPlacement(transform,Veto_log,"Veto",world_log,false,k);
	  } 
	  //Construct the windows
	  
	  G4double ICRminWin = 0.*mm;
	  G4double ICRmaxWin = (radiusGasTube/2.)*mm;  //diameter/2.
	  G4double halfHeightWin = (0.0508/2.)*mm; // 2 mils thick
	  StartAngle = 0.*deg;
	  SpanAngle = 360.*deg;
	  
	  G4Material *Kapton = G4Material:: GetMaterial("G4_KAPTON");

	  G4Tubs* solidWindow = new G4Tubs("KaptonWindow",ICRminWin,ICRmaxWin,halfHeightWin,StartAngle,SpanAngle);
	  
	  G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow,        // Solid
							     Kapton,            // Material
							     "KaptonWindow",0,0,0);   // Name, other properties 0.
	  
	  G4ThreeVector positionWindowFront = G4ThreeVector(0.0,0.0,-.5*cm-(gasTube_Length-halfHeightWin)); // w.r.t. Ion Chamber
	  
	  G4RotationMatrix* WindowRot = new G4RotationMatrix();
	  WindowRot->rotateY(0.*deg);
	  
	  //G4VPhysicalVolume* physiWindowFront     = new G4PVPlacement(WindowRot,                // no rotation
	  //					 positionWindowFront, // position defined above
	  //					 logicWindow,    // its logical volume
	  //					 "KaptonWindowFront",       // its name
	  //					 gasVolume_log,      // its mother  volume
	  //					 false,              // no boolean operations
	  //					 0);                 // copy number
	  
	  G4ThreeVector positionWindowBack = G4ThreeVector(0.0,0.0,(gasTube_Length-halfHeightWin)); // w.r.t. Ion Chamber
	  
	  //G4VPhysicalVolume* physiWindowBack     = new G4PVPlacement(WindowRot,                // no rotation
	  //					 positionWindowBack, // position defined above
	  //					 logicWindow,    // its logical volume
	  ///"Kap/tonWindowBack",       // its name
	  //				 gasVolume_log,      // its mother  volume
	  //					 false,              // no boolean operations
	  //				 1);                 // copy number
	  
	  G4VisAttributes* WindowVisAtt
	    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));  //Blue
	  WindowVisAtt->SetForceSolid(true);
	  logicWindow->SetVisAttributes(WindowVisAtt);
	  
	  //Adding Aluminized mylar cathode for the  energy loss
	  
	  G4double cathode_Radius = radiusGasTube-15*mm;
	  G4double cathode_thickness= 3e-3*mm;
	  pttoMaterial= G4Material::GetMaterial("Almylar");
	  
	  //G4Tubs *solidCathode=new G4Tubs("cathode",0*mm,cathode_Radius,cathode_thickness,StartAngle,SpanAngle);
	  //G4LogicalVolume *logicCathode=new G4LogicalVolume(solidCathode,pttoMaterial,"Cathode",0,0,0); 
	  //G4ThreeVector positionCathode = G4ThreeVector(0.0,0.0,2*cm-(gasTube_Length)); // w.r.t. Ion Chamber
	  
	  G4VisAttributes* CathodeVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));  //Red
	  CathodeVisAtt->SetForceSolid(true);
	  //logicCathode->SetVisAttributes(CathodeVisAtt);
	  
	  //G4VPhysicalVolume *physiCathode = new G4PVPlacement(0, positionCathode,
	  //							 logicCathode,"Cathode",gasVolume_log,false,0);
	}
	
	if(degraderIncludedFlag =="on"){
	  
	  G4double DegraderSize = 8.*cm;
	  
	  G4VSolid *solidDegrader = new G4Box("AlDegrader",DegraderSize/2.,DegraderSize/2.,degraderThickness/2.);
	  
	  pttoMaterial= G4Material::GetMaterial("G4_Al");
	  
	  G4LogicalVolume *logicDegrader = new G4LogicalVolume(solidDegrader,pttoMaterial,"AlDegrader",0,0,0);
	  
	  G4ThreeVector positionDegrader = degraderPosition;
	  
	  G4RotationMatrix* DegRotation = new G4RotationMatrix();
	  
	  DegRotation->rotateY(degraderAngle);
	  
	  G4VisAttributes* DegraderVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));  //Red
	  DegraderVisAtt->SetForceSolid(true);
	  logicDegrader->SetVisAttributes(DegraderVisAtt);
	  
	  
	  G4VPhysicalVolume *physiDegrader = new G4PVPlacement(DegRotation, positionDegrader,
							       logicDegrader,"AlDegrader",world_log,false,0);
	  
	  logicDegrader->SetSensitiveDetector(theDegraderSD);
	  
	}
	
	
	//gasVolume_log->SetSensitiveDetector(theIonChamberSD);
	
	if(innerPadRadius>0){
	  InnerPad_log->SetSensitiveDetector(theProtonDetectorSD);
	}
	OuterPad_log->SetSensitiveDetector(theProtonDetectorSD);
	Veto_log->SetSensitiveDetector(theVetoSD);
	
	
	return world_phys;


}

void ProtonDetectorConstruction::UpdateGeometry(){

	 // Updates Proton Detector
	  //
	//Construct();
	G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());

}

void ProtonDetectorConstruction::PrintDetectorParameters(){

	 //
	  // Prints Gas volume detector parameters. To be filled
	  //

	  G4cout << "##################################################################"
		 << G4endl
		 << "##  ActarSimGasDetectorConstruction::PrintDetectorParameters() ###"
		 << G4endl
		 << " The gas volume is a " ;
	  if(detectorGeometry == "box")
	    G4cout << "box; its parameters are:" << G4endl;
	  if(detectorGeometry == "tube")
	    G4cout << "tube; its parameters are:" << G4endl;
//	  G4cout << " The gas material is: " << gasMaterial  << G4endl;
	  if(detectorGeometry == "box")
	    G4cout << " The gasBox size is : " << XboxLength/cm << "x" << YboxLength/cm
		   << "x" << XboxLength/cm << " cm3 " << G4endl << G4endl ;
	  if(detectorGeometry == "tube")
	    G4cout << " The gasTube parameters are: " << G4endl
		   << " radiusGasTube = " <<  radiusGasTube
		   << ",  lengthGasTube = " <<  lengthGasTube << G4endl ;

	  G4cout << "##################################################################"
		 << G4endl;


}



void ProtonDetectorConstruction::SetGasMaterial (G4String mat) {
  //
  // Sets the material the gas is made of
  //
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) gasMaterial = pttoMaterial;
}
