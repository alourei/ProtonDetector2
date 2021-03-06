/*
 * ProtonDetectorPhysicsList.cc
 *
 *  Created on: Jan 30, 2014
 *      Author: perezlou
 */

#include <ProtonDetectorPhysicsList.hh>

#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

// Particle Constructors

#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

// Home-Made Processes

#include "General_1pDecay.hh"
#include "BTRBetaPlus_CustomDecayAl23.hh"

#include "DPLBetaPlusDecay.hh"
#include "DPL1pdecay.hh"
#include "DPLalphadecay.hh"


  // Hadronic/Nucleus ion EM Proton Processes
#include "G4EmProcessOptions.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4eBremsstrahlung.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

// Strong force process Models ("Low" Energy only)
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

// Hadron Processes (Add models above to these processes)
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

// Proton only processes
#include "G4ProtonInelasticProcess.hh"
#include "G4LEProtonInelastic.hh"

#include "G4RadioactiveDecay.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsListHelper.hh"
#include "G4NucleusLimits.hh"


ProtonDetectorPhysicsList::ProtonDetectorPhysicsList() {
;
}

ProtonDetectorPhysicsList::~ProtonDetectorPhysicsList() {
;
}

void ProtonDetectorPhysicsList::ConstructParticle(){

	// In this method, static member functions should be called
	// for all particles which you want to use.
	// This ensures that objects of these particle types will be
	// created in the program.

	G4Geantino::GeantinoDefinition();

	//Gammas and Photons

	G4Gamma::GammaDefinition();

	//Leptons
	G4LeptonConstructor pLeptons;
	pLeptons.ConstructParticle();

	//Baryons
	G4BaryonConstructor pBaryons;
	pBaryons.ConstructParticle();

	//Construct Generic Ions

	// Construct all light nuclei!

	G4Deuteron::DeuteronDefinition();
	G4Triton::TritonDefinition();
	G4Alpha::AlphaDefinition();
	G4He3::He3Definition();

	G4IonConstructor pIons;
	pIons.ConstructParticle();

	//mesons - construct all mesons (strong interactions)
	G4MesonConstructor pMesons;
	pMesons.ConstructParticle();

}

void ProtonDetectorPhysicsList::ConstructProcess(){
  // Define physics processes

  AddTransportation();
  ConstructEM();
  ConstructNuclear();


  //G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

}

void ProtonDetectorPhysicsList::ConstructEM(){

	// Initialize process manager

	  G4ProcessManager *pManager = 0;

	  // The Definitions below taken from DemonEMPhysics.cc
	  // Gamma physics
	  pManager = G4Gamma::Gamma()->GetProcessManager();
	  pManager->AddDiscreteProcess(&thePhotoEffect);
	  pManager->AddDiscreteProcess(&theComptonEffect);
	  pManager->AddDiscreteProcess(&thePairProduction);

	  // Electron Physics
	  pManager = G4Electron::Electron()->GetProcessManager();

	  pManager->AddProcess(&theElectronMultipleScattering, -1, 1, 1);
	  pManager->AddProcess(&theElectronIonisation,         -1, 2, 2);
	  pManager->AddProcess(&theElectronBremsStrahlung,     -1, 3, 3);


	  //Positron Physics
	  pManager = G4Positron::Positron()->GetProcessManager();

	  pManager->AddProcess(&thePositronMultipleScattering, -1, 1, 1);
	  pManager->AddProcess(&thePositronIonisation,         -1, 2, 2);
	  pManager->AddProcess(&thePositronBremsStrahlung,     -1, 3, 3);
	  pManager->AddProcess(&theAnnihilation,                0,-1, 4);

}

void ProtonDetectorPhysicsList::ConstructNuclear(){

	  // Proton nuclear scattering processes

	  G4ProcessManager *pManager = 0;

	  pManager = G4Proton::Proton()->GetProcessManager();

	  // Add proton EM processes
	  pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	  pManager->AddProcess(new G4hIonisation,       -1,2,2);


	  // ////// Nuclei /////////////////////////////////////////////////

	  // Add Deuteron Physics (EM Stopping Power only)
	  	pManager = G4Deuteron::Deuteron()->GetProcessManager();

	          pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	  	pManager->AddProcess(new G4ionIonisation,       -1,2,2);

	  // Add Triton Physics (EM Stopping Power only)
	  	pManager = G4Triton::Triton()->GetProcessManager();

	          pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	  	pManager->AddProcess(new G4ionIonisation,       -1,2,2);

	  // Add Helium 3 Physics (EM Stopping Power only)
	  	pManager = G4He3::He3()->GetProcessManager();

	          pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	  	pManager->AddProcess(new G4ionIonisation,       -1,2,2);


	  // Add Alpha Physics (EM Stopping Power only)
	  	pManager = G4Alpha::Alpha()->GetProcessManager();

	  	pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	  	pManager->AddProcess(new G4ionIonisation, -1,2,2);

	   // Add Generic Ion Physics (EM and Nuclear Elastic, Stopping Power, etc.)
	       pManager = G4GenericIon::GenericIon()->GetProcessManager();

	   // Add EM Processes to Generic Ions

	       pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	       pManager->AddProcess(new G4ionIonisation,       -1,2,2);

	       // To add BetaDecay, use SetParentNucleus() method

	       // Al23 Beta-Decay ----------------------------------------------------
	       // 2/17/2010 -> Qvalue here is g.s. to g.s. as given in ENSDF.
	       // correction due to electron-positron pair production is
	       // taken in account in the process (see PostStepDoIt)

	       //G4String BetaDecayAl23ProcessName = "BTRBetaPlusDecayAl23";
	       //BTRBetaPlus_CustomDecayAl23* theAl23BetaPlusDecay = new BTRBetaPlus_CustomDecayAl23(BetaDecayAl23ProcessName);
	       //theAl23BetaPlusDecay->SetParentNucleus(23,13);
	       //theAl23BetaPlusDecay->SetBetaDecayDataAl23("Al23_BetaIntensity.dat","Branch",12.243*MeV,470.*ms);
	       //pManager->AddRestProcess(theAl23BetaPlusDecay);

	       //23Al
	        // G4String BetaDecayProcessName = "DPLBetaPlusDecay";
	        // DPLBetaPlus_Decay* theBetaPlusDecay = new DPLBetaPlus_Decay(BetaDecayProcessName);
	        // theBetaPlusDecay->SetParentNucleus(23,13);
	        // theBetaPlusDecay->SetQvalue_Daughter(4056.1*keV);
	        // theBetaPlusDecay->SetBetaDecayData("Al23_BetaIntensity.dat","Branch",12.243*MeV,470.*ms);
	        // pManager->AddRestProcess(theBetaPlusDecay);

	       //20Mg
	       //G4String BetaDecayProcessName = "DPLBetaPlusDecay";
	        //DPLBetaPlus_Decay* theBetaPlusDecay = new DPLBetaPlus_Decay(BetaDecayProcessName);
	        //theBetaPlusDecay->SetParentNucleus(20,12);
	        //theBetaPlusDecay->SetQvalue_Daughter(13892.5*keV);
	        //theBetaPlusDecay->SetBetaDecayData("Mg20_BetaIntensity.dat","Branch",10710*keV,470.*ms);
	        //pManager->AddRestProcess(theBetaPlusDecay);
	       
	       //31Cl
	       // G4String BetaDecayProcessName = "DPLBetaPlusDecay";
	       // DPLBetaPlus_Decay* theBetaPlusDecay = new DPLBetaPlus_Decay(BetaDecayProcessName);
	       // theBetaPlusDecay->SetParentNucleus(31,17);
	       // theBetaPlusDecay->SetBetaDecayData("Cl31_BetaIntensity.dat","Branch",11.9858*MeV,150.*ms);
	       // theBetaPlusDecay->SetQvalue_Daughter(5.39802*MeV);
	       // pManager->AddRestProcess(theBetaPlusDecay,1);
		
		//25Si
	        // G4String BetaDecayProcessName = "DPLBetaPlusDecay";
	        // DPLBetaPlus_Decay* theBetaPlusDecay = new DPLBetaPlus_Decay(BetaDecayProcessName);
	        // theBetaPlusDecay->SetParentNucleus(25,14);
	        // theBetaPlusDecay->SetQvalue_Daughter(4276.6*keV);
		// theBetaPlusDecay->SetBetaDecayData("Si25_BetaIntensity.dat","Branch",12.743*MeV,220.*ms);
		// pManager->AddRestProcess(theBetaPlusDecay);

		//32Ar
	        G4String BetaDecayProcessName = "DPLBetaPlusDecay";
	        DPLBetaPlus_Decay* theBetaPlusDecay = new DPLBetaPlus_Decay(BetaDecayProcessName);
	        theBetaPlusDecay->SetParentNucleus(32,18);
	        theBetaPlusDecay->SetQvalue_Daughter(12680.9*keV);
		theBetaPlusDecay->SetBetaDecayData("Ar32_BetaIntensity.dat","Branch",11.1343*MeV,100.5*ms);
		pManager->AddRestProcess(theBetaPlusDecay);

	       G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
	        radioactiveDecay->SetHLThreshold(-1.*s);
	        radioactiveDecay->SetICM(true);                //Internal Conversion
	        radioactiveDecay->SetARM(false);               //Atomic Rearangement

	        const G4int aMin=29;
	        const G4int aMax=30;
	        const G4int zMin=1;
	        const G4int zMax=15;

	       G4cout<<"Radioactive decay parameters "<<radioactiveDecay->GetNucleusLimits().GetAMin()<<" "<<radioactiveDecay->GetNucleusLimits().GetAMax()
		     <<" "<<radioactiveDecay->GetNucleusLimits().GetZMin()<<" "<<radioactiveDecay->GetNucleusLimits().GetZMax()<<G4endl;



	       G4NucleusLimits limits(aMin, aMax, zMin, zMax);

	       radioactiveDecay->SetNucleusLimits(limits);


	       G4cout<<"Radioactive decay parameters "<<radioactiveDecay->GetNucleusLimits().GetAMin()<<" "<<radioactiveDecay->GetNucleusLimits().GetAMax()
		     <<" "<<radioactiveDecay->GetNucleusLimits().GetZMin()<<" "<<radioactiveDecay->GetNucleusLimits().GetZMax()<<G4endl;

	       //pManager->AddRestProcess(radioactiveDecay);


	       //G4String pDecayProcessName = "1pDecay23Mg";
	       //General_1pDecay* the1pDecay = new General_1pDecay(pDecayProcessName);
	       //the1pDecay->SetParentNucleus(23,12);
	       //pManager->AddDiscreteProcess(the1pDecay);


		// G4String pDecayProcessName = "DPL1pDecay31Cl";
	        // DPL_1pDecay* the1pDecay = new DPL_1pDecay(pDecayProcessName);
		// the1pDecay->SetParentNucleus(31,16);
		// the1pDecay->SetProtonSeparationEnergy(6.1309*MeV);
		//the1pDecay->SetDeltaSeparationEnergy(.5*MeV);
		//the1pDecay->SetDeltaSeparationEnergy(.0*MeV);
		//the1pDecay->SetProtonDecayData("ProtonEnergies_31Cl.dat");
		//pManager->AddDiscreteProcess(the1pDecay,1);



	       // G4String pDecayProcessName = "DPL1pDecay25Si";
	       // DPL_1pDecay* the1pDecay = new DPL_1pDecay(pDecayProcessName);
	       // the1pDecay->SetParentNucleus(25,13);
	       // the1pDecay->SetProtonSeparationEnergy(2.2716*MeV);
	       // the1pDecay->SetDeltaSeparationEnergy(.5*MeV);
	       // the1pDecay->SetDeltaSeparationEnergy(.0*MeV);
	       // the1pDecay->SetProtonDecayData2("ProtonEnergies_25Si.dat");
	       // pManager->AddDiscreteProcess(the1pDecay,1);


	       G4String pDecayProcessName = "DPL1pDecay32Ar";
	       DPL_1pDecay* the1pDecay = new DPL_1pDecay(pDecayProcessName);
	       the1pDecay->SetParentNucleus(32,17);
	       the1pDecay->SetProtonSeparationEnergy(1.5811*MeV);
	       the1pDecay->SetDeltaSeparationEnergy(.5*MeV);
	       the1pDecay->SetDeltaSeparationEnergy(.0*MeV);
	       the1pDecay->SetProtonDecayData2("ProtonEnergies.dat");
	       pManager->AddDiscreteProcess(the1pDecay,1);


	       // G4String pDecayProcessName = "DPL1pDecay";
	       //  DPL_1pDecay* the1pDecay = new DPL_1pDecay(pDecayProcessName);
	       //  the1pDecay->SetParentNucleus(20,11);
	       //  the1pDecay->SetProtonSeparationEnergy(2.190*MeV);
	         //the1pDecay->SetDeltaSeparationEnergy(.5*MeV);
		 //  the1pDecay->SetDeltaSeparationEnergy(.0*MeV);
	       // the1pDecay->SetExcitationEnergy(4.033*MeV);
	       // the1pDecay->SetProtonDecayData("ProtonEnergies_20Mg.dat");
	         //pManager->AddDiscreteProcess(the1pDecay,1);
	       //	G4String aDecayProcessName = "DPLalphaDecay";
	       // DPL_alphaDecay* thealphaDecay = new DPL_alphaDecay(aDecayProcessName);
	       //thealphaDecay->SetParentNucleus(19,10);
	       //thealphaDecay->SetAlphaSeparationEnergy(3.5285*MeV);
		 //pManager->AddDiscreteProcess(thealphaDecay,1);

	  //----------------------------------------------------------------
	  // Em Process options ---- from example TestEM5
	  // Need to set these last to apply to all particles!
	  // Matched to 23Al W1 beta tail data (best result).

	       // Updated 03/24/2010 - Forces 5um step size!

	       G4EmProcessOptions emOptions;

	       //Multiple scattering
	       emOptions.SetMscStepLimitation(fMinimal);
	       emOptions.SetMscLateralDisplacement(true);
	       emOptions.SetLossFluctuations(true);
	       emOptions.SetMscStepLimitation(fUseDistanceToBoundary);//default=fUseSafety

	       //physics tables
	       emOptions.SetMinEnergy(1.*eV);
	       //emOptions.SetMaxEnergy(2.5*GeV);
	       emOptions.SetMaxEnergy(20.0*MeV);
	       emOptions.SetDEDXBinning(15000);
	       emOptions.SetLambdaBinning(15000);

	       //energy loss
	       emOptions.SetLinearLossLimit(0.005);
	       emOptions.SetStepFunction(0.005, 0.005*mm);
	       emOptions.SetRandomStep(true);

	       //ionization
	       //emOptions.SetSubCutoff(false);

}

void ProtonDetectorPhysicsList::SetCuts(){

	  // suppress error messages even in case e/gamma/proton do not exist
	  G4int temp = GetVerboseLevel();
	  SetVerboseLevel(2);
	  //
	  //  Sets the cut on the physics interaction calculations
	  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
	  //   the default cut value for all particle types
	  SetCutsWithDefault();

	  //G4double em_cuts = 0.1*mm;
	  G4double em_cuts = 10.0*cm;   // Avoid Delta-ray production
	  SetCutValue(em_cuts,"gamma");
	  SetCutValue(em_cuts,"e-");
	  SetCutValue(em_cuts,"e+");

	  // Retrieve verbose level
	  SetVerboseLevel(temp);

}


