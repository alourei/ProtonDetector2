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
// $Id: MargotPhysicsList.cc,v 1.6 2007/02/14 17:47:21 roeder Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Modified by Brian Roeder LPC Caen on 02/14/2007
// email - roeder@lpccaen.in2p3.fr
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
//
// 11/16/2006 - Added em processes only, for protons, etc.
//            - Added "ConstructEM()" Process to "ConstructProcess()"
//            - Added "G4Proton::ProtonDefinition()" to add proton 
// 12/12/2006 - Adding hadronic processes for proton and neutron
//
// 02/14/2007 - Copied from SerraPhysicsList
//            - Uses GEANT4 Elastic Scattering for Neutrons (no data)
//            - Uses GEANT4 Inelastic Scattering with JENDL Inelastic data set
//

#include "MargotPhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

// Particle Constructors

#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

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

// Neutron only processes

// Most important part of program is "ConstructProcess()"
// Add EM physics processes in Construct EM.


MargotPhysicsList::MargotPhysicsList()
{;}

MargotPhysicsList::~MargotPhysicsList()
{;}

void MargotPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4Geantino::GeantinoDefinition();
  
  // Particles
  // gamma

  G4Gamma::GammaDefinition();
 
  // leptons
  G4LeptonConstructor pLeptons;
  pLeptons.ConstructParticle();

  //baryons - Construct all baryons (protons and neutrons)
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

void MargotPhysicsList::ConstructProcess()
{
  // Define physics processes

  AddTransportation();
  ConstructEM();
  ConstructNuclear();
}

void MargotPhysicsList::ConstructEM()
{
  // We declare here the particles and interactions for EM processes
  // Baryons and other hadrons defined in ConstructNuclear()  

// Proton EM processes -> Defined below in "Construct Nuclear" Process

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

// Add Nuclear Processes

void MargotPhysicsList::ConstructNuclear()
{
  // Proton nuclear scattering processes

  G4ProcessManager *pManager = 0;
  
   pManager = G4Proton::Proton()->GetProcessManager();
   /*
// Create Pointers for Proton nuclear processes

   G4HadronElasticProcess* theLE_ElasticProcess = new G4HadronElasticProcess();
   G4HadronFissionProcess* theLE_FissionProcess = new G4HadronFissionProcess();
   G4ProtonInelasticProcess* theLE_Proton_InelasticProcess = new G4ProtonInelasticProcess();


   G4LElastic* theLE_Elastic = new G4LElastic();
   G4LFission* theLE_Fission = new G4LFission();
   G4LEProtonInelastic* theLEProton_Inelastic = new G4LEProtonInelastic();

 // Add proton processes
       	theLE_ElasticProcess->RegisterMe(theLE_Elastic);
	theLE_FissionProcess->RegisterMe(theLE_Fission);
	theLE_Proton_InelasticProcess->RegisterMe(theLEProton_Inelastic);

	pManager->AddDiscreteProcess(theLE_ElasticProcess);
	pManager->AddDiscreteProcess(theLE_Proton_InelasticProcess);
	pManager->AddProcess(theLE_FissionProcess);
   */		
	// Add proton EM processes
	pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	pManager->AddProcess(new G4hIonisation,       -1,2,2);

  
// Neutron nuclear scattering processes
//

 pManager = G4Neutron::Neutron()->GetProcessManager();
	
// Create Pointers for Neutron nuclear processes

 			     
// ////// Nuclei ///////////////////////////////////////////////// 
	
// Add Deuteron Physics (EM Stopping Power only)
	pManager = G4Deuteron::Deuteron()->GetProcessManager();
   
	pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	pManager->AddProcess(new G4ionIonisation, -1,2,2);

// Add Triton Physics (EM Stopping Power only)
	pManager = G4Triton::Triton()->GetProcessManager();

	pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	pManager->AddProcess(new G4ionIonisation, -1,2,2);
      
// Add Helium 3 Physics (EM Stopping Power only)
	pManager = G4He3::He3()->GetProcessManager();
        	
	pManager->AddProcess(new G4hMultipleScattering,-1,1,1); 
	pManager->AddProcess(new G4ionIonisation, -1,2,2);
	
		
// Add Alpha Physics (EM Stopping Power only)
	pManager = G4Alpha::Alpha()->GetProcessManager();

	pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	pManager->AddProcess(new G4ionIonisation, -1,2,2);

	// 6-16-2009 - Note for geant4.9.2 and later versions - BTR
	// The ionIonisation GEANT4 process was changed such that it can
	// only be used with Generic Ions as defined in GEANT4 and does not
	// work with "Homemade" particles. To overcome this problem, I call
	// exotic ions as generic ions (not redefined) and then add homemade physics 
	// processes to all ions. I then make the process have infinite mean free path
	// for all ions except the ones I'm interested in. I hope to report 
	// this problem to the GEANT4 people, but I'm not sure they will
	// fix it...

 // Add Generic Ion Physics (EM and Nuclear Elastic, Stopping Power, etc.)
	pManager = G4GenericIon::GenericIon()->GetProcessManager();
	
 // Add EM Processes to Generic Ions
       
	pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
	pManager->AddProcess(new G4ionIonisation, -1,2,2);
         
     //----------------------------------------------------------------
     // Em Process options ---- from example TestEM5

     	// Updated 11/29/2009 - Forces 10um step size!
	G4EmProcessOptions emOptions;

	//Multiple scattering
     emOptions.SetMscStepLimitation(fMinimal);
     emOptions.SetMscLateralDisplacement(true);
     emOptions.SetLossFluctuations(true);
     emOptions.SetMscStepLimitation(fUseDistanceToBoundary);//default=fUseSafety
  
     //physics tables
     emOptions.SetMinEnergy(1.*eV);    
     emOptions.SetMaxEnergy(2.5*GeV);  
     emOptions.SetDEDXBinning(15000);  
     emOptions.SetLambdaBinning(15000);  

     //energy loss
     emOptions.SetLinearLossLimit(1.e-9);
     emOptions.SetStepFunction(0.001, 10.*um);
     emOptions.SetRandomStep(true);
   
     //ionization
     emOptions.SetSubCutoff(false);
    
}

void MargotPhysicsList::SetCuts()
{
  // suppress error messages even in case e/gamma/proton do not exist            
  G4int temp = GetVerboseLevel();                                                
  SetVerboseLevel(2);
  //
  //  Sets the cut on the physics interaction calculations                                                       
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault(); 

  //G4double em_cuts = 1.0*mm;
  G4double em_cuts = 10.*cm; // Use 10cm to avoid generating delta-rays
  SetCutValue(em_cuts,"gamma");
  SetCutValue(em_cuts,"e-");
  SetCutValue(em_cuts,"e+");

  // Retrieve verbose level
  SetVerboseLevel(temp);  
}

