/*
 * ProtonDetectorPhysicsList.hh
 *
 *  Created on: Jan 30, 2014
 *      Author: perezlou
 */

#ifndef PROTONDETECTORPHYSICSLIST_HH_
#define PROTONDETECTORPHYSICSLIST_HH_

#include <G4VUserPhysicsList.hh>
#include <globals.hh>

// Included Physics Lists

// EM processes
// Hadrons
#include <G4eMultipleScattering.hh>
#include "G4hIonisation.hh"

#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"

#include "G4StepLimiter.hh"

//Low Energy EM leptons and Boson Processes (following ExNO2)
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

// Strong force processes ("Low" Energy only)
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"


class ProtonDetectorPhysicsList: public G4VUserPhysicsList {
public:
	ProtonDetectorPhysicsList();
	~ProtonDetectorPhysicsList();

protected:
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

    // these methods Construct physics processes and register them
    // construct proton em scattering and energy-loss processes
    // construct Low Energy Scattering processes for Neutrons
      void ConstructEM();
      void ConstructNuclear();

    //EM process pointers
      // Gamma physics
      G4PhotoElectricEffect thePhotoEffect;
      G4ComptonScattering theComptonEffect;
      G4GammaConversion thePairProduction;

      // Electron physics
      G4eMultipleScattering theElectronMultipleScattering;
      G4eIonisation theElectronIonisation;
      G4eBremsstrahlung theElectronBremsStrahlung;

      //Positron physics
      G4eMultipleScattering thePositronMultipleScattering;
      G4eIonisation thePositronIonisation;
      G4eBremsstrahlung thePositronBremsStrahlung;
      G4eplusAnnihilation theAnnihilation;


};

#endif /* PROTONDETECTORPHYSICSLIST_HH_ */
