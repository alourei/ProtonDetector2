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
// $Id: MargotPrimaryGeneratorAction.cc,v 1.6 2009/06/16 17:47:21 roeder Exp $
// GEANT4 tag $Name: geant4-09-02-patch-01 $
//
// Modified by Brian Roeder LPC Caen on 6/16/2009
// email - roeder@lpccaen.in2p3.fr
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
// 
// 12/13/06 - Added code in "Generate Primaries" to generate a random 
//            energy particle and save the initial energy in a file
// 12/13/06 - Pencil Beam of neutrons at random energy in x-direction
//
// 01/13/06 - Random Energy generator now defined in function randnum() 
//            in header "randnum.h" - Generates a time seeded random number
//            for Random Particle Energy  

#include "MargotPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

// need the below for random theta angle source (from Demon)
 
#include "Randomize.hh"
#include "G4UnitsTable.hh"


// 19 Jan 07 - C++ data io goes to MargotDataRecordTree class 
// through "MargotDataOut" pointer (assigned to "MargotPointer")


MargotPrimaryGeneratorAction::MargotPrimaryGeneratorAction(G4String particleName, G4String SourceType, G4double BeamEng=5.0*MeV) : 
Source(SourceType)
{
  MargotDataOutPG = MargotDataRecordTree::MargotPointer;
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  BeamName = particleName;
  particle_energy = BeamEng;
  dataevent = 0;
}

MargotPrimaryGeneratorAction::~MargotPrimaryGeneratorAction()
{
  delete particleGun;
}


void MargotPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  dataevent++;
  G4int Z,A;

  if(BeamName == "Al23")
    {
      Z=13;
      A=23;
    }
  else if(BeamName == "Cl31")
    {
      Z=17;
      A=31;
    }
  else if(BeamName == "P27")
    {
      Z=15;
      A=27;
    }

  
  if(BeamName == "geantino" || BeamName == "proton")
    {particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(BeamName));}
  else
    {particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));}
   
  //****************************************************
  // Select type of beam
  //****************************************************
  // Source = "test" -> Beam on z-axis at beam energy given in the main()
  // Source = "diffuse" -> Pencil Beam spread out into a diffuse beam spot (11/15/2008)
  // Source = "pencil" -> Pencil beam with +- 5% energy variance
  // Source = "RealBeam" -> 5mm Gaussian in x-y beam spot with 1 deg. divergence
  // Source = "iso"    -> Isotropic Neutron Source
  // Source = "conic"  -> Conic Neutron Source ("semi-isotropic")

  particleGun->SetParticleEnergy(particle_energy);

  G4double momentum_x = 0;   // Default momentum is z-direction for
  G4double momentum_y = 0;   // Pencil beams
  G4double momentum_z = 1;

  G4double theta = 0.; 
  G4double phi= 0.;
  
  const G4double Pi = CLHEP::pi;
  const G4double twopi = 2*Pi;
  G4ThreeVector theStartPos(0.,0.,-25.*cm);
  
  particleGun->SetParticlePosition(theStartPos);

  if(Source == "test")
    {
     particleGun->SetParticlePosition(G4ThreeVector(theStartPos));
     momentum_x = 0.;
     momentum_y = 0.;
     momentum_z = 1.;
     MargotDataOutPG->senddataPG(particle_energy,theta,phi);
    }
else if(Source == "diffuse")
    {
      //Straight Pencil Beam with a diffuse beam spot. 
      //Particle starts at a random position
      //in x-y plane.
      momentum_x = 0.;
      momentum_y = 0.;
      momentum_z = 1.;

      G4double theRadius = 5.0*mm;
      theRadius *= 22.5;          // correction factor 
    
      // Since diffuse beam is forward focused in small theta,
      // Need Pi-small number shaped like cosine to get circular beam

     G4double ran = G4UniformRand();
     theta = acos(1.-0.001*ran);

     ran = G4UniformRand();
     phi = twopi*ran; // Flat in phi

     G4double x_pos = cos(phi)*sin(theta);
     G4double y_pos = sin(phi)*sin(theta);

      theStartPos[0] = theRadius*x_pos;
      theStartPos[1] = theRadius*y_pos;
      theStartPos[2] = -25.*cm;

      particleGun->SetParticlePosition(theStartPos);
    }
  else if(Source == "pencil")
    {
     //position of source of particles
     //Straight Pencil Beam
     momentum_x = 0.;
     momentum_y = 0.;
     momentum_z = 1.;
     particleGun->SetParticlePosition(theStartPos);

     // Now create Beam Energy Spread of +- 0.5 % !

     G4double EnergySpread = 0.01*particle_energy;
     G4double LowLimit = particle_energy-0.5*EnergySpread;
     G4double DeltaEnergy = EnergySpread*G4UniformRand();
     G4double EnergyOut = LowLimit+DeltaEnergy;
     particleGun->SetParticleEnergy(EnergyOut);
     MargotDataOutPG->senddataPG(EnergyOut,theta,phi);
    }
  else if(Source == "RealBeam")
    {
      // 6/29/09 - Calculated that:
      // Coffin slits +/- 1.0cm = +/- 2.5% Energy Spread
      //              +/- 1.5cm = +/- 3.7% Energy Spread
      //              +/- 2.0cm = +/- 5.0% Energy Spread
      // Assuming +/- 1.25% /cm in Energy (-0.64%/cm in Mom.)
      // Spread at coffin slits - BTR (Pg. 18, Sims. logbook). 

      // 5 mm gaussian random starting position in x-y
      // 1 deg. divergence (from starting point)

      // 11/29/09 - Adjusted for deltaP/P = +-0.25% -> deltaE/E = +-0.5%
      // Coffin slits at  ~ +/- 0.4 cm total!
      
      // Set energy spread (as above)
     G4double EnergySpread = 0.01*particle_energy;  // +-0.5% Eng Spread
     G4double LowLimit = particle_energy-0.5*EnergySpread;
     G4double DeltaEnergy = EnergySpread*G4UniformRand();
     G4double EnergyOut = LowLimit+DeltaEnergy;
     particleGun->SetParticleEnergy(EnergyOut);
     
      // Set starting position (random gaussian on x-y)
     G4double FWHM = 5.0*mm;
     G4double conv_factor = 2.35482;   // 2*sqrt(2*log(2)) - converts FWHM to Dev.
     G4double sigma = FWHM/conv_factor;
     theStartPos[0] = (CLHEP::RandGauss::shoot(0.0,sigma))*mm;
     theStartPos[1] = (CLHEP::RandGauss::shoot(0.0,sigma))*mm;
     theStartPos[2] = -9.5*cm;
     particleGun->SetParticlePosition(theStartPos);

     // Set Beam Divergence 
     G4double BeamDiv = 1.*deg;
     G4double ran; 
     ran = G4UniformRand();
 
     // isotropic in cosine for OpenAngle (in degrees)
     theta = acos(1.+(cos(BeamDiv)-1.)*ran);  

     ran = G4UniformRand();
     phi = twopi*ran; // Can do a similar trick with phi.

     momentum_x = cos(phi)*sin(theta);  // source in z-direction
     momentum_y = sin(phi)*sin(theta);
     momentum_z = cos(theta);
     
     MargotDataOutPG->senddataPG(EnergyOut,theta,phi);
    }
   else if(Source == "iso" || Source == "conic")   // Conic or isotropic source
   {
     //particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
     particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 99.9365*mm));

     // theta taken from DemonPrimaryGeneratorAction.cc 
     // OpenAngle gives how many angles are considered in theta
     // OpenAngle = 180.0*deg = isotropic
     // OpenAngle set in this sim to cover face of det MO at 1.75m distance (BTR)
     G4double OpenAngle;

     if(Source == "iso")
       { OpenAngle = 180.0*deg; /* converts automatically to rads! */}
     else if (Source == "conic")
       {
	 G4double z_axis_distance = 1.75*m;
	 G4double rad_det = 8.0*cm;
	 OpenAngle = atan(rad_det/z_axis_distance);
       }
     else
       {
	 G4cout << "********** Error!! No OpenAngle Defined in PrimaryGenerator!" << G4endl;
	 exit(0);
       }
     if(dataevent == 1)
       {G4cout << "OpenAngle = " << OpenAngle/deg << G4endl;}

     G4double ran; 
     ran = G4UniformRand();
 
     // isotropic in cosine for OpenAngle (in degrees)
     theta = acos(1.+(cos(OpenAngle)-1.)*ran);  

     ran = G4UniformRand();
     phi = twopi*ran; // Flat in phi

     momentum_x = cos(phi)*sin(theta);  // source in z-direction
     momentum_y = sin(phi)*sin(theta);
     momentum_z = cos(theta);
   }
 

      // Set particle direction
      G4ThreeVector v(momentum_x,momentum_y,momentum_z);

      // MargotDataOutPG->senddataPG(particle_energy,theta,phi);

      particleGun->SetParticleMomentumDirection(v);
      particleGun->GeneratePrimaryVertex(anEvent);   
}


