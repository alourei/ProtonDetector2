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
// $Id: MargotPrimaryGeneratorAction.cc,v 1.6 2006/11/15 17:47:21 roeder Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
// Modified by Brian Roeder LPC Caen on 11/15/2006
// email - roeder@lpccaen.in2p3.fr
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
//

#ifndef MargotPrimaryGeneratorAction_h
#define MargotPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4LorentzVector.hh"

#include "G4VUserPrimaryGeneratorAction.hh"

#include "MargotDataRecordTree.hh"

class G4ParticleGun;
class G4Event;

class MargotPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
  MargotPrimaryGeneratorAction(){;}
  MargotPrimaryGeneratorAction(G4String pName, G4String SourceType, G4double BeamEng);
  ~MargotPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
   
  private:
    G4ParticleGun* particleGun;
    G4int dataevent;
    G4double particle_energy;
    MargotDataRecordTree* MargotDataOutPG;
    G4int seed;
    G4String Source;
    G4String BeamName;
   
};

#endif


