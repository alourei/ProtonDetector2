/*
 * ProtonDetectorGeantHit.hh
 *
 *  Created on: Dec 6, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORGEANTHIT_HH_
#define PROTONDETECTORGEANTHIT_HH_

#include <G4VHit.hh>

#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"


class ProtonDetectorGeantHit: public G4VHit {
public:
	ProtonDetectorGeantHit();
	~ProtonDetectorGeantHit();
	ProtonDetectorGeantHit(const ProtonDetectorGeantHit&);

	const ProtonDetectorGeantHit& operator=(const ProtonDetectorGeantHit&);
	G4int operator==(const ProtonDetectorGeantHit&) const;

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	void Draw();
	void Print();
	void PrinttoFile();

	//Getters snd Setters

	void SetTrackID(G4int track){ trackID = track; }
	void SetParentID(G4int track){ parentID = track; }
	void SetEdep(G4double de){ edep = de; }
	void SetParticleCharge(G4double pc){particleCharge=pc;}
	void SetParticleMass(G4double pm){particleMass=pm;}
	void SetParticleID(G4int pi){particleID=pi;}
	void SetPrePos(G4ThreeVector xyz){ prePos = xyz; }
	void SetPostPos(G4ThreeVector xyz){ postPos = xyz; }
	void SetPreToF(G4double Time){ preToF = Time; }
	void SetPostToF(G4double Time){ postToF = Time; }
	void SetDetName(G4String Name){ detName = Name; }
	void SetDetID(G4int id){ detID = id; }
	void SetStepLength(G4double len){ stepLength = len; }
	void SetProcessName(G4String name){ processName = name; }

       void  SetParticleKineticEnergy(G4double energy){particleKineticEnergy=energy;}  

  
	G4int         GetTrackID(){ return trackID; }
	G4int         GetParentID(){ return parentID; }
	G4double      GetEdep(){ return edep; }
	G4double      GetParticleCharge(){return particleCharge;}
	G4double      GetParticleMass(){return particleMass;}
	G4int         GetParticleID(){return particleID;}
	G4ThreeVector GetPrePos(){ return prePos; }
	G4ThreeVector GetPostPos(){ return postPos; }
	G4String      GetDetName(){ return detName; }
	G4int         GetDetID(){ return detID; }
	G4double      GetPreToF(){ return preToF; }
	G4double      GetPostToF(){ return postToF; }
	G4double      GetStepLength(){ return stepLength; }
        G4String      GetProcessName(){return processName;}

        G4double      GetParticleKineticEnergy(){return particleKineticEnergy;}  




private:
	G4int         trackID;
	G4int         parentID;
	G4double      edep;             //energy deposited
	G4double      particleCharge;   // charge of the particle
	G4double      particleKineticEnergy;   // Kinetic energy of the particle
	G4double      particleMass;     // mass of the particle
	G4int         particleID;       // particle ID according to the GDP-coding
	G4ThreeVector postPos;          //position after step
	G4ThreeVector prePos;           //position before step
	G4String      detName;          //detector where energy is deposited
	G4int         detID;            //
	G4double      preToF;           //time before step
	G4double      postToF;          //time after step
        G4double      stepLength;       //length of the step
        G4String      processName;      //process Name

};

typedef G4THitsCollection<ProtonDetectorGeantHit> ProtonDetectorGeantHitsCollection;

extern G4Allocator<ProtonDetectorGeantHit> ProtonDetectorGeantHitAllocator;

inline void* ProtonDetectorGeantHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) ProtonDetectorGeantHitAllocator.MallocSingle();
  return aHit;
}

inline void ProtonDetectorGeantHit::operator delete(void *aHit)
{
  ProtonDetectorGeantHitAllocator.FreeSingle((ProtonDetectorGeantHit*) aHit);
}


#endif /* PROTONDETECTORGEANTHIT_HH_ */
