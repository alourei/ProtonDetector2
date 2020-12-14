//==========================================================================
// PerfectDetHit.hh
// Based on DemonScintHit.hh and ExN04TrackerHit.hh
//
// Written/Modified by: Brian Roeder, LPC Caen 02/14/07
//                      email - roeder@lpccaen.in2p3.fr
//
// Purpose: Defines hit collection and values to be stored by PerfectSD.cc
//
//==========================================================================
//
// - See UM Hits Presentation and JLab Scoring 2 Talk for more info.


#ifndef PerfectDetHit_h
#define PerfectDetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


class PerfectDetHit : public G4VHit
{
  public:

  PerfectDetHit();        // Constructor
  ~PerfectDetHit();       // Destructor

  /* //From ExN04TrackerHit.cc - Not needed... 
  PerfectDetHit(const PerfectDetHit &right);
  const PerfectDetHit& operator=(const PerfectDetHit &right);
  G4int operator==(const PerfectDetHit &right) const;
  */

  inline void * operator new(size_t);
  inline void operator delete(void *aHit);
  

  void Draw();       // needed to Draw Hit in Visualization
  void Print();      // Print out Event (if not invoked in EndofEventAction)

  private:
  // These are the Data Storage objects of the class.
  // Set/Stored and Accessed through Public Methods below.

  G4double KinEng;           // Records Neutron Energy in Hit
  G4ThreeVector pos;       // Records Hit Position in Detector
  G4double TOF;            // Records Time of Flight (TOF) of particle hit

  G4String ParticleName;   // Records particle ID in reactions 
  G4String CreatorProcess; // Records the "Process" creating the Track in Hit

  G4double thePDGCharge;   // Records the PDGCharge of the particle
  G4double thePDGMass;

  G4int DetCopyNumber;
  G4String DetName;

 public:

  // Set Data Methods
      void SetKEDep(G4double e)
      { KinEng = e; }
      void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      void SetPDGCharge(G4double theCharge)
      { thePDGCharge = theCharge; }
      void SetPDGMass(G4double theMass)
      { thePDGMass = theMass; }
      void SetTOF(G4double hittime)
      { TOF = hittime; }
      void SetParticleName(G4String theParticleName)
      {ParticleName = theParticleName;}
      void SetParticleProcess(G4String theProcess)
      {CreatorProcess = theProcess;}
      void SetDetCopyNumber(G4int theCopy)
      {DetCopyNumber = theCopy;}
      void SetDetName(G4String theName)
      {DetName = theName;}
     
  // Get Data Methods (For EndofEventAction and DataRecordTree)

      G4double GetKEDep()
      { return KinEng; }
      G4ThreeVector GetPos()
      { return pos; }
      G4double GetPDGCharge()
      {return thePDGCharge;}
      G4double GetPDGMass()
      {return thePDGMass;}
      G4double GetTOF()
      { return TOF; }
      G4String GetParticleName()
      { return ParticleName; }
      G4String GetParticleProcess()
      { return CreatorProcess; }
      G4int GetDetCopyNumber()
      { return DetCopyNumber; }
      G4String GetDetName()
      { return DetName; }
};

typedef G4THitsCollection<PerfectDetHit> PerfectDetHitsCollection;

// These functions setup the memory functions for the Hit

extern G4Allocator<PerfectDetHit> PerfectDetHitAllocator;

// Create the Hit

inline void* PerfectDetHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *)PerfectDetHitAllocator.MallocSingle();
  return aHit;
}

// Delete the Hit (at the end of the event, for example)

inline void PerfectDetHit::operator delete(void *aHit)
{
  PerfectDetHitAllocator.FreeSingle((PerfectDetHit*) aHit);
}

#endif
