/*
 * ProtonDetectorSD.hh
 *
 *  Created on: Dec 6, 2013
 *      Author: perezlou
 */

#ifndef VETOSD_HH_
#define VETOSD_HH_

#include <G4VSensitiveDetector.hh>
#include <ProtonDetectorGeantHit.hh>

class G4Step;
class G4HCofThisEvent;

class VetoSD: public G4VSensitiveDetector {
public:
	VetoSD(G4String name);
	~VetoSD();

	void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	void EndOfEvent(G4HCofThisEvent*);


private:

	ProtonDetectorGeantHitsCollection* hitsCollection; //Geant step-like hits collect.

};

#endif /* VETOSD_HH_ */
