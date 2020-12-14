/*
 * ProtonDetectorRunAction.hh
 *
 *  Created on: Dec 11, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORRUNACTION_HH_
#define PROTONDETECTORRUNACTION_HH_

#include "G4UserRunAction.hh"

#include <G4RunManager.hh>
#include "globals.hh"

class G4run;

class ProtonDetectorRunAction: public G4UserRunAction {
public:
	ProtonDetectorRunAction();
	~ProtonDetectorRunAction();

	void BeginOfRunAction(const G4Run*);
	void EndOfRunAction(const G4Run*);


};

#endif /* PROTONDETECTORRUNACTION_HH_ */
