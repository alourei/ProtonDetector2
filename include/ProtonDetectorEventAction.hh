/*
 * ProtonDetectorEventAction.hh
 *
 *  Created on: Dec 11, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTOREVENTACTION_HH_
#define PROTONDETECTOREVENTACTION_HH_

#include <G4UserEventAction.hh>
#include "globals.hh"
#include "G4ThreeVector.hh"


class ProtonDetectorEventAction: public G4UserEventAction {
private:
	  G4String  drawFlag;
	  G4int     printModulo;

public:

	ProtonDetectorEventAction();
	~ProtonDetectorEventAction();
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

	void SetDrawFlag   (G4String val)  {drawFlag = val;};
	void SetPrintModulo(G4int    val)  {printModulo = val;};

};

#endif /* PROTONDETECTOREVENTACTION_HH_ */
