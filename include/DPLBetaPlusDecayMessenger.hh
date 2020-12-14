#ifndef DPLBetaPlusDecayMessenger_h
#define DPLBetaPlusDecayMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DPLBetaPlusDecay;


class  DPLBetaPlusDecayMessenger: public G4UImessenger {
private:

  DPLBetaPlusDecay* theBetaDecay;
  G4UIdirectory*               DecayDir;

#endif
