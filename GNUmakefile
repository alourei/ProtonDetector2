# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------\

# The following allows GEANT4 to make ROOT files/histograms

include $(G4INSTALL)/config/architecture.gmk

# ROOT configuration

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)
ROOT_GLIBS  := $(shell root-config --glibs)

CXXFLAGS += $(ROOT_CFLAGS)
CPPFLAGS += $(ROOT_CFLAGS)
LDFLAGS  += $(ROOT_LIBS)
LDFLAGS  += $(ROOT_GLIBS)

# Now Compile Simulation

name := ProtonDetector_Implant
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
