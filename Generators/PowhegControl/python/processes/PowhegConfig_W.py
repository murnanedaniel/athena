# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

## @PowhegControl PowhegConfig_W
#  Powheg configuration for W subprocess
#
#  Authors: James Robinson  <james.robinson@cern.ch>
#           Daniel Hayden   <danhayden0@googlemail.com>
#           Stephen Bieniek <stephen.paul.bieniek@cern.ch>

#! /usr/bin/env python
from ..PowhegConfig_base import PowhegConfig_base

## Default Powheg configuration for W generation
#
#  Create a full configurable with all available Powheg options
class PowhegConfig_W(PowhegConfig_base) :

  def __init__( self, runArgs=None, opts=None ) :
    ## Constructor: set process-dependent executable path here
    super(PowhegConfig_W, self).__init__( runArgs, opts )
    self._powheg_executable += '/W/pwhg_main'

    ## Decorate with generic option sets
    self.add_parameter_set( 'CKM' )
    self.add_parameter_set( 'mass window' )
    self.add_parameter_set( 'running scale' )
    self.add_parameter_set( 'second generation quark mass' )
    self.add_parameter_set( 'sin**2 theta W' )
    self.add_parameter_set( 'single vector boson' )
    self.add_parameter_set( 'vector boson decay' )
    self.add_parameter_set( 'W ID' )

    ## Set optimised integration parameters
    self.ncall1  = 120000
    self.ncall2  = 250000
    self.nubound = 20000

    ## Override defaults
    self.mass_low        = 2.5
    self.mass_high       = 2.0 * self.beam_energy
    self.withdamp        = 1
    self.withsubtr       = 1
    self.withnegweights  = 1

    self.populate_default_strings()
