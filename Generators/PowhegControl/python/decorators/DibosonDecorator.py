# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

## @PowhegControl DibosonDecorator
#  Powheg runcard decorator for diboson settings
#
#  Authors: James Robinson  <james.robinson@cern.ch>

#! /usr/bin/env python

class DibosonDecorator(object) :

  ## Define decorator name string
  name = 'diboson'

  def __init__( self, decorated ) :
    ## Attach decorations to Powheg configurable
    decorated.run_card_decorators.append( self )
    self.decorated = decorated

    self.decorated.allowed_decay_modes = []

    self.decorated.add_phantom(   'decay_mode', None, default='{0}', desc='Diboson decay mode' )
    self.decorated.add_parameter( 'dronly', 0,        default='{0}', desc='(0:disabled; 1:enabled) include only double resonant diagrams' )


  def finalise( self ) :
    # __decay_mode = self.decorated.pop('decay_mode')
    if self.decorated.decay_mode not in self.decorated.allowed_decay_modes :
      self.decorated.logger.warning( 'Decay mode {0} not recognised!'.format( self.decorated.decay_mode) )
    ## Add entry for each decay mode
    for decay_mode in self.decorated.allowed_decay_modes :
      self.decorated.fix_parameter( decay_mode, [-1,1][decay_mode==self.decorated.decay_mode], default='{0}', desc='(default user-configured) Diboson decay mode, 1:enabled; -1:disabled' )
