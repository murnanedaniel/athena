# Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

from enum import Enum


class FlagEnum(Enum):
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError(f"Invalid comparison of {self.__class__} with {type(other)}")
        return self is other


class Format(FlagEnum):
    BS = 'BS'
    POOL = 'POOL'


class ProductionStep(FlagEnum):
    # steps should be added when needed
    Default = 'Default'
    Simulation = 'Simulation'
    PileUpPresampling = 'PileUpPresampling'
    Overlay = 'Overlay'
    FastChain = 'FastChain'
    Digitization = 'Digitization'
    Reconstruction = 'Reconstruction'
