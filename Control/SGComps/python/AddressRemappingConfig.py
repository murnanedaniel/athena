""" Helpers for configuring the AddressRemappingSvc

Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
"""

from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator


def AddressRemappingCfg(renameMaps=[]):
    """
    Creates a ComponentAccumulator instance containing the 
    AddressRemappingSvc and other needed services
    """

    from SGComps.SGCompsConf import AddressRemappingSvc, ProxyProviderSvc

    acc = ComponentAccumulator()
    svc = AddressRemappingSvc()
    svc.TypeKeyRenameMaps = renameMaps
    acc.addService(svc)
    acc.addService(ProxyProviderSvc(ProviderNames=["AddressRemappingSvc"])) 

    return acc


def InputRenameCfg(type, from_name, to_name):
    """ Add a new input renaming.

    For example:

      InputRenameCfg ("Foo", "foo", "bar")

    to rename the object of type Foo named `foo' to `bar'.

    May also be used to rename dynamic (but NOT static) auxiliary variables:

      InputRenameCfg ("Foo", "foo.d1", "foo.x1")

    If both are combined, write it like this:

      InputRenameCfg ("Foo", "foo", "bar")
      InputRenameCfg ("Foo", "foo.d1", "bar.x1")
    """

    return AddressRemappingCfg([ '%s#%s->%s' % (type, from_name, to_name) ])


def InputOverwriteCfg(from_type, from_name, to_type, to_name):
    """ Add a new type overwrite mapping. """

    return AddressRemappingCfg([ '%s#%s->%s#%s' % (from_type, from_name,
                                                   to_type, to_name) ])
