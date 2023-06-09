# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

import os.path
from RTTSException   import RTTCodingError

class RTTdict:
    __module__ = os.path.splitext(os.path.basename(__file__))[0]

    def __init__(self, rtt={}, user={}):
        self.rtt  = rtt
        self.user = user
        
    def __str__(self):
        out = ''
        for k,v in self.items():
            out += '%s => %s\n' % (k,v)
        return out
    
    def whoHasKey(self, k):
        if k in self.rtt.keys(): return self.rtt
        if k in self.user.keys(): return self.user
        return None
    
    def __getitem__(self, k):
        d = self.whoHasKey(k)
        if d is not None: return d[k] 
        m = 'KeyError: Unknown key %s' %k
        # raise RTTCodingError(m)
        raise KeyError(m)
            
    def __setitem__(self, k, v):
        d = self.whoHasKey(k)
        if d is not None: d[k] = v
        else:
            self.user[k] = v

    def __delitem__(self, k):
        if k in self.rtt.keys(): del self.rtt[k]
        if k in self.user.keys(): del self.user[k]
        m = 'KeyError: Unknown key %s' %k
        raise RTTCodingError(m)
    
    def get(self, k, default=None):
        try:
            return self.__getitem__(k)
        except KeyError:
            return default

    def __contains__(self, k):
        return k in self.rtt or k in self.user

    def has_key(self, k):
        return self.__contains__(k)

    def __len__(self):
        return len(self.rtt)+len(self.user)

    def keys(self):
        k = self.rtt.keys()
        k.extend(self.user.keys())
        return k

    def values(self):
        v = self.rtt.values()
        v.extend(self.user.values())
        return v

    def items(self):
        i = self.rtt.items()
        i.extend(self.user.items())
        return [(k,v) for k,v in i]
        
    def update(self, d):
        for k,v in d.items():
            self.__setitem__(k,v)

    def setdefault(self, k, v=None):
        self.__setitem__(k,v)
        return v

        
