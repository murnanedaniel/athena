# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

test_strings = [
    'simple([(10et, 0eta320)])',
    'simple([(10et, 0eta320)(20et, 0eta320)])',
    'or([] simple([(10et)]) simple([(20et)(40et)]))',
    'and([] simple([(10et)]) simple([(20et)]))',
    'not([] simple([(10et)]))',
    'not([] and([] simple([(40et, 0eta320)]) simple([(100et, 0eta320)])))',
    'or([] not([] simple([(40et, 0eta320)])) not([] simple([(100et, 0eta320)])))',
    'or([] and([] simple([(40et, 0eta320)]) simple([(40et, 0eta320)])) not([] simple([(100et, 0eta320)])))',

    'and([] simple([(50et)(70et)]) dijet([(900djmass, 26djdphi)]))',
    'and([]simple([(50et)(70et)])combgen([]dijet([(900djmass,26djdphi)]) simple([(10et)(20et)])))',
    'and([]simple([(81et)(81et)])combgen([(50et, eta100)]dijet([(26djdphi)])))',
    'simple([(70et,0eta240)(70et,0eta240)(70et,0eta240)(70et,0eta240)(70et,0eta240)])',
    'partgen([(20et,0eta320)]simple([(40et,0eta320)(50et,0eta320)])simple([(35et,0eta240)(55et,0eta240)]))',
    'simple([(10et, neta0)(20et, peta)])', # missing low value for eta - take default
    'simple([(100momwidth200)])', # jet moment condition
    
    # from HLT_j0_vbenfSEP30etSEP34mass35SEP50fbet_L1J20:
    'and([]simple([(30et,500neta)(30et,peta500)])combgen([(10et)]dijet([(34djmass,26djdphi)])simple([(10et)(20et)])))',
    'qjet([(34qjmass)])',
    """partgen([]simple([(neta)(peta)])
                 combgen([]
                         qjet([(qjmass)])
                         partgen([]
                                 combgen([] 
                                         dijet([(djmass)])
                                         simple([(10et)(11et)]))
                                 combgen([] 
                                         dijet([(djmass)])
                                         simple([(12et)(13et)])))))""",
]


if __name__ == '__main__':
    print('There are %d test cases' % len(test_strings))
