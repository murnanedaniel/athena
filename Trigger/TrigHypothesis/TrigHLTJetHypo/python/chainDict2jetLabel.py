# Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
from __future__ import print_function
from __future__ import absolute_import
import re

# substrings that cannot occur in any chainPartName for simple chains.
reject_substr = (
    #    'gsc',
    'ion',
    'dphi',
    'deta',
    'invm',
    '0i1',
    '1i2',
    'dphi',
    r'agg\d',)

reject_substr_res = re.compile(r'%s' % '|'.join(reject_substr))

def make_label(scenario, pattern, template, extra={}):

    r = re.compile(pattern)
    m = r.match(scenario)
    assert m, 'chainDict2jetlabel - pattern %s does not match scenario %s' % (
        pattern, scenario)

    argdict = m.groupdict()
    argdict.update(extra)

    label = template % argdict
    return label


def _select_simple_chainparts(chain_parts):
    """ Reject unsuported chain parts """

    for cp in chain_parts:
        if  reject_substr_res.search(cp['chainPartName']):
            return False

    return True


def _make_simple_label(chain_parts, leg_label):
    """Marshal information deom the selected chainParts to create a
    'simple' label. NOTE: THIS IS A SPECIAL CASE - IT DOES NOT DEPEND
    SOLELY ON THE HYPO SCENARIO.
    Argument leg_label is not used - rather the leg label is fouNd
    from the chain parts.
    """
    
    if not _select_simple_chainparts(chain_parts):
        msg = 'Jet Configuration error: '\
              'chain fails substring selection: not "simple" '

        raise NotImplementedError(msg)

    chainpartind = 0
    label = 'root([]'
    for cp in chain_parts:
        smcstr =  str(cp['smc'])
        jvtstr =  str(cp['jvt'])
        momstr =  str(cp['momCuts'])
        if smcstr == 'nosmc':
            smcstr = ''
        for i in range(int(cp['multiplicity'])):
            label += 'simple(['
            # condition_str = '(%set,%s,%s)' % (str(cp['threshold']),
            #                                  str(cp['etaRange']),
            #                                  smcstr,)
            condition_str = '(%set,%s' % (str(cp['threshold']),
                                              str(cp['etaRange']),)
            if smcstr: # Run 2 chains have "INF" in the SMC substring
                condition_str += ',%s' % smcstr.replace('INF','')
            if jvtstr:
                condition_str += ',%s' % jvtstr
            if momstr:
                if 'SEP' in momstr:
                    print('_cuts_from_momCuts(momstr):')
                    print(_cuts_from_momCuts(momstr))
                    for cut in _cuts_from_momCuts(momstr):
                        condition_str += ',%s' % cut
                else:
                    condition_str += ',%s' % momstr
            condition_str += ', leg{:0>3}'.format(chainpartind)
            if not condition_str.endswith(')'):
                condition_str += ')'
            label += condition_str
            label += '])'
        chainpartind += 1
    label += ')'
    return label


def _args_from_scenario(scenario):
    separator = 'SEP'
    
    args = scenario.split(separator)
    if len(args) > 1:
        return args[1:]
    return ''


def _cuts_from_momCuts(momCuts):
    separator = 'SEP'
    args      = momCuts.split(separator)
    if len(args) > 1:
        return args
    return ''


def _make_fbdjnoshared_label(chain_parts, leg_label):
    """Marshal information from the selected chainParts to create a
    fbdjnoshared (forward-backward and dijet, no jet sharing) label. 
    """

    assert len(chain_parts) == 1
    
    scenario = chain_parts[0]['hypoScenario']
    assert scenario.startswith('fbdjnoshared')

    # example scenario: fbdjnosharedSEP10etSEP20etSEP34massSEP50fbet
    # example label:
    # root([]
    #  simple([(50et, 500neta, leg000)])
    #  simple([(50et, peta500, leg000)])
    #
    #  dijet
    #  (
    #    [(34djmass, 26djdphi)]
    #    simple([(20et, 0eta320, leg000)])
    #    simple([(10et, 0eta320, leg000)])
    #   )
    # )

    pattern = r'^fbdjnosharedSEP'\
        r'(?P<j1etlo>\d*)et(?P<j1ethi>\d*)SEP'\
        r'(?P<j2etlo>\d*)et(?P<j2ethi>\d*)SEP'\
        r'(?P<masslo>\d*)mass(?P<masshi>\d*)SEP'\
        r'(?P<fbetlo>\d*)fbet(?P<fbethi>\d*)$'
        
        
    template = r'root([]'\
        r'simple([(%(fbetlo)set%(fbethi)s, 500neta, %(leg_label)s)])'\
        r'simple([(%(fbetlo)set%(fbethi)s, peta500, %(leg_label)s)])'\
        r'dijet([(%(masslo)sdjmass%(masshi)s, 26djdphi)]'\
        r'simple([(%(j1etlo)set%(j1ethi)s, 0eta320, %(leg_label)s)])'\
        r'simple([(%(j2etlo)set%(j2ethi)s, 0eta320, %(leg_label)s)])))'

    extra = {'leg_label': leg_label}
    
    label = make_label(scenario, pattern, template, extra)
    return label


def  _make_fbdjshared_label(chain_parts, leg_label):
    """example label for a 2-tree forest.
    The fbdjshared contains a dijet and forward backward jets, in separate 
    trees, to allow the fb jets to particoate in the dijet."""

    
    return """
    root([]
    simple([(50et, 500neta, %s)])
    simple([(50et, peta500, %s)])
    )
    root([]
    dijet
    (
    [(34djmass, 26djdphi)]
        simple([(10et, 0eta320, %s)])
        simple([(20et, 0eta320, %s)])
    ))""" % ((leg_label,) * 4)

    
def _make_dijet_label(chain_parts, leg_label):
    """dijet label. supports dijet cuts, and cuts on particpating jets
    Currently supported cuts:
    - dijet mass
    - jet1 et, eta
    - jet2 et, eta

    - default values are used for unspecified cuts.
    The cut set can be extended according to the pattern
    """

    assert len(chain_parts) == 1
    scenario = chain_parts[0]['hypoScenario']
    
    assert scenario.startswith('dijet')

    # example scenario: 'dijetSEP80j1etSEP0j1eta240SEP80j2etSEP0j2eta240SEP700djmass',

    pattern = r'^dijetSEP('\
    r'(?P<j1etlo>\d*)j1et(?P<j1ethi>\d*)SEP'\
    r'(?P<j1etalo>\d*)j1eta(?P<j1etahi>\d*)SEP'\
    r'(?P<j2etlo>\d*)j2et(?P<j2ethi>\d*)SEP'\
    r'(?P<j2etalo>\d*)j2eta(?P<j2etahi>\d*)SEP'\
    r'(?P<djmasslo>\d*)djmass(?P<djmasshi>\d*))$'

    template = 'root([] dijet(' \
        '[(%(djmasslo)sdjmass%(djmasshi)s, 26djdphi)]'\
        'simple([(%(j1etlo)set, %(j1etalo)seta%(j1etahi)s, %(leg_label)s)])'\
        'simple([(%(j2etlo)set, %(j2etalo)seta%(j2etahi)s, %(leg_label)s)])))'

    # example label:
    #    dijet([(700djmass)] simple([(80et, 0eta240, leg002)]) simple([(80et, 0eta240, leg002)])))

    extra = {'leg_label': leg_label}
    label = make_label(scenario, pattern, template, extra)

    return label


def _make_agg_label(chain_parts, leg_label):
    """agg label. cuts on aggregate quantities, and cuts on particpating jets
    Only partway migrated from pure ht to more general agg
    Currently supported cuts:
    - all jets: ht
    - all jets: et
    - all jets: eta

    - default values are used for unspecified cuts.
    The cut set can be extended according to the pattern
    """

    assert len(chain_parts) == 1, '_make_agg_label, no. of chain parts != 1'
    scenario = chain_parts[0]['hypoScenario']
    
    # assert scenario.startswith('agg'), '_make_agg_label(): scenario does not start with agg'

    # the scenario contains the  ht cut, and filter cuts.
    # all cuts thast do no start with 'ht are filter cuts

    pattern = r'^aggSEP(?P<htlo>\d*)ht(?P<hthi>\d*)SEP'\
        r'(?P<etlo>\d*)et(?P<ethi>\d*)SEP'\
        r'(?P<etalo>\d*)eta(?P<etahi>\d*)$'
    
    template = 'root([]agg([(%(htlo)sht, %(leg_label)s)'\
        '(%(etlo)sfltr:et)'\
        '(%(etalo)sfltr:eta%(etahi)s)]))'
    
    extra = {'leg_label': leg_label}
    label = make_label(scenario, pattern, template, extra)

    return label
    

def chainDict2jetLabel(chain_dict):
    """Entry point to this Module. Return a chain label according to the
    value of cp['hypoScenario'], where cp is an element of list/
    chainDict['chainPart']

    Due to historical reasons, the ;logic is as sollows:
    hypoScenatio   Action
    "simple"       examine all chain parts, to greate an EtaEt label.
    other          if len(chainParts) == 1 create correponding  chain label
                   if len(chainParts) > 1 create and of simple and other.
    """

    # suported scenarios. Caution! two keys in the router dict
    # must not share a common initial substring.
    router = {
        'simple': _make_simple_label,
        'agg':   _make_agg_label,
        'dijet': _make_dijet_label,
        'fbdjshared': _make_fbdjshared_label,
        'fbdjnoshared': _make_fbdjnoshared_label,
    }

    # check that no key is the initial susbstring of another key
    # such a case would break the code below.
    keys = sorted(router.keys(), key=len)
    for i in range(1, len(keys)):
        assert not (keys[i].startswith(keys[i-1]))

    # chain_part - scenario association
    cp_sorter = {}
    for k in router: cp_sorter[k] = []

    chain_parts = chain_dict['chainParts']
    for cp in chain_parts:
        if cp['signature'] in ('Jet', 'Bjet'): 
            for k in cp_sorter:
                if cp['hypoScenario'].startswith(k):
                    cp_sorter[k].append(cp)
                    break

    # obtain labels by scenario.
    labels = []

    leg_label = 'leg%03d' % (len(chain_parts) - 1)
    for k, chain_parts in cp_sorter.items():
        if chain_parts:
            labels.append(router[k](chain_parts, leg_label))

    assert labels
    nlabels = len(labels)
    return ''.join(labels)
    if nlabels == 1: return labels[0]
    if nlabels == 2:
        # two labels occur when combining simple and a non-simple scenario
        alabel = """\
all([]
    %s
    %s)""" % (tuple(labels))
        return alabel

    # more than 2 labels is not expected
    assert False

# module tests now in testChainDictMaker.py
