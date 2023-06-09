#! /usr/bin/env python

# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

# Run a Reco job and test key metadata in the output
#

import glob
import json
import subprocess
import os
import os.path
import sys
import unittest

from PyJobTransforms.trfLogger import msg
from PyJobTransforms.trfReports import pyJobReportToFileDict

sourceFiles = '/afs/cern.ch/work/g/graemes/ddm/pmb/data12_8TeV.00209109.physics_JetTauEtmiss.merge.RAW._lb0186._SFO-1._0001.1'

class Reco_tftest(unittest.TestCase):
    
    def test_runReco_tf(self):
        inputs = glob.glob(sourceFiles)
        self.assertEqual(len(inputs), 1)
        cmd = ['Reco_tf.py', '--inputBSFile']
        cmd.extend(inputs)
        cmd.extend(['--outputESDFile', 'my.ESD.pool.root', '--autoConfiguration', 'everything'])
        cmd.extend(['--outputAODFile', 'my.AOD.pool.root'])
        cmd.extend(['--outputHISTFile', 'my.HIST.root'])
        cmd.extend(['--outputTAGFile', 'my.TAG.pool.root'])
        cmd.extend(['--maxEvents', '10'])
        cmd.extend(['--preExec', 'rec.doTrigger=False'])  # This is temporary while trigger doesn't work in r19
        ## Event counting currently broken for multi-step transforms
        cmd.extend(['--checkEventCount', 'true'])
        msg.info('Will run this transform: {0}'.format(cmd))
        p = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, bufsize = 1)
        while p.poll() is None:
            line = p.stdout.readline()
            sys.stdout.write(line)
        # Hoover up remaining buffered output lines
        for line in p.stdout:
            sys.stdout.write(line)
        self.assertEqual(p.returncode, 0)
        
        # Now load metadata and test a few important values
        with open('jobReport.json') as jr:
            md = json.load(jr)
            self.assertEqual(isinstance(md, dict), True)
            dataDict = pyJobReportToFileDict(md)
            self.assertTrue('ESD' in dataDict.keys())
            self.assertTrue('AOD' in dataDict.keys())
            self.assertTrue('HIST' in dataDict.keys())
            self.assertTrue('TAG' in dataDict.keys())
            self.assertEqual(dataDict['ESD']['subFiles'][0]['nentries'], 10)
            self.assertEqual(dataDict['ESD']['subFiles'][0]['name'], 'my.ESD.pool.root')
            self.assertEqual(dataDict['AOD']['subFiles'][0]['nentries'], 10)
            self.assertEqual(dataDict['AOD']['subFiles'][0]['name'], 'my.AOD.pool.root')
            self.assertEqual(dataDict['HIST']['subFiles'][0]['nentries'], 10)
            self.assertEqual(dataDict['TAG']['subFiles'][0]['nentries'], 10)

if __name__ == '__main__':
    unittest.main()
