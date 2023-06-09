# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration

from DataQualityUtils.DQWebDisplayConfig import DQWebDisplayConfig
import os

dqconfig                = DQWebDisplayConfig()
dqconfig.config         = "TEST"

# Use this setup to test your new han configuration before committing it. Set "hcfg_dir" manually if the default doesn't work for you.
# You can change these settings in your local working copy of this file, but please do not commit the change to SVN.
try:
    hcfg_dir = os.environ["TestArea"] + "/DataQuality/DataQualityConfigurations/config"
    print "$TestArea is defined, using it."
except KeyError:
    hcfg_dir = os.getcwd()
    print "$TestArea is not defined, using $PWD instead."
print "Looking for collisions_*.hcfg files in %s" % (hcfg_dir)

dqconfig.hcfg           = "%s/cosmics_run.hcfg"       % (hcfg_dir)
dqconfig.hcfg_min10     = "%s/cosmics_minutes10.hcfg" % (hcfg_dir)
dqconfig.hcfg_min30     = "%s/cosmics_minutes30.hcfg" % (hcfg_dir)

dqconfig.histogramCache = "/afs/cern.ch/user/a/atlasdqm/w1/histogram_web_display_cache"
dqconfig.hanResultsDir  = "/afs/cern.ch/user/a/atlasdqm/dqmdisk/han_results/test"
dqconfig.htmlDir        = "/afs/cern.ch/user/a/atlasdqm/dqmdisk/www/test"
dqconfig.htmlWeb        = "http://atlasdqm.web.cern.ch/atlasdqm/test"
dqconfig.runlist        = "runlist_TEST.xml"
dqconfig.indexFile      = "results_TEST.html"
dqconfig.lockFile       = "DQWebDisplay_TEST.lock"
dqconfig.doHandi        = False
