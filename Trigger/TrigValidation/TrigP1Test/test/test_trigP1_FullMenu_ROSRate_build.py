#!/usr/bin/env python

# art-description: Same as full_menu test from TrigUpgradeTest, but with athenaHLT, and adding ROS simulation
# art-type: build                                                                  
# art-include: master/Athena                                                       

from TrigValTools.TrigValSteering import Test, ExecStep, CheckSteps, Step

ex = ExecStep.ExecStep()
ex.type = 'athenaHLT'
ex.job_options = 'TrigUpgradeTest/full_menu.py'
ex.input = 'data'
ex.args = '-c "doWriteESD=False" --ros2rob /cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/TrigP1Test/ATLASros2rob2018-r22format.py'
ex.perfmon = False # perfmon with athenaHLT doesn't work at the moment

ros2json = CheckSteps.InputDependentStep("RosRateToJson")
ros2json.executable = 'ros-hitstats-to-json.py'
ros2json.input_file = 'ros_hitstats_reject.txt'
ros2json.output_stream = Step.Step.OutputStream.STDOUT_ONLY

test = Test.Test()
test.art_type = 'build'
test.exec_steps = [ex]
test.check_steps = CheckSteps.default_check_steps(test)
test.check_steps.append(ros2json)

import sys
sys.exit(test.run())
