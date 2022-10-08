#!/bin/bash

# Pass in name and status
function die { echo $1: status $2 ;  exit $2; }

TEST_DIR=src/CalibPPS/TimingCalibration/test

cmsRun  $TEST_DIR/DiamondCalibrationWorker_cfg.py || die "HPTDC PCL failed at worker stage" $?
echo "HPTDC PCL worker succeeded"
cmsRun  $TEST_DIR/DiamondCalibrationHarvester_cfg.py || die "HPTDC PCL failed at harvester stage" $?
echo "HPTDC PCL harvester succeeded"
