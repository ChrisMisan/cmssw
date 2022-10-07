#!/bin/bash

# Pass in name and status
function die { echo $1: status $2 ;  exit $2; }



cmsRun DiamondCalibrationWorker_cfg.py || die "HPTDC PCL failed at worker stage" $?
echo "HPTDC PCL worker succeeded"
cmsRun DiamondCalibrationHarvester_cfg.py || die "HPTDC PCL failed at harvester stage" $?
echo "HPTDC PCL harvester succeeded"
