import FWCore.ParameterSet.Config as cms

from Validation.CTPPS.simu_config.year_2018_cff import *
profile_2018_preTS1=profile_2018.clone()

#LHCInfo
profile_2018_preTS1.ctppsLHCInfo.xangleBetaStarHistogramObject="2018_preTS1/h2_betaStar_vs_xangle"

# alignment
alignmentFile = "Validation/CTPPS/alignment/2018_preTS1.xml"
profile_2018_preTS1.ctppsRPAlignmentCorrectionsDataXML.MisalignedFiles = [alignmentFile]
profile_2018_preTS1.ctppsRPAlignmentCorrectionsDataXML.RealFiles = [alignmentFile]

# timing not available in this period
profile_2018_preTS1.ctppsDirectSimuData.timeResolutionDiamonds45 = "0.200"
profile_2018_preTS1.ctppsDirectSimuData.timeResolutionDiamonds56 = "0.200"

def UseConstantXangleBetaStar(process, xangle, betaStar):
  process.profile_2018_preTS1.ctppsLHCInfo.xangle = xangle
  process.profile_2018_preTS1.ctppsLHCInfo.betaStar = betaStar

def UseXangleBetaStarHistogram(process,f, obj):
  process.profile_2018_preTS1.ctppsLHCInfo.xangleBetaStarHistogramFile = f
  process.profile_2018_preTS1.ctppsLHCInfo.xangleBetaStarHistogramObject = obj
