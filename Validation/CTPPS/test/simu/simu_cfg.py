import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import *
process = cms.Process('CTPPSTest', Run2_2017)

# load configs
import Validation.CTPPS.simu_config.year_2017_postTS2_cff as config_2017_postTS2
import Validation.CTPPS.simu_config.year_2017_preTS2_cff as config_2017_preTS2
#import Validation.CTPPS.simu_config.year_2018_TS1_TS2_cff as config_2018_TS1_TS2
#import Validation.CTPPS.simu_config.year_2018_postTS2_cff as config_2018_postTS2


#process.load("Validation.CTPPS.simu_config.year_2018_postTS2_cff")
#process.load("Validation.CTPPS.simu_config.year_2018_preTS1_cff")
#process.load("Validation.CTPPS.simu_config.year_2018_TS1_TS2_cff")
process.load("Validation.CTPPS.simu_config.year_2017_postTS2_cff")
process.load("Validation.CTPPS.simu_config.year_2017_preTS2_cff")

config_2017_postTS2.UseConstantXangleBetaStar(process,140,0.3)
config_2017_preTS2.UseConstantXangleBetaStar(process,140,0.3)

#set profiles data
process.profile_2017_postTS2.L_i=1
process.profile_2017_preTS2.L_i=10

process.ctppsCompositeESSource.periods=[process.profile_2017_postTS2,process.profile_2017_preTS2]
process.ctppsCompositeESSource.generateEveryNEvents=1000
process.source.numberEventsInLuminosityBlock=1000


process.Timing=cms.Service("Timing",
	summaryOnly=cms.untracked.bool(True),
	useJobReport=cms.untracked.bool(True)				
)
# minimal logger settings
process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cout'),
  cout = cms.untracked.PSet(
    threshold = cms.untracked.string('WARNING')
  )
)

# number of events
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(100000)
)

# track distribution plotter
process.ctppsTrackDistributionPlotter = cms.EDAnalyzer("CTPPSTrackDistributionPlotter",
  tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),

  rpId_45_F = process.rpIds.rp_45_F,
  rpId_45_N = process.rpIds.rp_45_N,
  rpId_56_N = process.rpIds.rp_56_N,
  rpId_56_F = process.rpIds.rp_56_F,

  outputFile = cms.string("simu_tracks.root")
)


# reconstruction plotter
process.ctppsProtonReconstructionPlotter = cms.EDAnalyzer("CTPPSProtonReconstructionPlotter",
  tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
  tagRecoProtonsSingleRP = cms.InputTag("ctppsProtons", "singleRP"),
  tagRecoProtonsMultiRP = cms.InputTag("ctppsProtons", "multiRP"),

  rpId_45_F = process.rpIds.rp_45_F,
  rpId_45_N = process.rpIds.rp_45_N,
  rpId_56_N = process.rpIds.rp_56_N,
  rpId_56_F = process.rpIds.rp_56_F,

  association_cuts_45 = process.ctppsProtons.association_cuts_45,
  association_cuts_56 = process.ctppsProtons.association_cuts_56,

  outputFile = cms.string("simu_protons.root")
)


# processing path
process.p = cms.Path(
  process.generator
  * process.beamDivergenceVtxGenerator
  * process.ctppsDirectProtonSimulation

  * process.reco_local
  * process.ctppsProtons

  * process.ctppsTrackDistributionPlotter
  * process.ctppsProtonReconstructionPlotter
)
