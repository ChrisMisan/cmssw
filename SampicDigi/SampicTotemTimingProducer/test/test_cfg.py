import FWCore.ParameterSet.Config as cms

process = cms.Process("SampicDigi")

process.load('RecoPPS.Local.totemTimingLocalReconstruction_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("EmptySource")

process.totemTimingRawToDigi = cms.EDProducer('SampicTotemTimingProducer',
	sampicFilesVec=cms.vstring("/eos/user/k/kmisan/sampic/CMSSW_11_1_8/src/SampicDigi/SampicTotemTimingProducer/test/Ntuple_runsampic_159_runtelescope_636.root"),
	idsMapping = cms.VPSet(
		cms.PSet(detId = cms.vuint32(2054160384,2054553600), treeChId = cms.uint32(8)),
		cms.PSet(detId = cms.vuint32(2054164480,2054557696), treeChId = cms.uint32(9)),
		cms.PSet(detId = cms.vuint32(2054168576,2054561792), treeChId = cms.uint32(10)),
		cms.PSet(detId = cms.vuint32(2054172672,2054565888), treeChId = cms.uint32(11)),
		cms.PSet(detId = cms.vuint32(2054176768,2054569984), treeChId = cms.uint32(12)),
		cms.PSet(detId = cms.vuint32(2054180864,2054574080), treeChId = cms.uint32(13)),
		cms.PSet(detId = cms.vuint32(2054184960,2054578176), treeChId = cms.uint32(14)),
		cms.PSet(detId = cms.vuint32(2054189056,2054582272), treeChId = cms.uint32(15)),
		cms.PSet(detId = cms.vuint32(2054193152,2054586368), treeChId = cms.uint32(16)),
		cms.PSet(detId = cms.vuint32(2054197248,2054590464), treeChId = cms.uint32(17)),
		cms.PSet(detId = cms.vuint32(2054201344,2054594560), treeChId = cms.uint32(18)),
		cms.PSet(detId = cms.vuint32(2054205440,2054598656), treeChId = cms.uint32(19)),

		cms.PSet(detId = cms.vuint32(2054291456,2054422528), treeChId = cms.uint32(20)),
		cms.PSet(detId = cms.vuint32(2054295552,2054426624), treeChId = cms.uint32(21)),
		cms.PSet(detId = cms.vuint32(2054299648,2054430720), treeChId = cms.uint32(22)),
		cms.PSet(detId = cms.vuint32(2054303744,2054434816), treeChId = cms.uint32(23)),
		cms.PSet(detId = cms.vuint32(2054307840,2054438912), treeChId = cms.uint32(24)),
		cms.PSet(detId = cms.vuint32(2054311936,2054443008), treeChId = cms.uint32(25)),
		cms.PSet(detId = cms.vuint32(2054316032,2054447104), treeChId = cms.uint32(26)),
		cms.PSet(detId = cms.vuint32(2054320128,2054451200), treeChId = cms.uint32(27)),
		cms.PSet(detId = cms.vuint32(2054324224,2054455296), treeChId = cms.uint32(28)),
		cms.PSet(detId = cms.vuint32(2054328320,2054459392), treeChId = cms.uint32(29)),
		cms.PSet(detId = cms.vuint32(2054332416,2054463488), treeChId = cms.uint32(30)),
		cms.PSet(detId = cms.vuint32(2054336512,2054467584), treeChId = cms.uint32(31))

	)
)
process.load('Geometry.VeryForwardGeometry.geometryRPFromDD_2018_cfi')
process.totemTimingRecHits.timingCalibrationTag= cms.string('ppsTimingCalibrationESSource:TotemTimingCalibration')
process.ppsTimingCalibrationESSource = cms.ESSource('PPSTimingCalibrationESSource',
  calibrationFile = cms.FileInPath(''),
  subDetector = cms.uint32(1),
  appendToDataLabel = cms.string('TotemTimingCalibration')
)

process.ppsTimingCalibrationESSource.calibrationFile = cms.FileInPath('RecoPPS/Local/data/adjusted_offsets_withcal.cal.json')

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cout'),
  cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO')
  )
)
process.content = cms.EDAnalyzer("EventContentAnalyzer")


from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
process.totemTimingDQMSource = DQMEDAnalyzer('TotemTimingDQMSource',
    tagDigi = cms.InputTag("totemTimingRawToDigi", "TotemTiming"),
    tagFEDInfo = cms.InputTag("totemTimingRawToDigi", "TotemTiming"),
    tagRecHits = cms.InputTag("totemTimingRecHits"),
    tagTracks = cms.InputTag("totemTimingLocalTracks"),
    tagLocalTrack = cms.InputTag("totemRPLocalTrackFitter"),

    minimumStripAngleForTomography = cms.double(0),
    maximumStripAngleForTomography = cms.double(1),
    samplesForNoise = cms.untracked.uint32(6),

    verbosity = cms.untracked.uint32(10),
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"

process.p = cms.Path(process.totemTimingRawToDigi*
process.totemTimingLocalReconstruction*
process.totemTimingDQMSource)

#process.e = cms.EndPath(process.out)
process.end_path = cms.EndPath(
    process.dqmEnv +
    process.dqmSaver
)

process.schedule = cms.Schedule(
    process.p,
    process.end_path
)


