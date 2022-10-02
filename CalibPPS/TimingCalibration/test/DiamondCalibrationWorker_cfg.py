import FWCore.ParameterSet.Config as cms

process = cms.Process("worker")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('RecoPPS.Local.totemTimingLocalReconstruction_cff')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.a1 = cms.EDAnalyzer("StreamThingAnalyzer",
    product_to_get = cms.string('m1')
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_frozen_v4', '')
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")
process.load("RecoPPS.Configuration.recoCTPPS_cff")

#process.load('CondCore.CondDB.CondDB_cfi')
#process.CondDB.connect = 'sqlite_file:ppsDiamondTiming_calibration.sqlite' # SQLite input
#process.PoolDBESSource = cms.ESSource('PoolDBESSource',
#        process.CondDB,
#        DumpStats = cms.untracked.bool(True),
#        toGet = cms.VPSet(
#            cms.PSet(
#                record = cms.string('PPSTimingCalibrationRcd'),
#                tag = cms.string('DiamondTimingCalibration')
#        )
#    )
#)
    
# raw data source
#process.source = cms.Source('PoolSource',
#    fileNames = cms.untracked.vstring(
#    'file:/eos/cms/store/group/dpg_ctpps/comm_ctpps/AlcaReco/354332/00000/64f0826f-49e0-4876-9a4a-f28d5e97170d.root'
#    
#
#),
#)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/eos/cms/store/data/Run2022B/AlCaPPS/RAW/v1/000/355/207/00000/c23440f4-49c0-44aa-b8f6-f40598fb4705.root',
    

),
    #inputFileTransitionsEachEvent = cms.untracked.bool(True)
    #firstEvent = cms.untracked.uint64(10123456835)
)
################
#geometry
################
#process.load('Geometry.VeryForwardGeometry.geometryRPFromDD_2021_cfi')

process.load("CalibPPS.TimingCalibration.ppsTimingCalibrationPCLWorker_cfi")
process.DQMStore = cms.Service("DQMStore")

process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
    fileName = cms.untracked.string("worker_output.root")
)

process.load("CalibPPS.TimingCalibration.PPSDiamondSampicTimingCalibrationPCLWorker_cfi")

process.load("CalibPPS.TimingCalibration.ALCARECOPromptCalibProdPPSTimingCalib_cff")

process.ctppsPixelDigis.inputLabel = cms.InputTag("hltPPSCalibrationRaw")
process.ctppsDiamondRawToDigi.rawDataTag = cms.InputTag("hltPPSCalibrationRaw")
process.totemRPRawToDigi.rawDataTag = cms.InputTag("hltPPSCalibrationRaw")
process.totemTimingRawToDigi.rawDataTag = cms.InputTag("hltPPSCalibrationRaw")

process.path = cms.Path(
    #process.a1* 
    process.ctppsRawToDigi *
    process.recoCTPPS *
    process.ppsTimingCalibrationPCLWorker
    #process.PPSDiamondSampicTimingCalibrationPCLWorker
)

process.end_path = cms.EndPath(
    process.dqmOutput
)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
