import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RECODQM', Run3)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(0) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# import of standard configurations

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"

process.a1 = cms.EDAnalyzer("StreamThingAnalyzer",
    product_to_get = cms.string('m1')
)

process.wav = cms.EDAnalyzer('WaveformAnalyzer',
     tagDigi = cms.InputTag("totemTimingRawToDigi","TotemTiming"),
     timingCalibrationTag = cms.string("GlobalTag:TotemTimingCalibration"),
)


# raw data source
process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
    'file:/eos/cms/store/t0streamer/Minidaq/A/000/352/678/run352678_ls0001_streamA_StorageManager.dat',
    

),
    inputFileTransitionsEachEvent = cms.untracked.bool(True)
    #firstEvent = cms.untracked.uint64(10123456835)
)

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_relval', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_dataRun3_v2', '')

# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

# local RP reconstruction chain with standard settings
process.load("RecoPPS.Configuration.recoCTPPS_cff")
#process.load('Geometry.VeryForwardGeometry.geometryRPFromDD_2021_cfi')
# CTPPS DQM modules
process.load("DQM.CTPPS.ctppsDQM_cff")
process.ctppsDiamondDQMSource.excludeMultipleHits = cms.bool(True)
process.ctppsDiamondDQMSource.plotOnline = cms.untracked.bool(True)
process.ctppsDiamondDQMSource.plotOffline = cms.untracked.bool(False)
process.path = cms.Path(
    process.a1*
    process.ctppsRawToDigi *
    process.recoCTPPS*
    process.wav
    
    #process.diamondSampicDQMSourceOnline*
    #process.ctppsDQMOnlineHarvest
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("wav.root"),
      closeFileFast = cms.untracked.bool(True)
  )

process.end_path = cms.EndPath(
    process.dqmEnv +
    process.dqmSaver
)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
