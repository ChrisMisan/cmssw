import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

ALCARECOPPSCalMaxTracksFilter = hlt.hltHighLevel.clone(
  eventSetupPathsKey = 'PPSCalMaxTracks',
)

from EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff import *
from EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff import ctppsDiamondRawToDigi as _ctppsDiamondRawToDigi
from EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff import totemTimingRawToDigi as _totemTimingRawToDigi
from EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff import ctppsPixelDigis as _ctppsPixelDigis

ctppsDiamondRawToDigiAlCaRecoProducer  = _ctppsDiamondRawToDigi.clone(rawDataTag = cms.InputTag('hltPPSCalibrationRaw'))
totemTimingRawToDigiAlCaRecoProducer = _totemTimingRawToDigi.clone(rawDataTag = cms.InputTag('hltPPSCalibrationRaw') )
#totemTimingRawToDigiAlCaRecoProducer.RawToDigi.verbosity = 1
#totemTimingRawToDigiAlCaRecoProducer.RawToDigi.printErrorSummary = 1
ctppsPixelDigisAlCaRecoProducer = _ctppsPixelDigis.clone(inputLabel = cms.InputTag('hltPPSCalibrationRaw'))

ctppsRawToDigiTaskAlCaRecoProducer = cms.Task(
  ctppsDiamondRawToDigiAlCaRecoProducer,
  totemTimingRawToDigiAlCaRecoProducer,
  ctppsPixelDigisAlCaRecoProducer
)

ALCARECOPPSCalMaxTracksRaw2Digi = cms.Sequence(ctppsRawToDigiTaskAlCaRecoProducer)

from RecoPPS.Configuration.recoCTPPS_cff import *


###################################################################################################

from RecoPPS.Local.ctppsPixelClusters_cfi import ctppsPixelClusters as _ctppsPixelClusters
from RecoPPS.Local.ctppsPixelRecHits_cfi import ctppsPixelRecHits as _ctppsPixelRecHits
from RecoPPS.Local.ctppsPixelLocalTracks_cfi import ctppsPixelLocalTracks as _ctppsPixelLocalTracks

ctppsPixelClustersAlCaRecoProducer = _ctppsPixelClusters.clone(tag = 'ctppsPixelDigisAlCaRecoProducer')
ctppsPixelRecHitsAlCaRecoProducer = _ctppsPixelRecHits.clone(RPixClusterTag = 'ctppsPixelClustersAlCaRecoProducer')
ctppsPixelLocalTracksAlCaRecoProducer = _ctppsPixelLocalTracks.clone(tag = 'ctppsPixelRecHitsAlCaRecoProducer')

ctppsPixelLocalReconstructionTaskAlCaRecoProducer = cms.Task(
    ctppsPixelClustersAlCaRecoProducer,
    ctppsPixelRecHitsAlCaRecoProducer,
    ctppsPixelLocalTracksAlCaRecoProducer
)

###################################################################################################

from RecoPPS.Local.ctppsDiamondRecHits_cfi import ctppsDiamondRecHits as _ctppsDiamondRecHits
from RecoPPS.Local.ctppsDiamondLocalTracks_cfi import ctppsDiamondLocalTracks as _ctppsDiamondLocalTracks
ctppsDiamondRecHitsAlCaRecoProducer = _ctppsDiamondRecHits.clone(digiTag = cms.InputTag("ctppsDiamondRawToDigiAlCaRecoProducer","TimingDiamond"))
ctppsDiamondLocalTracksAlCaRecoProducer = _ctppsDiamondLocalTracks.clone(recHitsTag = 'ctppsDiamondRecHitsAlCaRecoProducer')

ctppsDiamondLocalReconstructionTaskAlCaRecoProducer = cms.Task(
    ctppsDiamondRecHitsAlCaRecoProducer,
    ctppsDiamondLocalTracksAlCaRecoProducer
)

###################################################################################################

from RecoPPS.Local.totemTimingRecHits_cfi import totemTimingRecHits as _totemTimingRecHits
from RecoPPS.Local.diamondSampicLocalTracks_cfi import diamondSampicLocalTracks as _diamondSampicLocalTracks
#from RecoPPS.Local.totemTimingLocalTracks_cfi import totemTimingLocalTracks as _totemTimingLocalTracks

totemTimingRecHitsAlCaRecoProducer = _totemTimingRecHits.clone(digiTag = 'totemTimingRawToDigiAlCaRecoProducer')
diamondSampicLocalTracksAlCaRecoProducer = _diamondSampicLocalTracks.clone(recHitsTag = 'totemTimingRecHitsAlCaRecoProducer')
#totemTimingLocalTracksAlCaRecoProducer = _totemTimingLocalTracks.clone( recHitsTag= 'totemTimingRecHits')

diamondSampicLocalReconstructionTaskAlCaRecoProducer = cms.Task(
		totemTimingRecHitsAlCaRecoProducer,
    diamondSampicLocalTracksAlCaRecoProducer
    )

###################################################################################################

from RecoPPS.Local.ctppsLocalTrackLiteProducer_cff import ctppsLocalTrackLiteProducer as _ctppsLocalTrackLiteProducer

ctppsLocalTrackLiteProducerAlCaRecoProducer = _ctppsLocalTrackLiteProducer.clone(
    includeStrips = False,
    includeDiamonds = False,  # TODO switch to true
    includePixels = True,
    tagDiamondTrack = 'ctppsDiamondLocalTracksAlCaRecoProducer',
    tagPixelTrack = 'ctppsPixelLocalTracksAlCaRecoProducer'
)

###################################################################################################

from RecoPPS.ProtonReconstruction.ctppsProtons_cff import ctppsProtons as _ctppsProtons

ctppsProtonsAlCaRecoProducer = _ctppsProtons.clone(tagLocalTrackLite = 'ctppsLocalTrackLiteProducerAlCaRecoProducer')

###################################################################################################


recoPPSTaskAlCaRecoProducer = cms.Task(
    #ctppsDiamondLocalReconstructionTaskAlCaRecoProducer ,
    #diamondSampicLocalReconstructionTaskAlCaRecoProducer ,
    ctppsPixelLocalReconstructionTaskAlCaRecoProducer ,
    ctppsLocalTrackLiteProducerAlCaRecoProducer ,
    ctppsProtonsAlCaRecoProducer
)

recoPPSSequenceAlCaRecoProducer = cms.Sequence(recoPPSTaskAlCaRecoProducer)

seqALCARECOPPSCalMaxTracksReco = cms.Sequence( ALCARECOPPSCalMaxTracksFilter  + ALCARECOPPSCalMaxTracksRaw2Digi  + recoPPSSequenceAlCaRecoProducer)
#seqALCARECOPPSCalMaxTracksReco = cms.Sequence( ALCARECOPPSCalMaxTracksFilter  + ALCARECOPPSCalMaxTracksRaw2Digi)
