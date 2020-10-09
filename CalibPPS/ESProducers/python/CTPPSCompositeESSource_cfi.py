import FWCore.ParameterSet.Config as cms

ctppsCompositeESSource = cms.ESSource("CTPPSCompositeESSource",
lhcInfoLabel=cms.string(""),
opticsLabel=cms.string(""),
generateEveryNEvents = cms.uint32(1),
seed=cms.uint32(1),
periods=cms.VPSet(),
verbosity=cms.untracked.uint32(1)
)
