import FWCore.ParameterSet.Config as cms

from Validation.CTPPS.simu_config.year_2016_cff import *

profile_2016_preTS2=cms.PSet(
  L_i=cms.double(1),
  #LHCInfo
  ctppsLHCInfo = cms.PSet(
	xangle=cms.double(-1),
	betaStar=cms.double(-1),
  	beamEnergy = cms.double(6500),  # GeV
  	xangleBetaStarHistogramFile=cms.string(default_xangle_beta_star_file),
  	xangleBetaStarHistogramObject=cms.string("2016_preTS2/h2_betaStar_vs_xangle")
  ),

  #Optics
  ctppsOpticalFunctions = cms.PSet(

  	opticalFunctions = cms.VPSet(
    		cms.PSet( xangle = cms.double(185), fileName = cms.FileInPath("CalibPPS/ESProducers/data/optical_functions/2016_preTS2/version2/185urad.root") )
  	),

  	scoringPlanes = cms.VPSet(
	      # z in cm
	      cms.PSet( rpId = cms.uint32(0x76100000), dirName = cms.string("XRPH_C6L5_B2"), z = cms.double(-20382.6) ),  # RP 002, strip
	      cms.PSet( rpId = cms.uint32(0x76180000), dirName = cms.string("XRPH_D6L5_B2"), z = cms.double(-21255.1) ),  # RP 003, strip
	      cms.PSet( rpId = cms.uint32(0x77100000), dirName = cms.string("XRPH_C6R5_B1"), z = cms.double(+20382.6) ),  # RP 102, strip
	      cms.PSet( rpId = cms.uint32(0x77180000), dirName = cms.string("XRPH_D6R5_B1"), z = cms.double(+21255.1) ),  # RP 103, strip
  	)
  ),
  #geometry
  xmlIdealGeometry=cms.PSet(
	geomXMLFiles = totemGeomXMLFiles + ctppsDiamondGeomXMLFiles + ctppsUFSDGeomXMLFiles + ctppsPixelGeomXMLFiles,
	rootNodeName = cms.string('cms:CMSE')

  ),
  #alignment
  ctppsRPAlignmentCorrectionsDataXML=cms.PSet(
	MeasuredFiles=cms.vstring(),
	RealFiles=cms.vstring("Validation/CTPPS/alignment/2016_preTS2.xml"),
	MisalignedFiles=cms.vstring("Validation/CTPPS/alignment/2016_preTS2.xml")
  ),

  #direct simu data
  ctppsDirectSimuData=cms.PSet(
	useEmpiricalApertures=cms.bool(True),
	empiricalAperture45=cms.string("3.76296E-05+(([xi]<0.117122)*0.00712775+([xi]>=0.117122)*0.0148651)*([xi]-0.117122)"),
	empiricalAperture56=cms.string("1.85954E-05+(([xi]<0.14324)*0.00475349+([xi]>=0.14324)*0.00629514)*([xi]-0.14324)"),
	timeResolutionDiamonds45=cms.string("0.200"),
	timeResolutionDiamonds56=cms.string("0.200"),
	useTimeEfficiencyCheck=cms.bool(False),
	effTimePath=cms.string(""),
	effTimeObject45=cms.string(""),
	effTimeObject56=cms.string("")
  )
)

profile_2016_preTS2.xmlIdealGeometry.geomXMLFiles.append("Geometry/VeryForwardData/data/2016_ctpps_15sigma_margin0/RP_Dist_Beam_Cent.xml")

from CalibPPS.ESProducers.ctppsInterpolatedOpticalFunctionsESSource_cfi import *
ctppsInterpolatedOpticalFunctionsESSource.lhcInfoLabel = ""

def UseConstantXangleBetaStar(process, xangle, betaStar):
  process.profile_2016_preTS2.ctppsLHCInfo.xangle = xangle
  process.profile_2016_preTS2.ctppsLHCInfo.betaStar = betaStar

def UseXangleBetaStarHistogram(process,f, obj):
  process.profile_2016_preTS2.ctppsLHCInfo.xangleBetaStarHistogramFile = f
  process.profile_2016_preTS2.ctppsLHCInfo.xangleBetaStarHistogramObject = obj


