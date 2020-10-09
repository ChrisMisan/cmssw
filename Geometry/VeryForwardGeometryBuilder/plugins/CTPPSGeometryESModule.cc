/****************************************************************************
*
* Authors:
*  Jan Kaspar (jan.kaspar@gmail.com)
*  Dominik Mierzejewski <dmierzej@cern.ch>
*
****************************************************************************/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDMaterial.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "DetectorDescription/Core/interface/DDRotationMatrix.h"

#include "CondFormats/PPSObjects/interface/CTPPSRPAlignmentCorrectionsData.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSDetId/interface/TotemTimingDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "CondFormats/AlignmentRecord/interface/RPRealAlignmentRecord.h"
//#include "CondFormats/AlignmentRecord/interface/CTPPSRPAlignmentCorrectionsDataRcd.h"

#include "CondFormats/AlignmentRecord/interface/RPMisalignedAlignmentRecord.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/VeryForwardMisalignedGeometryRecord.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/DetGeomDesc.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSDDDNames.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometryESModuleCommon.h"

#include <regex>

/**
 * \brief Builds ideal, real and misaligned geometries.
 *
 * First, it creates a tree of DetGeomDesc from DDCompView. For real and misaligned geometries,
 * it applies alignment corrections (RPAlignmentCorrections) found in corresponding ...GeometryRecord.
 *
 * Second, it creates CTPPSGeometry from DetGeoDesc tree.
 **/
class CTPPSGeometryESModule : public edm::ESProducer {
public:
  CTPPSGeometryESModule(const edm::ParameterSet&);
  ~CTPPSGeometryESModule() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  template <typename ALIGNMENT_REC>
  struct GDTokens {
  explicit GDTokens(edm::ESConsumesCollector&& iCC)
      : idealGDToken_{iCC.consumesFrom<DetGeomDesc, IdealGeometryRecord>(edm::ESInputTag())},
      alignmentToken_{iCC.consumesFrom<CTPPSRPAlignmentCorrectionsData, ALIGNMENT_REC>(edm::ESInputTag())} {}
    const edm::ESGetToken<DetGeomDesc, IdealGeometryRecord> idealGDToken_;
    const edm::ESGetToken<CTPPSRPAlignmentCorrectionsData, ALIGNMENT_REC> alignmentToken_;
  };

  template <typename REC>
  std::unique_ptr<DetGeomDesc> produceGD(IdealGeometryRecord const&,
                                         const std::optional<REC>&,
                                         GDTokens<REC> const&,
                                         const char*);

  std::unique_ptr<DetGeomDesc> produceIdealGD(const IdealGeometryRecord&);
  std::unique_ptr<DetGeomDesc> produceRealGD(const VeryForwardRealGeometryRecord&);
  std::unique_ptr<CTPPSGeometry> produceRealTG(const VeryForwardRealGeometryRecord&);

  std::unique_ptr<DetGeomDesc> produceMisalignedGD(const VeryForwardMisalignedGeometryRecord&);
  std::unique_ptr<CTPPSGeometry> produceMisalignedTG(const VeryForwardMisalignedGeometryRecord&);

  const unsigned int verbosity_;
  std::unique_ptr<CTPPSGeometryESModuleCommon> ctppsGeometryESModuleCommon;
  const edm::ESGetToken<DDCompactView, IdealGeometryRecord> compactViewToken_;

  const GDTokens<RPRealAlignmentRecord> gdRealTokens_;
  const GDTokens<RPMisalignedAlignmentRecord> gdMisTokens_;

  const edm::ESGetToken<DetGeomDesc, VeryForwardRealGeometryRecord> dgdRealToken_;
  const edm::ESGetToken<DetGeomDesc, VeryForwardMisalignedGeometryRecord> dgdMisToken_;
};

//----------------------------------------------------------------------------------------------------

CTPPSGeometryESModule::CTPPSGeometryESModule(const edm::ParameterSet& iConfig)
    : verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity")),
      compactViewToken_{setWhatProduced(this, &CTPPSGeometryESModule::produceIdealGD)
                            .consumes<DDCompactView>(edm::ESInputTag(
                                "" /*optional module label */, iConfig.getParameter<std::string>("compactViewTag")))},
      gdRealTokens_{setWhatProduced(this, &CTPPSGeometryESModule::produceRealGD)},
      gdMisTokens_{setWhatProduced(this, &CTPPSGeometryESModule::produceMisalignedGD)},
      dgdRealToken_{
          setWhatProduced(this, &CTPPSGeometryESModule::produceRealTG).consumes<DetGeomDesc>(edm::ESInputTag())},
      dgdMisToken_{
          setWhatProduced(this, &CTPPSGeometryESModule::produceMisalignedTG).consumes<DetGeomDesc>(edm::ESInputTag())} {
  ctppsGeometryESModuleCommon=std::make_unique<CTPPSGeometryESModuleCommon>(iConfig);
}

void CTPPSGeometryESModule::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<unsigned int>("verbosity", 1);
  desc.add<std::string>("compactViewTag", std::string());
  descriptions.add("DoodadESSource", desc);
}

//----------------------------------------------------------------------------------------------------
template <typename REC>
std::unique_ptr<DetGeomDesc> CTPPSGeometryESModule::produceGD(IdealGeometryRecord const& iIdealRec,
                                                              std::optional<REC> const& iAlignRec,
                                                              GDTokens<REC> const& iTokens,
                                                              const char* name) {
  // get the input GeometricalDet
  auto const& idealGD = iIdealRec.get(iTokens.idealGDToken_);

  // load alignments
  edm::ESHandle<CTPPSRPAlignmentCorrectionsData> alignments;
  if (iAlignRec) {
    alignments = iAlignRec->getHandle(iTokens.alignmentToken_);
  }

  if (alignments.isValid()) {
    if (verbosity_)
      edm::LogVerbatim(name) << ">> " << name << " > Real geometry: " << alignments->getRPMap().size() << " RP and "
                             << alignments->getSensorMap().size() << " sensor alignments applied.";
  } else {
    if (verbosity_)
      edm::LogVerbatim(name) << ">> " << name << " > Real geometry: No alignments applied.";
  }

  DetGeomDesc* newGD = nullptr;
  ctppsGeometryESModuleCommon->applyAlignments(idealGD, alignments.product(), newGD);
  return std::unique_ptr<DetGeomDesc>(newGD);
}
//----------------------------------------------------------------------------------------------------
std::unique_ptr<DetGeomDesc> CTPPSGeometryESModule::produceIdealGD(const IdealGeometryRecord& iRecord) {
  // get the DDCompactView from EventSetup
  auto const& cpv = iRecord.get(compactViewToken_);

  // create DDFilteredView and apply the filter
  DDPassAllFilter filter;
  DDFilteredView fv(cpv, filter);

  // conversion to DetGeomDesc structure
  auto root = std::make_unique<DetGeomDesc>(&fv);
  ctppsGeometryESModuleCommon->buildDetGeomDesc(&fv, root.get());

  // construct the tree of DetGeomDesc
  return root;
}

//----------------------------------------------------------------------------------------------------

std::unique_ptr<DetGeomDesc> CTPPSGeometryESModule::produceRealGD(const VeryForwardRealGeometryRecord& iRecord) {
  return produceGD(iRecord.getRecord<IdealGeometryRecord>(),
                   iRecord.tryToGetRecord<RPRealAlignmentRecord>(),
                   gdRealTokens_,
                   "CTPPSGeometryESModule::produceRealGD");
}

//----------------------------------------------------------------------------------------------------

std::unique_ptr<DetGeomDesc> CTPPSGeometryESModule::produceMisalignedGD(
    const VeryForwardMisalignedGeometryRecord& iRecord) {
  return produceGD(iRecord.getRecord<IdealGeometryRecord>(),
                   iRecord.tryToGetRecord<RPMisalignedAlignmentRecord>(),
                   gdMisTokens_,
                   "CTPPSGeometryESModule::produceMisalignedGD");
}

//----------------------------------------------------------------------------------------------------

std::unique_ptr<CTPPSGeometry> CTPPSGeometryESModule::produceRealTG(const VeryForwardRealGeometryRecord& iRecord) {
  auto const& gD = iRecord.get(dgdRealToken_);

  return std::make_unique<CTPPSGeometry>(&gD);
}

//----------------------------------------------------------------------------------------------------

std::unique_ptr<CTPPSGeometry> CTPPSGeometryESModule::produceMisalignedTG(
    const VeryForwardMisalignedGeometryRecord& iRecord) {
  auto const& gD = iRecord.get(dgdMisToken_);

  return std::make_unique<CTPPSGeometry>(&gD);
}

DEFINE_FWK_EVENTSETUP_MODULE(CTPPSGeometryESModule);
