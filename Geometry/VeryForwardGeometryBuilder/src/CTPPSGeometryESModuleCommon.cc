#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometryESModuleCommon.h"


CTPPSGeometryESModuleCommon::CTPPSGeometryESModuleCommon(const edm::ParameterSet& iConfig)
    : verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity")){}

//----------------------------------------------------------------------------------------------------

void CTPPSGeometryESModuleCommon::applyAlignments(const DetGeomDesc& idealGD,
                                            const CTPPSRPAlignmentCorrectionsData* alignments,
                                            DetGeomDesc*& newGD) {
  newGD = new DetGeomDesc(idealGD);
  std::deque<const DetGeomDesc*> buffer;
  std::deque<DetGeomDesc*> bufferNew;
  buffer.emplace_back(&idealGD);
  bufferNew.emplace_back(newGD);

  while (!buffer.empty()) {
    const DetGeomDesc* sD = buffer.front();
    DetGeomDesc* pD = bufferNew.front();
    buffer.pop_front();
    bufferNew.pop_front();

    const std::string name = pD->name();

    // Is it sensor? If yes, apply full sensor alignments
    if (name == DDD_TOTEM_RP_SENSOR_NAME || name == DDD_CTPPS_DIAMONDS_SEGMENT_NAME ||
        name == DDD_CTPPS_UFSD_SEGMENT_NAME || name == DDD_CTPPS_PIXELS_SENSOR_NAME ||
        std::regex_match(name, std::regex(DDD_TOTEM_TIMING_SENSOR_TMPL))) {
      unsigned int plId = pD->geographicalID();

      if (alignments) {
        const auto& ac = alignments->getFullSensorCorrection(plId);
        pD->applyAlignment(ac);
      }
    }

    // Is it RP box? If yes, apply RP alignments
    if (name == DDD_TOTEM_RP_RP_NAME || name == DDD_CTPPS_DIAMONDS_RP_NAME || name == DDD_CTPPS_PIXELS_RP_NAME ||
        name == DDD_TOTEM_TIMING_RP_NAME) {
      unsigned int rpId = pD->geographicalID();

      if (alignments) {
        const auto& ac = alignments->getRPCorrection(rpId);
        pD->applyAlignment(ac);
      }
    }

    // create and add children
    for (unsigned int i = 0; i < sD->components().size(); i++) {
      const DetGeomDesc* sDC = sD->components()[i];
      buffer.emplace_back(sDC);

      // create new node with the same information as in sDC and add it as a child of pD
      DetGeomDesc* cD = new DetGeomDesc(*sDC);
      pD->addComponent(cD);

      bufferNew.emplace_back(cD);
    }
  }
}

//----------------------------------------------------------------------------------------------------

void CTPPSGeometryESModuleCommon::buildDetGeomDesc(DDFilteredView* fv, DetGeomDesc* gd) {
  // try to dive into next level
  if (!fv->firstChild())
    return;

  // loop over siblings in the level
  do {
    // create new DetGeomDesc node and add it to the parent's (gd) list
    DetGeomDesc* newGD = new DetGeomDesc(fv);

    const std::string name = fv->logicalPart().name().name();

    // strip sensors
    if (name == DDD_TOTEM_RP_SENSOR_NAME) {
      const std::vector<int>& copy_num = fv->copyNumbers();
      // check size of copy numubers array
      if (copy_num.size() < 3)
        throw cms::Exception("DDDTotemRPContruction")
            << "size of copyNumbers for strip sensor is " << copy_num.size() << ". It must be >= 3.";

      // extract information
      const unsigned int decRPId = copy_num[copy_num.size() - 3];
      const unsigned int arm = decRPId / 100;
      const unsigned int station = (decRPId % 100) / 10;
      const unsigned int rp = decRPId % 10;
      const unsigned int detector = copy_num[copy_num.size() - 1];
      newGD->setGeographicalID(TotemRPDetId(arm, station, rp, detector));
    }

    // strip and pixels RPs
    else if (name == DDD_TOTEM_RP_RP_NAME || name == DDD_CTPPS_PIXELS_RP_NAME) {
      unsigned int decRPId = fv->copyno();

      // check if it is a pixel RP
      if (decRPId >= 10000) {
        decRPId = decRPId % 10000;
        const unsigned int armIdx = (decRPId / 100) % 10;
        const unsigned int stIdx = (decRPId / 10) % 10;
        const unsigned int rpIdx = decRPId % 10;
        newGD->setGeographicalID(CTPPSPixelDetId(armIdx, stIdx, rpIdx));
      } else {
        const unsigned int armIdx = (decRPId / 100) % 10;
        const unsigned int stIdx = (decRPId / 10) % 10;
        const unsigned int rpIdx = decRPId % 10;
        newGD->setGeographicalID(TotemRPDetId(armIdx, stIdx, rpIdx));
      }
    }

    else if (std::regex_match(name, std::regex(DDD_TOTEM_TIMING_SENSOR_TMPL))) {
      const std::vector<int>& copy_num = fv->copyNumbers();
      // check size of copy numbers array
      if (copy_num.size() < 4)
        throw cms::Exception("DDDTotemRPContruction")
            << "size of copyNumbers for TOTEM timing sensor is " << copy_num.size() << ". It must be >= 4.";

      const unsigned int decRPId = copy_num[copy_num.size() - 4];
      const unsigned int arm = decRPId / 100, station = (decRPId % 100) / 10, rp = decRPId % 10;
      const unsigned int plane = copy_num[copy_num.size() - 2], channel = copy_num[copy_num.size() - 1];
      newGD->setGeographicalID(TotemTimingDetId(arm, station, rp, plane, channel));
    }

    else if (name == DDD_TOTEM_TIMING_RP_NAME) {
      const unsigned int arm = fv->copyno() / 100, station = (fv->copyno() % 100) / 10, rp = fv->copyno() % 10;
      newGD->setGeographicalID(TotemTimingDetId(arm, station, rp));
    }

    // pixel sensors
    else if (name == DDD_CTPPS_PIXELS_SENSOR_NAME) {
      const std::vector<int>& copy_num = fv->copyNumbers();
      // check size of copy numubers array
      if (copy_num.size() < 4)
        throw cms::Exception("DDDTotemRPContruction")
            << "size of copyNumbers for pixel sensor is " << copy_num.size() << ". It must be >= 4.";

      // extract information
      const unsigned int decRPId = copy_num[copy_num.size() - 4] % 10000;
      const unsigned int arm = decRPId / 100;
      const unsigned int station = (decRPId % 100) / 10;
      const unsigned int rp = decRPId % 10;
      const unsigned int detector = copy_num[copy_num.size() - 2] - 1;
      newGD->setGeographicalID(CTPPSPixelDetId(arm, station, rp, detector));
    }

    // diamond/UFSD sensors
    else if (name == DDD_CTPPS_DIAMONDS_SEGMENT_NAME || name == DDD_CTPPS_UFSD_SEGMENT_NAME) {
      const std::vector<int>& copy_num = fv->copyNumbers();

      const unsigned int id = copy_num[copy_num.size() - 1];
      const unsigned int arm = copy_num[1] - 1;
      const unsigned int station = 1;
      const unsigned int rp = 6;
      const unsigned int plane = (id / 100);
      const unsigned int channel = id % 100;

      newGD->setGeographicalID(CTPPSDiamondDetId(arm, station, rp, plane, channel));
    }

    // diamond/UFSD RPs
    else if (name == DDD_CTPPS_DIAMONDS_RP_NAME) {
      const std::vector<int>& copy_num = fv->copyNumbers();

      // check size of copy numubers array
      if (copy_num.size() < 2)
        throw cms::Exception("DDDTotemRPContruction")
            << "size of copyNumbers for diamond RP is " << copy_num.size() << ". It must be >= 2.";

      const unsigned int arm = copy_num[1] - 1;
      const unsigned int station = 1;
      const unsigned int rp = 6;

      newGD->setGeographicalID(CTPPSDiamondDetId(arm, station, rp));
    }

    // add component
    gd->addComponent(newGD);

    // recursion
    buildDetGeomDesc(fv, newGD);
  } while (fv->nextSibling());

  // go a level up
  fv->parent();
}

