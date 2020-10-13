/****************************************************************************
*
* Authors:
*  Jan Kaspar (jan.kaspar@gmail.com)
*  Christopher Misan (krzysmisan@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"
#include "CondFormats/PPSObjects/interface/LHCOpticalFunctionsSetCollection.h"
#include "CondFormats/DataRecord/interface/CTPPSOpticsRcd.h"
#include "GeometryReaders/XMLIdealGeometryESSource/interface/GeometryConfiguration.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DetectorDescription/Parser/interface/DDLParser.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDRoot.h"
#include "DetectorDescription/Core/interface/DDMaterial.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "DetectorDescription/Core/interface/DDRotationMatrix.h"
#include "DetectorDescription/Core/src/Material.h"
#include "DetectorDescription/Core/src/Solid.h"
#include "DetectorDescription/Core/src/LogicalPart.h"
#include "DetectorDescription/Core/src/Specific.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/DetGeomDesc.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSDDDNames.h"
#include "CondFormats/PPSObjects/interface/CTPPSRPAlignmentCorrectionsData.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSDetId/interface/TotemTimingDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "Geometry/Records/interface/VeryForwardMisalignedGeometryRecord.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "CondFormats/PPSObjects/interface/CTPPSRPAlignmentCorrectionsDataSequence.h"
#include "CondFormats/PPSObjects/interface/CTPPSRPAlignmentCorrectionsMethods.h"
#include "CondFormats/AlignmentRecord/interface/CTPPSRPAlignmentCorrectionsDataRcd.h"  // this used to be RPMeasuredAlignmentRecord.h
#include "CondFormats/AlignmentRecord/interface/RPRealAlignmentRecord.h"
#include "CondFormats/AlignmentRecord/interface/RPMisalignedAlignmentRecord.h"
#include "CondFormats/PPSObjects/interface/CTPPSDirectSimuData.h"
#include "CondFormats/DataRecord/interface/CTPPSDirectSimuDataRcd.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/JamesRandom.h"
#include <memory>
#include <regex>
#include <vector>
#include <string>
#include <map>
#include <set>
#include "TRandom.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "CalibPPS/ESProducers/plugins/CTPPSRPAlignmentCorrectionsDataESSourceXMLCommon.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometryESModuleCommon.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/DetGeomDescBuilder.h"

//----------------------------------------------------------------------------------------------------

class CTPPSCompositeESSource : public edm::ESProducer, public edm::EventSetupRecordIntervalFinder {
public:
  CTPPSCompositeESSource(const edm::ParameterSet &);
  static void fillDescriptions(edm::ConfigurationDescriptions &);
  std::unique_ptr<LHCInfo> produceLhcInfo(const LHCInfoRcd &);
  std::unique_ptr<LHCOpticalFunctionsSetCollection> produceOptics(const CTPPSOpticsRcd &);
  std::shared_ptr<DDCompactView> produceGeom(const IdealGeometryRecord &);
  std::shared_ptr<CTPPSRPAlignmentCorrectionsData> produceMeasured(const CTPPSRPAlignmentCorrectionsDataRcd &);
  std::shared_ptr<CTPPSRPAlignmentCorrectionsData> produceReal(const RPRealAlignmentRecord &);
  std::shared_ptr<CTPPSRPAlignmentCorrectionsData> produceMisaligned(const RPMisalignedAlignmentRecord &);
  std::shared_ptr<DetGeomDesc> produceIdealGD( const IdealGeometryRecord& );
  std::shared_ptr<CTPPSGeometry> produceRealTG( const VeryForwardRealGeometryRecord& );
  std::shared_ptr<CTPPSGeometry> produceMisalignedTG( const VeryForwardMisalignedGeometryRecord& );
  std::unique_ptr<CTPPSDirectSimuData> produceDirectSimuData(const CTPPSDirectSimuDataRcd &);

private:
  void setIntervalFor(const edm::eventsetup::EventSetupRecordKey &,
                      const edm::IOVSyncValue &,
                      edm::ValidityInterval &) override;

  struct FileInfo {
    double m_xangle;
    std::string m_fileName;
  };

  struct RPInfo {
    std::string m_dirName;
    double m_scoringPlaneZ;
  };

  struct Entry {
    std::vector<FileInfo> m_fileInfo;
    std::unordered_map<unsigned int, RPInfo> m_rpInfo;
  };

  template <typename T>
  struct BinData {
    double min, max;
    T data;
  };

  struct profileData{
    //lhcInfo
    double xangle,betaStar;
    TH2D *h_xangle_beta_star;
    std::vector<BinData<std::pair<double,double>>> xangleBetaStarBins;
    double m_beamEnergy;

    //optics
    LHCOpticalFunctionsSetCollection lhcOptical;

    //geometry
    std::shared_ptr<DDCompactView> ddCompactView;
    std::shared_ptr<DetGeomDesc> idealGD;
    std::shared_ptr<DetGeomDesc> misalignedGD;
    std::shared_ptr<DetGeomDesc> realGD;
    std::shared_ptr<CTPPSGeometry> misalignedTG;
    std::shared_ptr<CTPPSGeometry> realTG;

    //alignmentacMisaligned
    std::shared_ptr<CTPPSRPAlignmentCorrectionsData> acMeasured, acReal, acMisaligned;

    //data pass
    CTPPSDirectSimuData directSimuData;
  };

  static void applyAlignments(const DetGeomDesc&, const CTPPSRPAlignmentCorrectionsData*, DetGeomDesc*&);
  void buildDetGeomDesc(DDFilteredView* fv, DetGeomDesc* gd);
  void buildLhcInfo(const edm::ParameterSet& profile,profileData& pData);
  void buildOptics(const edm::ParameterSet& profile,profileData& pData);
  void buildCompactView(const edm::ParameterSet& profile,profileData& pData);
  void buildDirectSimuData(const edm::ParameterSet& profile,profileData& pData);


  std::string lhcInfoLabel_;
  std::string opticsLabel_;
  unsigned int m_generateEveryNEvents_;
  std::unique_ptr<CLHEP::HepRandomEngine> m_engine_;
  unsigned int verbosity_;
  std::vector<profileData> profiles_;
  std::vector<BinData<profileData&>> weights_;
  profileData* currentProfile_;
  std::unique_ptr<CTPPSGeometryESModuleCommon> ctppsGeometryESModuleCommon;

};

//----------------------------------------------------------------------------------------------------

CTPPSCompositeESSource::CTPPSCompositeESSource(const edm::ParameterSet &conf)
    : lhcInfoLabel_(conf.getParameter<std::string>("lhcInfoLabel")),

    opticsLabel_(conf.getParameter<std::string>("opticsLabel")),

    m_generateEveryNEvents_(conf.getParameter<unsigned int>("generateEveryNEvents")),

    m_engine_(new CLHEP::HepJamesRandom(conf.getParameter<unsigned int>("seed"))),
    verbosity_(conf.getUntrackedParameter<unsigned int>("verbosity")){
      
  gRandom->SetSeed(0);

  double s=0;
  for(const auto &profile:conf.getParameter<std::vector<edm::ParameterSet>>("periods")){
    profiles_.push_back(profileData());
    auto& pData=profiles_.back();
    s+=profile.getParameter<double>("L_i");

    auto ctppsRPAlignmentCorrectionsDataXMLpSet=profile.getParameter<edm::ParameterSet>("ctppsRPAlignmentCorrectionsDataXML");
    ctppsRPAlignmentCorrectionsDataXMLpSet.addUntrackedParameter("verbosity",verbosity_);
    CTPPSRPAlignmentCorrectionsDataESSourceXMLCommon ctppsRPAlignmentCorrectionsDataESSourceXMLCommon(ctppsRPAlignmentCorrectionsDataXMLpSet);
  if(ctppsRPAlignmentCorrectionsDataXMLpSet.getParameter<std::vector<std::string> >("MisalignedFiles").empty())
      pData.acMeasured=std::make_shared<CTPPSRPAlignmentCorrectionsData>(ctppsRPAlignmentCorrectionsDataESSourceXMLCommon.acsMeasured[0].second);
    pData.acReal=std::make_shared<CTPPSRPAlignmentCorrectionsData>(ctppsRPAlignmentCorrectionsDataESSourceXMLCommon.acsReal[0].second);
    pData.acMisaligned=std::make_shared<CTPPSRPAlignmentCorrectionsData>(ctppsRPAlignmentCorrectionsDataESSourceXMLCommon.acsMisaligned[0].second);

    buildLhcInfo(profile,pData);
    buildOptics(profile,pData);
    buildCompactView(profile,pData);
    buildDirectSimuData(profile,pData);
    
    
    auto const& cpv = *(pData.ddCompactView);
    pData.idealGD= std::move(detgeomdescbuilder::buildDetGeomDescFromCompactView(cpv));
  
  
    pData.misalignedGD=std::move(ctppsGeometryESModuleCommon->applyAlignments(*(pData.idealGD), pData.acMisaligned.get()));
    pData.misalignedTG=std::make_shared<CTPPSGeometry>(pData.misalignedGD.get(),verbosity_);



    pData.realGD=std::move(ctppsGeometryESModuleCommon->applyAlignments(*(pData.idealGD), pData.acReal.get()));
    pData.realTG=std::make_shared<CTPPSGeometry>(pData.realGD.get(),verbosity_);



    //calc xangle/beta*
    if(pData.xangle<0){
      double sum = 0.;
      for (int bi = 1; bi <= pData.h_xangle_beta_star->GetNcells(); ++bi){
        double val=pData.h_xangle_beta_star->GetBinContent(bi);
        sum+=val;
      }
      
      double cw=0;
      for (int x = 1; x <= pData.h_xangle_beta_star->GetNbinsX(); ++x) 
        for (int y = 1; y <= pData.h_xangle_beta_star->GetNbinsY(); ++y) {
          double sample=pData.h_xangle_beta_star->GetBinContent(pData.h_xangle_beta_star->GetBin(x,y));
          if(sample>=1){
            sample/=sum;
            pData.xangleBetaStarBins.push_back({cw,cw+sample,std::pair<double,double>(pData.h_xangle_beta_star->GetXaxis()->GetBinCenter(x),pData.h_xangle_beta_star->GetYaxis()->GetBinCenter(y))});
            //edm::LogVerbatim("CTPPSCompositeESSource::produceLhcInfo")<<"x "<<pData.h_xangle_beta_star->GetXaxis()->GetBinCenter(x)<<" y "<<pData.h_xangle_beta_star->GetYaxis()->GetBinCenter(y)<<" w "<<sample<<std::endl;
            cw+=sample;
          }
        }
      delete pData.h_xangle_beta_star;
    }
  }
  //calc weights_
  double cw = 0.;
  int counter=0;
  for(const auto &profile:conf.getParameter<std::vector<edm::ParameterSet>>("periods")){
    double w=profile.getParameter<double>("L_i")/s;
    weights_.push_back({cw, cw + w, profiles_[counter]});
    counter++;
    cw+=w;
  }

  setWhatProduced(this, &CTPPSCompositeESSource::produceLhcInfo,edm::es::Label(lhcInfoLabel_));
  setWhatProduced(this, &CTPPSCompositeESSource::produceOptics,edm::es::Label(opticsLabel_));
  setWhatProduced(this, &CTPPSCompositeESSource::produceGeom,edm::es::Label("XMLIdealGeometryESSource_CTPPS"));
  setWhatProduced(this, &CTPPSCompositeESSource::produceIdealGD);
  setWhatProduced(this, &CTPPSCompositeESSource::produceRealTG);
  setWhatProduced(this, &CTPPSCompositeESSource::produceMisalignedTG);
  setWhatProduced(this, &CTPPSCompositeESSource::produceMeasured);
  setWhatProduced(this, &CTPPSCompositeESSource::produceReal);
  setWhatProduced(this, &CTPPSCompositeESSource::produceMisaligned);
  setWhatProduced(this, &CTPPSCompositeESSource::produceDirectSimuData);
  findingRecord<LHCInfoRcd>();
  findingRecord<CTPPSOpticsRcd>();
  findingRecord<IdealGeometryRecord>();
  findingRecord<CTPPSRPAlignmentCorrectionsDataRcd>();
  findingRecord<RPRealAlignmentRecord>();
  findingRecord<RPMisalignedAlignmentRecord>();
  findingRecord<CTPPSDirectSimuDataRcd>();

}
//----------------------------------------------------------------------------------------------------
void CTPPSCompositeESSource::buildDirectSimuData(const edm::ParameterSet& profile,profileData& pData){
  const auto& ctppsDirectSimuData=profile.getParameter<edm::ParameterSet>("ctppsDirectSimuData");
  pData.directSimuData.setUseEmpiricalApertures(ctppsDirectSimuData.getParameter<bool>("useEmpiricalApertures"));
  pData.directSimuData.setEmpiricalAperture45(ctppsDirectSimuData.getParameter<std::string>("empiricalAperture45"));
  pData.directSimuData.setEmpiricalAperture56(ctppsDirectSimuData.getParameter<std::string>("empiricalAperture56"));
  pData.directSimuData.setTimeResolutionDiamonds45(ctppsDirectSimuData.getParameter<std::string>("timeResolutionDiamonds45"));
  pData.directSimuData.setTimeResolutionDiamonds56(ctppsDirectSimuData.getParameter<std::string>("timeResolutionDiamonds56"));
  pData.directSimuData.setUseTimeEfficiencyCheck(ctppsDirectSimuData.getParameter<bool>("useTimeEfficiencyCheck"));
  pData.directSimuData.setEffTimePath(ctppsDirectSimuData.getParameter<std::string>("effTimePath"));
  pData.directSimuData.setEffTimeObject45(ctppsDirectSimuData.getParameter<std::string>("effTimeObject45"));
  pData.directSimuData.setEffTimeObject56(ctppsDirectSimuData.getParameter<std::string>("effTimeObject56"));

}
//----------------------------------------------------------------------------------------------------
void CTPPSCompositeESSource::buildCompactView(const edm::ParameterSet& profile,profileData& pData){
  const auto& xmlIdealGeometry=profile.getParameter<edm::ParameterSet>("xmlIdealGeometry");

  DDName ddName(xmlIdealGeometry.getParameter<std::string>("rootNodeName"));
  DDLogicalPart rootNode(ddName);
  DDRootDef::instance().set(rootNode);
  std::shared_ptr<DDCompactView> returnValue(new DDCompactView(rootNode));
  DDLParser parser(*returnValue);  //* parser = DDLParser::instance();
  parser.getDDLSAX2FileHandler()->setUserNS(xmlIdealGeometry.getUntrackedParameter<bool>("userControlledNamespace", false));
  int result2 = parser.parse(*(std::make_unique<GeometryConfiguration>(xmlIdealGeometry)));
  if (result2 != 0)
    throw cms::Exception("DDException") << "DDD-Parser: parsing failed!";

  // after parsing the root node should be valid!

  if (!rootNode.isValid()) 
    throw cms::Exception("Geometry") << "There is no valid node named \"" << xmlIdealGeometry.getParameter<std::string>("rootNodeName") << "\"";
  
  returnValue->lockdown();
  pData.ddCompactView=returnValue;
}
//----------------------------------------------------------------------------------------------------
void CTPPSCompositeESSource::buildOptics(const edm::ParameterSet& profile,profileData& pData){
  const auto& ctppsOpticalFunctions=profile.getParameter<edm::ParameterSet>("ctppsOpticalFunctions");
  std::vector<FileInfo> fileInfo;
  for (const auto &pset : ctppsOpticalFunctions.getParameter<std::vector<edm::ParameterSet>>("opticalFunctions")) {
    const double &xangle = pset.getParameter<double>("xangle");
    const std::string &fileName = pset.getParameter<edm::FileInPath>("fileName").fullPath();
    fileInfo.push_back({xangle, fileName});
  }

  std::unordered_map<unsigned int, RPInfo> rpInfo;
  for (const auto &pset : ctppsOpticalFunctions.getParameter<std::vector<edm::ParameterSet>>("scoringPlanes")) {
    const unsigned int rpId = pset.getParameter<unsigned int>("rpId");
    const std::string dirName = pset.getParameter<std::string>("dirName");
    const double z = pset.getParameter<double>("z");
    const RPInfo entry = {dirName, z};
    rpInfo.emplace(rpId, entry);
  }

  Entry entry({fileInfo, rpInfo});

  for (const auto &fi : entry.m_fileInfo) {
    std::unordered_map<unsigned int, LHCOpticalFunctionsSet> xa_data;

    for (const auto &rpi : entry.m_rpInfo) {
      LHCOpticalFunctionsSet fcn(fi.m_fileName, rpi.second.m_dirName, rpi.second.m_scoringPlaneZ);
      xa_data.emplace(rpi.first, std::move(fcn));
    }

    pData.lhcOptical.emplace(fi.m_xangle, xa_data);
  }
}
//----------------------------------------------------------------------------------------------------
void CTPPSCompositeESSource::buildLhcInfo(const edm::ParameterSet& profile,profileData& pData){
  //LHCInfo
  const auto& ctppsLHCInfo=profile.getParameter<edm::ParameterSet>("ctppsLHCInfo");
  pData.xangle=ctppsLHCInfo.getParameter<double>("xangle");
  pData.betaStar=ctppsLHCInfo.getParameter<double>("betaStar");

  pData.m_beamEnergy=ctppsLHCInfo.getParameter<double>("beamEnergy");

  if(pData.xangle>0)
    return;

  edm::FileInPath fip(ctppsLHCInfo.getParameter<std::string>("xangleBetaStarHistogramFile").c_str());
  std::unique_ptr<TFile> f_in(TFile::Open(fip.fullPath().c_str()));
  if (!f_in)
    throw cms::Exception("PPS") << "Cannot open input file '" << ctppsLHCInfo.getParameter<std::string>("xangleBetaStarHistogramFile") << "'.";

  pData.h_xangle_beta_star = (TH2D *)f_in->Get(ctppsLHCInfo.getParameter<std::string>("xangleBetaStarHistogramObject").c_str());
  if (!pData.h_xangle_beta_star)
    throw cms::Exception("PPS") << "Cannot load input object '" << ctppsLHCInfo.getParameter<std::string>("xangleBetaStarHistogramObject") << "'.";

  pData.h_xangle_beta_star->SetDirectory(0);
  
}

//----------------------------------------------------------------------------------------------------

void CTPPSCompositeESSource::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("lhcInfoLabel", "")->setComment("label of the LHCInfo record");
  desc.add<std::string>("opticsLabel", "")->setComment("label of the optics record");
  desc.add<unsigned int>("seed", 1)->setComment("random seed");
  desc.add<unsigned int>("generateEveryNEvents", 1)->setComment("how often to generate new xangle");
  desc.addUntracked<unsigned int>("verbosity",0);

  edm::ParameterSetDescription desc_profile;
  std::vector<edm::ParameterSet> vp;
  desc_profile.add<double>("L_i", 0.)->setComment("integrated luminosity");

  //lhcInfo
  edm::ParameterSetDescription desc_profile_ctppsLHCInfo;
  desc_profile_ctppsLHCInfo.add<double>("xangle", 0.)->setComment("constant xangle");
  desc_profile_ctppsLHCInfo.add<double>("betaStar", 0.)->setComment("constant betaStar");
  desc_profile_ctppsLHCInfo.add<double>("beamEnergy", 0.)->setComment("beam energy");
  desc_profile_ctppsLHCInfo.add<std::string>("xangleBetaStarHistogramFile", "")->setComment("ROOT file with xangle/beta* distribution");
  desc_profile_ctppsLHCInfo.add<std::string>("xangleBetaStarHistogramObject", "")->setComment("xangle distribution object in the ROOT file");
  desc_profile.add<edm::ParameterSetDescription>("ctppsLHCInfo",desc_profile_ctppsLHCInfo);

  //optics
  edm::ParameterSetDescription desc_profile_ctppsOpticalFunctions;
  edm::ParameterSetDescription of_desc;
  of_desc.add<double>("xangle")->setComment("half crossing angle value in urad");
  of_desc.add<edm::FileInPath>("fileName")->setComment("ROOT file with optical functions");
  std::vector<edm::ParameterSet> of;
  desc_profile_ctppsOpticalFunctions.addVPSet("opticalFunctions", of_desc, of)
      ->setComment("list of optical functions at different crossing angles");

  edm::ParameterSetDescription sp_desc;
  sp_desc.add<unsigned int>("rpId")->setComment("associated detector DetId");
  sp_desc.add<std::string>("dirName")->setComment("associated path to the optical functions file");
  sp_desc.add<double>("z")->setComment("longitudinal position at scoring plane/detector");
  std::vector<edm::ParameterSet> sp;
  desc_profile_ctppsOpticalFunctions.addVPSet("scoringPlanes", sp_desc, sp)->setComment("list of sensitive planes/detectors stations");
  desc_profile.add<edm::ParameterSetDescription>("ctppsOpticalFunctions",desc_profile_ctppsOpticalFunctions);

  //geometry
  edm::ParameterSetDescription desc_profile_xmlIdealGeometry;
  desc_profile_xmlIdealGeometry.add<std::vector<std::string>>("geomXMLFiles");
  desc_profile_xmlIdealGeometry.add<std::string>("rootNodeName");
  desc_profile.add<edm::ParameterSetDescription>("xmlIdealGeometry",desc_profile_xmlIdealGeometry);

  //alignment
  edm::ParameterSetDescription desc_profile_ctppsRPAlignmentCorrectionsDataXML;
  desc_profile_ctppsRPAlignmentCorrectionsDataXML.add<std::vector<std::string>>("MeasuredFiles");
  desc_profile_ctppsRPAlignmentCorrectionsDataXML.add<std::vector<std::string>>("RealFiles");
  desc_profile_ctppsRPAlignmentCorrectionsDataXML.add<std::vector<std::string>>("MisalignedFiles");
  desc_profile.add<edm::ParameterSetDescription>("ctppsRPAlignmentCorrectionsDataXML",desc_profile_ctppsRPAlignmentCorrectionsDataXML);

  //direct simu data
  edm::ParameterSetDescription desc_profile_ctppsDirectSimuData;
  desc_profile_ctppsDirectSimuData.add<bool>("useEmpiricalApertures",false);
  desc_profile_ctppsDirectSimuData.add<std::string>("empiricalAperture45");
  desc_profile_ctppsDirectSimuData.add<std::string>("empiricalAperture56");
  desc_profile_ctppsDirectSimuData.add<std::string>("timeResolutionDiamonds45");
  desc_profile_ctppsDirectSimuData.add<std::string>("timeResolutionDiamonds56");

  desc_profile_ctppsDirectSimuData.add<bool>("useTimeEfficiencyCheck",false);
  desc_profile_ctppsDirectSimuData.add<std::string>("effTimePath");
  desc_profile_ctppsDirectSimuData.add<std::string>("effTimeObject45");
  desc_profile_ctppsDirectSimuData.add<std::string>("effTimeObject56");


  desc_profile.add<edm::ParameterSetDescription>("ctppsDirectSimuData",desc_profile_ctppsDirectSimuData);

  desc.addVPSet("periods", desc_profile,vp)->setComment("profiles");

  

  descriptions.add("CTPPSCompositeESSource", desc);
}

//----------------------------------------------------------------------------------------------------

void CTPPSCompositeESSource::setIntervalFor(const edm::eventsetup::EventSetupRecordKey &key,
                                                      const edm::IOVSyncValue &iosv,
                                                      edm::ValidityInterval &oValidity) {
  edm::EventID beginEvent = iosv.eventID();
  edm::EventID endEvent(beginEvent.run(), beginEvent.luminosityBlock(), beginEvent.event() + m_generateEveryNEvents_);
  oValidity = edm::ValidityInterval(edm::IOVSyncValue(beginEvent), edm::IOVSyncValue(endEvent));

  if (strcmp(key.name(), "LHCInfoRcd") == 0){
    const double u = CLHEP::RandFlat::shoot(m_engine_.get(), 0., 1.);

    for (const auto &d : weights_) {
      if (d.min <= u && u <= d.max) {
        currentProfile_ = &d.data;
        break;
      }
    }
  }

}
//----------------------------------------------------------------------------------------------------
std::shared_ptr<CTPPSGeometry> CTPPSCompositeESSource::produceRealTG( const VeryForwardRealGeometryRecord& ){
  return currentProfile_->realTG;
}
//----------------------------------------------------------------------------------------------------
std::shared_ptr<CTPPSGeometry> CTPPSCompositeESSource::produceMisalignedTG( const VeryForwardMisalignedGeometryRecord& ){
  return currentProfile_->misalignedTG;
}
//----------------------------------------------------------------------------------------------------
std::shared_ptr<DetGeomDesc> CTPPSCompositeESSource::produceIdealGD(const IdealGeometryRecord& iRecord){
  return currentProfile_->idealGD;
}
//----------------------------------------------------------------------------------------------------
std::unique_ptr<CTPPSDirectSimuData> CTPPSCompositeESSource::produceDirectSimuData(const CTPPSDirectSimuDataRcd &){
  return std::make_unique<CTPPSDirectSimuData>(currentProfile_->directSimuData);
}
//----------------------------------------------------------------------------------------------------
std::shared_ptr<CTPPSRPAlignmentCorrectionsData> CTPPSCompositeESSource::produceMeasured(
    const CTPPSRPAlignmentCorrectionsDataRcd &iRecord) {
  return currentProfile_->acMeasured;
}
//----------------------------------------------------------------------------------------------------
std::shared_ptr<CTPPSRPAlignmentCorrectionsData> CTPPSCompositeESSource::produceReal(
    const RPRealAlignmentRecord &iRecord) {
  return currentProfile_->acReal;
}
//----------------------------------------------------------------------------------------------------
std::shared_ptr<CTPPSRPAlignmentCorrectionsData> CTPPSCompositeESSource::produceMisaligned(
    const RPMisalignedAlignmentRecord &iRecord) {
  return currentProfile_->acMisaligned;
}
//----------------------------------------------------------------------------------------------------
std::shared_ptr<DDCompactView> CTPPSCompositeESSource::produceGeom(const IdealGeometryRecord &){
  return currentProfile_->ddCompactView;
}
//----------------------------------------------------------------------------------------------------
std::unique_ptr<LHCOpticalFunctionsSetCollection> CTPPSCompositeESSource::produceOptics(const CTPPSOpticsRcd &) {
  return std::make_unique<LHCOpticalFunctionsSetCollection>(currentProfile_->lhcOptical);
}
//----------------------------------------------------------------------------------------------------
std::unique_ptr<LHCInfo> CTPPSCompositeESSource::produceLhcInfo(const LHCInfoRcd &) {
  auto lhcInfo = std::make_unique<LHCInfo>();
  //lhcInfo
  if(currentProfile_->xangle<0){
    const double u = CLHEP::RandFlat::shoot(m_engine_.get(), 0., 1.);
    for (const auto &d : currentProfile_->xangleBetaStarBins) {
      if (d.min <= u && u <= d.max) {
        currentProfile_->xangle = d.data.first;
        currentProfile_->betaStar=d.data.second;
        break;
      }
    }
  }

  lhcInfo->setEnergy(currentProfile_->m_beamEnergy);
  lhcInfo->setCrossingAngle(currentProfile_->xangle);
  lhcInfo->setBetaStar(currentProfile_->betaStar);
  edm::LogVerbatim("CTPPSCompositeESSource::produceLhcInfo")<<"new xangle: "<<currentProfile_->xangle<<" betaStar: "<<currentProfile_->betaStar<<std::endl;
  return lhcInfo;
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_EVENTSETUP_SOURCE(CTPPSCompositeESSource);