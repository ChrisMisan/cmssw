// -*- C++ -*-
//
// Package:    SampicDigi/WaveformAnalyzer
// Class:      WaveformAnalyzer
//
/**\class WaveformAnalyzer WaveformAnalyzer.cc SampicDigi/WaveformAnalyzer/plugins/WaveformAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Krzysztof Misan
//         Created:  Mon, 28 Jun 2021 16:08:10 GMT
//
//
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

// system include files
#include <memory>

// user include files
#include <map>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CTPPSDigi/interface/TotemTimingDigi.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "DataFormats/CTPPSReco/interface/TotemTimingRecHit.h"

#include "RecoPPS/Local/interface/TotemTimingConversions.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "CondFormats/DataRecord/interface/PPSTimingCalibrationRcd.h"

#include "TH1.h"
#include "TGraph.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class WaveformAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit WaveformAnalyzer(const edm::ParameterSet&);
  ~WaveformAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void initHistograms( const CTPPSDiamondDetId&);
  // ----------member data ---------------------------
  edm::EDGetTokenT< edm::DetSetVector<TotemTimingDigi> > tokenDigi_;
  edm::ESGetToken<PPSTimingCalibration, PPSTimingCalibrationRcd> timingCalibrationToken_;


  std::map< CTPPSDiamondDetId, TFileDirectory > maindir_map_;
  //TotemTimingConversions* conv;
  TotemTimingConversions* conv2;

  static const double SAMPIC_MAX_NUMBER_OF_SAMPLES;
  static const double SAMPIC_ADC_V;
  static const double SAMPIC_SAMPLING_PERIOD_NS;
  int ct;
};

//
// constants, enums and typedefs
//
const double    WaveformAnalyzer::SAMPIC_MAX_NUMBER_OF_SAMPLES = 64;
const double    WaveformAnalyzer::SAMPIC_ADC_V = 1./256;
const double    WaveformAnalyzer::SAMPIC_SAMPLING_PERIOD_NS = 1./7.8;
//
// static data member definitions
//

//
// constructors and destructor
//
WaveformAnalyzer::WaveformAnalyzer(const edm::ParameterSet& iConfig)
    : tokenDigi_( consumes< edm::DetSetVector<TotemTimingDigi> >( iConfig.getParameter<edm::InputTag>( "tagDigi" ) )),
      timingCalibrationToken_(esConsumes<PPSTimingCalibration, PPSTimingCalibrationRcd>(
          edm::ESInputTag(iConfig.getParameter<std::string>("timingCalibrationTag")))){
  
  usesResource("TFileService");
  //now do what ever initialization is needed
}

WaveformAnalyzer::~WaveformAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

void WaveformAnalyzer::initHistograms(const CTPPSDiamondDetId& detId)
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  if (maindir_map_.find(detId) == maindir_map_.end())
  {
    std::string dirName;
    detId.rpName(dirName, CTPPSDiamondDetId::nPath);
    std::string chName;
    detId.channelName(chName, CTPPSDiamondDetId::nFull);

    // create directory for the detector, if not already done
    maindir_map_[detId] = fs->mkdir( dirName);

  }

}

// ------------ method called for each event  ------------
void WaveformAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  edm::ESHandle<PPSTimingCalibration> hTimingCalib = iSetup.getHandle(timingCalibrationToken_);
  conv2 =new TotemTimingConversions(1. / 7.695, false, *hTimingCalib);
  edm::Handle< edm::DetSetVector<TotemTimingDigi> > timingDigi;
  iEvent.getByToken( tokenDigi_, timingDigi );
  for (const auto& digis : *timingDigi){
    const CTPPSDiamondDetId detId( digis.detId() );

    if (maindir_map_.find(detId) == maindir_map_.end())
      initHistograms( detId );
    std::string chName;
    detId.channelName(chName, CTPPSDiamondDetId::nFull);
    // std::string samples_graph_name(chName);
    // samples_graph_name.insert(0, "_cal_");
    // samples_graph_name.insert(0, std::to_string(ct));
    // samples_graph_name.insert(0, "samples_graph_");
    // TGraph *graph_hndl = maindir_map_[ detId ].make<TGraph>();
    // graph_hndl->SetName(samples_graph_name.c_str());
    // graph_hndl->SetTitle(samples_graph_name.c_str());
    // std::vector<float> voltage;
    
    std::string samples_graph_name2(chName);
    samples_graph_name2.insert(0, "_noncal_");
    samples_graph_name2.insert(0, std::to_string(ct));
    samples_graph_name2.insert(0, "samples_graph_");
    TGraph *graph_hndl2 = maindir_map_[ detId ].make<TGraph>();
    graph_hndl2->SetName(samples_graph_name2.c_str());
    graph_hndl2->SetTitle(samples_graph_name2.c_str());
    std::vector<float> voltage2;
    ct++;


    for (const auto& digi : digis )
    {
      //detId_calFileId_map_[detId] = calFileId(digi.getHardwareBoardId(), digi.getHardwareSampicId(), digi.getHardwareChannelId());
      // // Do stuff on DIGIs
      //voltage = conv->getVoltSamples(digi);
      voltage2 = conv2->voltSamples(digi);
      int i=0;
      // for (auto sampIt = voltage.begin(); sampIt != voltage.end(); ++sampIt){
      //   // std::cout << *sampIt << ", ";
      //   graph_hndl->SetPoint(i, SAMPIC_SAMPLING_PERIOD_NS*(i), *sampIt  );
      //   i++;
      // }
      // // std::cout << "\n";
      
      // i=0;
      for (auto sampIt = voltage2.begin(); sampIt != voltage2.end(); ++sampIt){
        graph_hndl2->SetPoint(i, SAMPIC_SAMPLING_PERIOD_NS*(i), *sampIt  );
        i++;
      }
    }

  }


}

// ------------ method called once each job just before starting event loop  ------------
void WaveformAnalyzer::beginJob() {
  //conv = new TotemTimingConversions(false);//accepts path to the calibration file
  ct=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void WaveformAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void WaveformAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WaveformAnalyzer);
