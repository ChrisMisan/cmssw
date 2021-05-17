// -*- C++ -*-
//
// Package:    SampicDigi/SampicTotemTimingProducer
// Class:      SampicTotemTimingProducer
//
/**\class SampicTotemTimingProducer SampicTotemTimingProducer.cc SampicDigi/SampicTotemTimingProducer/plugins/SampicTotemTimingProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Krzysztof Misan
//         Created:  Sun, 07 Mar 2021 14:42:52 GMT
//
//

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/CTPPSDigi/interface/TotemTimingDigi.h"
#include "DataFormats/CTPPSDetId/interface/TotemTimingDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

const int MAX_SAMPIC_CHANNELS=32;  
const int SAMPIC_SAMPLES=64;

//
// class declaration
//

class SampicTotemTimingProducer : public edm::stream::EDProducer<> {
public:
  explicit SampicTotemTimingProducer(const edm::ParameterSet&);
  ~SampicTotemTimingProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  std::vector<std::string> sampicFilesVec;
  std::unordered_map<unsigned int, std::vector<unsigned int>> detid_vs_chid_;
  TChain *inputTree;
  std::stringstream ssLog;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
using namespace edm;
using namespace std;

SampicTotemTimingProducer::SampicTotemTimingProducer(const edm::ParameterSet& iConfig)
: sampicFilesVec(iConfig.getParameter<std::vector<std::string>>("sampicFilesVec"))
{
  //register your products
/* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
*/

  for (const auto& id_map : iConfig.getParameter<std::vector<edm::ParameterSet>>("idsMapping"))
    detid_vs_chid_[id_map.getParameter<unsigned int>("treeChId")] =id_map.getParameter<vector<unsigned int>>("detId");

  inputTree = new TChain("desy");
	for ( uint i = 0; i < sampicFilesVec.size() ; i++)
		inputTree->Add(sampicFilesVec[i].c_str());


  produces<DetSetVector<TotemTimingDigi>>("TotemTiming");

  //now do what ever other initialization is needed
}

SampicTotemTimingProducer::~SampicTotemTimingProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
  delete inputTree;
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void SampicTotemTimingProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {



 /* edm::FileInPath fip(sampicFilesVec.c_str());
  std::unique_ptr<TFile> f_in(TFile::Open(fip.fullPath().c_str()));
  if (!f_in)
    throw cms::Exception("PPS") << "Cannot open input file '" << sampicFilesVec << "'.";



  TTree* tree = (TTree*)f_in->Get("desy");*/

  int eventNum=iEvent.id().event();
  auto digi = std::make_unique<edm::DetSetVector<TotemTimingDigi>>();

  uint num_samples;
  double trigger_time;
  int sample_channel[MAX_SAMPIC_CHANNELS];
  int sample_cellInfoForOrderedData[MAX_SAMPIC_CHANNELS];
  uint sample_timestampA[MAX_SAMPIC_CHANNELS];
  uint sample_timestampB[MAX_SAMPIC_CHANNELS];
  ULong64_t sample_timestampFPGA[MAX_SAMPIC_CHANNELS];
  double sample_amp[MAX_SAMPIC_CHANNELS][SAMPIC_SAMPLES];

  inputTree->SetBranchAddress("num_samples", &num_samples);
  inputTree->SetBranchAddress("trigger_time", &trigger_time);
  inputTree->SetBranchAddress("sample_channel", sample_channel);
  inputTree->SetBranchAddress("sample_timestampFPGA", sample_timestampFPGA);
  inputTree->SetBranchAddress("sample_timestampA", sample_timestampA);
  inputTree->SetBranchAddress("sample_timestampB", sample_timestampB);
  inputTree->SetBranchAddress("sample_cellInfoForOrderedData", sample_cellInfoForOrderedData);
  inputTree->SetBranchAddress("sample_ampl", sample_amp);

  inputTree->GetEntry(eventNum);
      
  int bunchNumber = ((int) trigger_time/25 ) % 3564;
  int orbitNumber = (int) (trigger_time / 88924.45);


  
  if(num_samples==0){
    //edm::DetSet<TotemTimingDigi>& digis_for_detid = digi->find_or_insert(detid_vs_chid_.at(8)[0]);
    //TotemTimingDigi digiZero;
    //digis_for_detid.push_back(digiZero);
    iEvent.put(std::move(digi),"TotemTiming");
    return;
  }
  TotemTimingEventInfo eventInfoTmp(0,
                                                  sample_timestampFPGA[0],//l1ATimestamp
                                                  bunchNumber,
                                                  orbitNumber,
                                                  eventNum,
                                                  1,//
                                                  0,//l1ALatency
                                                  64,//
                                                  0,//
                                                  1);

  for(uint i=0;i<num_samples;i++){

    unsigned short ch_id=sample_channel[i];//sample_channel[rand()%num_samples];

    std::vector<uint8_t> samples;

    for(int y=0;y<SAMPIC_SAMPLES;y++){
      samples.push_back((int)(sample_amp[i][y]/1.2*256));
      ssLog<<(int)(sample_amp[i][y]/1.2*256)<<endl;}
        

    if(ch_id<8){
      //edm::DetSet<TotemTimingDigi>& digis_for_detid = digi->find_or_insert(detid_vs_chid_.at(8)[0]);
      //TotemTimingDigi digiZero;
      //digiZero.setEventInfo(eventInfoTmp);
      //digis_for_detid.push_back(digiZero);
      continue;
    }

    TotemTimingDigi digiTmp(0,
                                        sample_timestampFPGA[i],
                                        sample_timestampA[i],
                                        sample_timestampB[i],
                                        sample_cellInfoForOrderedData[i],
                                        samples,
                                        eventInfoTmp);
        

    auto vec=detid_vs_chid_.at(ch_id);
    for(auto id:vec){
      CTPPSDiamondDetId detId(id);
      edm::DetSet<TotemTimingDigi>& digis_for_detid = digi->find_or_insert(detId);
      digis_for_detid.push_back(digiTmp);
    }

  }
  //ssLog<<"EventNum3: "<<eventNum<<" "<<digi->begin()->begin()->timestampA()<<std::endl;
  iEvent.put(std::move(digi),"TotemTiming");
  
  if(eventNum==10000)
    edm::LogInfo("SampicTotemTimingProducer") << ssLog.str();
/*
  uint num_samples;
  double trigger_time;
  //event scope
  tree->SetBranchAddress("num_samples", &num_samples);
  tree->SetBranchAddress("trigger_time", &trigger_time);
  tree->GetEntry(eventNum);

  
  if(num_samples>0){

    //waveform/channel scope
    int* sample_channel=new int[num_samples];
    int* sample_cellInfoForOrderedData=new int[num_samples];
    uint* sample_timestampA=new uint[num_samples];
    uint* sample_timestampB=new uint[num_samples];
    ULong64_t* sample_timestampFPGA=new ULong64_t[num_samples];
    tree->SetBranchAddress("sample_channel", sample_channel);
    tree->SetBranchAddress("sample_timestampFPGA", sample_timestampFPGA);
    tree->SetBranchAddress("sample_timestampA", sample_timestampA);
    tree->SetBranchAddress("sample_timestampB", sample_timestampB);
    tree->SetBranchAddress("sample_cellInfoForOrderedData", sample_cellInfoForOrderedData);

    //sample scope
    auto sample_amp=new double[num_samples][64];
    tree->SetBranchAddress("sample_ampl", sample_amp);

    tree->GetEntry(eventNum);
  
    int bunchNumber = ((int) trigger_time/25 ) % 3564;
    int orbitNumber = (int) (trigger_time / 88924.45);

    

    for(uint i=0;i<num_samples;i++){
      std::vector<uint8_t>* samples=new std::vector<uint8_t>();
      for(int y=0;y<64;y++)
        samples->push_back((int)(sample_amp[i][y]/1.2)*256);
      

      TotemTimingEventInfo eventInfoTmp(sample_channel[i],
                                                trigger_time,
                                                bunchNumber,
                                                orbitNumber,
                                                eventNum,
                                                1,//
                                                1,//
                                                64,//
                                                0,//
                                                1);
      
      TotemTimingDigi digiTmp(sample_channel[i],
                                      sample_timestampFPGA[i],
                                      sample_timestampA[i],
                                      sample_timestampB[i],
                                      sample_cellInfoForOrderedData[i],
                                      *samples,
                                      eventInfoTmp);
      
      

      TotemTimingDetId detId(detid_vs_chid_.at(sample_channel[i]));

      DetSet<TotemTimingDigi> &digiDetSet = digi.find_or_insert(detId);
      digiDetSet.push_back(digiTmp);
    
    }
  }

  iEvent.put(make_unique<DetSetVector<TotemTimingDigi>>(digi),"TotemTiming");
  f_in->Close();


 This is an event example
  //Read 'ExampleData' from the Event
  ExampleData const& in = iEvent.get(inToken_);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  iEvent.put(std::make_unique<ExampleData2>(in));


 this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
*/
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void SampicTotemTimingProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void SampicTotemTimingProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
SampicTotemTimingProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
SampicTotemTimingProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SampicTotemTimingProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SampicTotemTimingProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SampicTotemTimingProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters

  edm::ParameterSetDescription desc;
  desc.add<vector<std::string>>("sampicFilesVec")->setComment("path to sampic root data");

  edm::ParameterSetDescription idmap_valid;
  idmap_valid.add<unsigned int>("treeChId", 0)->setComment("HW id as retrieved from tree");
  idmap_valid.add<vector<unsigned int>>("detId")->setComment("mapped CTPPSDiamondDetId's for this channel");

  desc.addVPSet("idsMapping",idmap_valid);

  descriptions.add("SampicTotemTimingProducer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SampicTotemTimingProducer);
