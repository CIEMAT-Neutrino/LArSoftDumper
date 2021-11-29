////////////////////////////////////////////////////////////////////////
// Class:       NNDumper
// Plugin Type: analyzer (Unknown Unknown)
// File:        NNDumper_module.cc
//
// Generated at Wed Nov 24 08:30:17 2021 by Jose Ignacio Crespo Anadon using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "larcore/Geometry/Geometry.h"

//Charge Hits
#include "lardataobj/RecoBase/Hit.h"

//Optical Hits
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "lardataobj/RecoBase/OpHit.h"

//MC Truth
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "TTimeStamp.h"

//Sim::Photons (truth photons arriving at Optical channels from G4)
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "TH2F.h"
#include "TTree.h"


class NNDumper;


class NNDumper : public art::EDAnalyzer {
public:
  explicit NNDumper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NNDumper(NNDumper const&) = delete;
  NNDumper(NNDumper&&) = delete;
  NNDumper& operator=(NNDumper const&) = delete;
  NNDumper& operator=(NNDumper&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void extractImage(TH2F* h2, std::vector<float> &img );
  void resetImage(TH2F* h2, std::vector<float> &img );

private:

  // Declare member data here.

  std::string fHitProducer; // Module label that created the optical hits
  std::string fOpHitProducer; // Module label that created the hits
  std::string fMCTruthProducer; // Module label that created the hits

  TTree* fOutTree; // Output tree

  bool fDebug; // Debug flag

  //----Charge Hits
  float fTPCReadoutWindowSize; // Number of samples in TPC readout window
  std::vector<float> fTPCWireNumbers; // Number of wires in TPC
  float fTPCDownsampleTicks; // Downsampling factor for TPC ticks
  float fTPCDownsampleWires; // Downsampling factor for TPC wires
  TH2F* fhDisplayHitU; // Plane U hit position
  TH2F* fhDisplayHitV; // Plane V hit position
  TH2F* fhDisplayHitY; // Plane Y hit position

  // To do: store 3 plane information as elements of vector after checking we can unpack them later
  std::vector<float> fTPCImgDataU; // TPC image data
  uint fTPCImgWidthU; // TPC image width
  uint fTPCImgHeightU; // TPC image height
  std::vector<float> fTPCImgDataV; // TPC image data
  uint fTPCImgWidthV; // TPC image width
  uint fTPCImgHeightV; // TPC image height
  std::vector<float> fTPCImgDataY; // TPC image data
  uint fTPCImgWidthY; // TPC image width
  uint fTPCImgHeightY; // TPC image height

  //----OpHits
  float fOpReadoutWindowLength; // Length of PDS readout window in ns
  float fOpReadoutSample; // Sample size in ns
  float fOpDownsampleTicks; // Downsample factor for optical ticks
  float fOpChannelBins; // Number of optical channels
  TH2F* fhDisplayOpHit0; // Optical hit position in TPC 0
  TH2F* fhDisplayOpHit1; // Optical hit position in TPC 1
  std::vector<float> fOpImgData0; // Optical image data in TPC 0
  std::vector<float> fOpImgData1; // Optical image data in TPC 1
  uint fOpImgWidth; //  Optical image width
  uint fOpImgHeight; // Optical image height
  opdet::sbndPDMapAlg pdsMap; // Optical channels map
  // To do: get from fcl file
  std::vector<int> fOpChannels = pdsMap.getChannelsOfType("pmt_coated"); // List of selected optical channels

  //----MCTruth (interaction info)
  int fInteractionType;

  //----SimPhotons
  int NChannels = fOpChannels.size();
  std::vector<int> fSimPhotons; //vector of zeros for the channels
  // JICA: why not remove the two lines above and use the simpler version
  //std::vector<int> fSimPhotons(fOpChannels.size(), 0);
};


NNDumper::NNDumper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fHitProducer = p.get<std::string>("HitProducer");
  fOpHitProducer = p.get<std::string>("OpHitProducer");
  fMCTruthProducer = p.get<std::string>("MCTruthProducer");
  fTPCReadoutWindowSize = p.get<float>("TPCReadoutWindowSize");
  fTPCWireNumbers = p.get< std::vector<float> >("TPCWireNumbers");
  fTPCDownsampleTicks = p.get<float>("TPCDownsampleTicks");
  fTPCDownsampleWires = p.get<float>("TPCDownsampleWires");
  fDebug = p.get<bool>("Debug",false);
  fOpReadoutWindowLength = p.get<float>("OpReadoutWindowLength");
  fOpReadoutSample = p.get<float>("OpReadoutSample");
  fOpDownsampleTicks = p.get<float>("OpDownsampleTicks");
  fOpChannelBins = p.get<float>("OpChannelBins");
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void NNDumper::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  //*** Process TPC information ***

  art::InputTag hit_tag { fHitProducer };
  auto const& hit_handle = e.getValidHandle< std::vector<recob::Hit> >(hit_tag);

  for( auto const& hit : *hit_handle ) {
    unsigned int wire = hit.WireID().Wire;
    float tick = hit.PeakTime();
    float integral = hit.Integral();

    // Stitch the two TPC drift volumes together
    if( hit.WireID().TPC == 0 ) tick = tick - fTPCReadoutWindowSize;
    else if( hit.WireID().TPC == 1 ) tick = -(tick - fTPCReadoutWindowSize);

    if( hit.View() == 0 ) fhDisplayHitU->Fill( wire, tick, integral );
    else if( hit.View() == 1 ) fhDisplayHitV->Fill( wire, tick, integral );
    else if( hit.View() == 2 ) fhDisplayHitY->Fill( wire, tick, integral );
  }

  extractImage(fhDisplayHitU, fTPCImgDataU);
  extractImage(fhDisplayHitV, fTPCImgDataV);
  extractImage(fhDisplayHitY, fTPCImgDataY);

  //*** Process Optical information ***
  art::Handle<std::vector<recob::OpHit>> ophits_h;
  e.getByLabel(fOpHitProducer.c_str(), ophits_h);
  for(auto const& oph : *ophits_h){
    auto ch            = oph.OpChannel();
    auto peak_time_abs = oph.PeakTimeAbs();
    // peak_time     = oph.PeakTime();
    // width         = oph.Width();
    // area          = oph.Area();
    // amplitude     = oph.Amplitude();
    auto PE            = oph.PE();

    if(pdsMap.pdType(ch) == "pmt_coated"){
      // Get index of vector
      std::vector<int>::iterator it = std::find(fOpChannels.begin(), fOpChannels.end(), ch);
      int index = std::distance(fOpChannels.begin(), it);
      // SBND specific: even channels belong to TPC 0, odd channels to TPC 1
      if( ch%2 == 0 ) fhDisplayOpHit0->Fill( index/2, peak_time_abs, PE );
      else fhDisplayOpHit1->Fill( index/2, peak_time_abs, PE );
    }
  }

  extractImage(fhDisplayOpHit0, fOpImgData0);
  extractImage(fhDisplayOpHit1, fOpImgData1);

  //*** Process Interaction(MCTruth) information ***
  art::Handle<std::vector<simb::MCTruth> > mctruths;
  e.getByLabel(fMCTruthProducer.c_str(), mctruths);
  for (auto const& truth : *mctruths){
    const simb::MCNeutrino& nu = truth.GetNeutrino();
    fInteractionType = nu.InteractionType();
  }

  //*** Process Truth photons on optical channels (SimPhotons) information ***
  auto const& photon_handles = e.getMany<std::vector<sim::SimPhotonsLite>>();
  for (const art::Handle<std::vector<sim::SimPhotonsLite>> &opdetHandle : photon_handles) {
    // this now tells you if light collection is reflected
    // const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");

    for (auto const& litesimphotons : (*opdetHandle)){
      const unsigned ch = litesimphotons.OpChannel;
      const std::string pdtype = pdsMap.pdType(ch);
      std::map<int, int> const& photonMap = litesimphotons.DetectedPhotons;
    
      if(pdtype == "pmt_coated"){

        for (auto const& photonMember : photonMap){
          auto NPhotons = photonMember.second;
          // auto tphoton = photonMember.first;//maybe in future work!
          std::vector<int>::iterator it = std::find(fOpChannels.begin(), fOpChannels.end(), ch);
          int index = std::distance(fOpChannels.begin(), it);
          fSimPhotons[index] += NPhotons;
        }
      }
    }
  }
  fOutTree->Fill();

  resetImage(fhDisplayHitU, fTPCImgDataU);
  resetImage(fhDisplayHitV, fTPCImgDataV);
  resetImage(fhDisplayHitY, fTPCImgDataY);
  resetImage(fhDisplayOpHit0, fOpImgData0);
  resetImage(fhDisplayOpHit1, fOpImgData1);
  std::fill(fSimPhotons.begin(), fSimPhotons.end(), 0);//reset simphotons vector
}

void NNDumper::beginJob()
{
  // Implementation of optional member function here.
  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us.
  art::ServiceHandle<art::TFileService> tfs;
  fSimPhotons = std::vector<int>(NChannels);
  // Use the ROOT TH2 to create and downsample the "images"

  // To do: store 3 plane information as elements of vector
  fhDisplayHitU = new TH2F("hDisplayHitU","Plane U Evt. hit pos.; Wire; Tick; Charge",
			   int(fTPCWireNumbers[0]/fTPCDownsampleWires), 0., fTPCWireNumbers[0],
			   int(2*fTPCReadoutWindowSize/fTPCDownsampleTicks), -fTPCReadoutWindowSize, fTPCReadoutWindowSize);
  fhDisplayHitV = new TH2F("hDisplayHitV","Plane V Evt. hit pos.; Wire; Tick; Charge",
			   int(fTPCWireNumbers[1]/fTPCDownsampleWires), 0., fTPCWireNumbers[1],
			   int(2*fTPCReadoutWindowSize/fTPCDownsampleTicks), -fTPCReadoutWindowSize, fTPCReadoutWindowSize);
  fhDisplayHitY = new TH2F("hDisplayHitY","Plane Y Evt. hit pos.; Wire; Tick; Charge",
			   int(fTPCWireNumbers[2]/fTPCDownsampleWires), 0., fTPCWireNumbers[2],
			   int(2*fTPCReadoutWindowSize/fTPCDownsampleTicks), -fTPCReadoutWindowSize, fTPCReadoutWindowSize);

  // To do: store 2 optical plane information as elements of vector
  fhDisplayOpHit0 = new TH2F("hDisplayOpHit0","TPC 0 Evt. OpHit pos.; Channel; Time (#mus); PEs",
			    fOpChannelBins, 0., fOpChannelBins,
			    int(fOpReadoutWindowLength/(fOpReadoutSample*fOpDownsampleTicks)), 0., fOpReadoutWindowLength/1000.); // in mus
  fhDisplayOpHit1 = new TH2F("hDisplayOpHit1","TPC 1 Evt. OpHit pos.; Channel; Time (#mus); PEs",
			    fOpChannelBins, 0., fOpChannelBins,
			    int(fOpReadoutWindowLength/(fOpReadoutSample*fOpDownsampleTicks)), 0., fOpReadoutWindowLength/1000.); // in mus

  // Output tree
  fOutTree = tfs->make<TTree>("evttree","NNDumper output tree");
  // To do: store 3 TPC plane information as elements of vector after checking we can unpack them later
  fOutTree->Branch("TPCImgDataU", "std::vector<float>", &fTPCImgDataU);
  fOutTree->Branch("TPCImgWidthU", &fTPCImgWidthU, "TPCImgWidthU/I");
  fOutTree->Branch("TPCImgHeightU", &fTPCImgHeightU, "TPCImgHeightU/I");
  fOutTree->Branch("TPCImgDataV", "std::vector<float>", &fTPCImgDataV);
  fOutTree->Branch("TPCImgWidthV", &fTPCImgWidthV, "TPCImgWidthV/I");
  fOutTree->Branch("TPCImgHeightV", &fTPCImgHeightV, "TPCImgHeightV/I");
  fOutTree->Branch("TPCImgDataY", "std::vector<float>", &fTPCImgDataY);
  fOutTree->Branch("TPCImgWidthY", &fTPCImgWidthY, "TPCImgWidthY/I");
  fOutTree->Branch("TPCImgHeightY", &fTPCImgHeightY, "TPCImgHeightY/I");
  // To do: store 2 optical plane information as elements of vector after checking we can unpack them later
  fOutTree->Branch("OpImgData0", "std::vector<float>", &fOpImgData0);
  fOutTree->Branch("OpImgData1", "std::vector<float>", &fOpImgData1);
  fOutTree->Branch("OpImgWidth", &fOpImgWidth, "OpImgWidth/I");
  fOutTree->Branch("OpImgHeight", &fOpImgHeight, "OpImgHeight/I");

  fOutTree->Branch("InteractionType", &fInteractionType, "InteractionType/I");

  fOutTree->Branch("SimPhotonsData", "std::vector<int>", &fSimPhotons);

  fTPCImgWidthU = fhDisplayHitU->GetNbinsX();
  fTPCImgHeightU = fhDisplayHitU->GetNbinsY();
  fTPCImgDataU.resize(fTPCImgWidthU*fTPCImgHeightU, 0); // fill with zeroes

  fTPCImgWidthV = fhDisplayHitV->GetNbinsX();
  fTPCImgHeightV = fhDisplayHitV->GetNbinsY();
  fTPCImgDataV.resize(fTPCImgWidthV*fTPCImgHeightV, 0); // fill with zeroes

  fTPCImgWidthY = fhDisplayHitY->GetNbinsX();
  fTPCImgHeightY = fhDisplayHitY->GetNbinsY();
  fTPCImgDataY.resize(fTPCImgWidthY*fTPCImgHeightY, 0); // fill with zeroes

  fOpImgWidth = fhDisplayOpHit0->GetNbinsX();
  fOpImgHeight = fhDisplayOpHit0->GetNbinsY();
  fOpImgData0.resize(fOpImgWidth*fOpImgHeight, 0); // fill with zeroes
  // Must have the same size as above
  fOpImgData1.resize(fOpImgWidth*fOpImgHeight, 0); // fill with zeroes
}

void NNDumper::endJob()
{
  // Implementation of optional member function here.
  delete fhDisplayHitU;
  delete fhDisplayHitV;
  delete fhDisplayHitY;
  delete fhDisplayOpHit0;
  delete fhDisplayOpHit1;
}

void NNDumper::extractImage(TH2F* h2, std::vector<float>& img)
{
  size_t index = 0;
  for(int biny = 1; biny < h2->GetNbinsY() + 1; biny++){
    for(int binx = 1; binx < h2->GetNbinsX() + 1; binx++){
      img.at(index) = h2->GetBinContent(binx, biny);
      index++;
    }
  }
}

void NNDumper::resetImage(TH2F* h2, std::vector<float>& img)
{
  if(fDebug) h2->Write(); // write to output file for debugging
  h2->Reset();
  std::fill(img.begin(), img.end(), 0);;
}

DEFINE_ART_MODULE(NNDumper)
