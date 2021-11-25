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

#include "lardataobj/RecoBase/Hit.h"

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

  std::string fHitProducer; // Module label that created the hits
  float fTPCReadoutWindowSize; // Number of samples in TPC readout window
  std::vector<float> fTPCWireNumbers; // Number of wires in TPC
  float fTPCDownsampleTicks; // Downsampling factor for TPC ticks
  float fTPCDownsampleWires; // Downsampling factor for TPC wires

  TH2F* fhDisplayHitU; // Plane U hit position
  TH2F* fhDisplayHitV; // Plane V hit position
  TH2F* fhDisplayHitY; // Plane Y hit position

  TTree* fOutTree; // Output tree

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

};


NNDumper::NNDumper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fHitProducer = p.get<std::string>("HitProducer");
  fTPCReadoutWindowSize = p.get<float>("TPCReadoutWindowSize");
  fTPCWireNumbers = p.get< std::vector<float> >("TPCWireNumbers");
  fTPCDownsampleTicks = p.get<float>("TPCDownsampleTicks");
  fTPCDownsampleWires = p.get<float>("TPCDownsampleWires");
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void NNDumper::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  //*** Process TPC information ***

  art::InputTag hit_tag { fHitProducer };
  auto const& hit_handle = e.getValidHandle< std::vector<recob::Hit> >(hit_tag);

  for( auto const& hit : *hit_handle ){
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

  fOutTree->Fill();

  resetImage(fhDisplayHitU, fTPCImgDataU);
  resetImage(fhDisplayHitV, fTPCImgDataV);
  resetImage(fhDisplayHitY, fTPCImgDataY);

}

void NNDumper::beginJob()
{
  // Implementation of optional member function here.
  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us.
  art::ServiceHandle<art::TFileService> tfs;

  // Use the ROOT TH2 to create and downsample the TPC "image"

  // To do: store 3 plane information as elements of vector
  fhDisplayHitU = new TH2F("hDisplayHitU","Plane U Evt. hit pos.; Wire; Tick; Entries", 
			   int(fTPCWireNumbers[0]/fTPCDownsampleWires), 0., fTPCWireNumbers[0], 
			   int(2*fTPCReadoutWindowSize/fTPCDownsampleTicks), -fTPCReadoutWindowSize, fTPCReadoutWindowSize);
  fhDisplayHitV = new TH2F("hDisplayHitV","Plane V Evt. hit pos.; Wire; Tick; Entries",
			   int(fTPCWireNumbers[1]/fTPCDownsampleWires), 0., fTPCWireNumbers[1], 
			   int(2*fTPCReadoutWindowSize/fTPCDownsampleTicks), -fTPCReadoutWindowSize, fTPCReadoutWindowSize);
  fhDisplayHitY = new TH2F("hDisplayHitY","Plane Y Evt. hit pos.; Wire; Tick; Entries",
			   int(fTPCWireNumbers[2]/fTPCDownsampleWires), 0., fTPCWireNumbers[2], 
			   int(2*fTPCReadoutWindowSize/fTPCDownsampleTicks), -fTPCReadoutWindowSize, fTPCReadoutWindowSize);

  // Output tree
  fOutTree = tfs->make<TTree>("evttree","NNDumper output tree");
  // To do: store 3 plane information as elements of vector after checking we can unpack them later
  fOutTree->Branch("TPCImgDataU", "std::vector<float>", &fTPCImgDataU);
  fOutTree->Branch("TPCImgWidthU", &fTPCImgWidthU, "TPCImgWidthU/I");
  fOutTree->Branch("TPCImgHeightU", &fTPCImgHeightU, "TPCImgHeightU/I");
  fOutTree->Branch("TPCImgDataV", "std::vector<float>", &fTPCImgDataV);
  fOutTree->Branch("TPCImgWidthV", &fTPCImgWidthV, "TPCImgWidthV/I");
  fOutTree->Branch("TPCImgHeightV", &fTPCImgHeightV, "TPCImgHeightV/I");
  fOutTree->Branch("TPCImgDataY", "std::vector<float>", &fTPCImgDataY);
  fOutTree->Branch("TPCImgWidthY", &fTPCImgWidthY, "TPCImgWidthY/I");
  fOutTree->Branch("TPCImgHeightY", &fTPCImgHeightY, "TPCImgHeightY/I");

  fTPCImgWidthU = fhDisplayHitU->GetNbinsX();
  fTPCImgHeightU = fhDisplayHitU->GetNbinsY();
  fTPCImgDataU.resize(fTPCImgWidthU*fTPCImgHeightU, 0); // fill with zeroes

  fTPCImgWidthV = fhDisplayHitV->GetNbinsX();
  fTPCImgHeightV = fhDisplayHitV->GetNbinsY();
  fTPCImgDataV.resize(fTPCImgWidthV*fTPCImgHeightV, 0); // fill with zeroes

  fTPCImgWidthY = fhDisplayHitY->GetNbinsX();
  fTPCImgHeightY = fhDisplayHitY->GetNbinsY();
  fTPCImgDataY.resize(fTPCImgWidthY*fTPCImgHeightY, 0); // fill with zeroes

}

void NNDumper::endJob()
{
  // Implementation of optional member function here.
  delete fhDisplayHitU;
  delete fhDisplayHitV;
  delete fhDisplayHitY;
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
  //h2->Write(); // write to output file for debugging
  h2->Reset();
  std::fill(img.begin(), img.end(), 0);;
}

DEFINE_ART_MODULE(NNDumper)
