#include "SusyEvent.h"
#include "../CommonCode/ObjectVars.h"
#include "../CommonCode/ObjectSelector.h"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <set>
#include <cmath>

class UnclusteredMETCorr {
public:
  UnclusteredMETCorr();
  ~UnclusteredMETCorr();
  bool initialize(char const*);
  void addInput(char const*);
  bool run();
  void clearInput();
  bool finalize();

private:
  TChain input_;

  TTree* output_;
};

UnclusteredMETCorr::UnclusteredMETCorr() :
  input_("susyTree"),
  output_(0)
{
}

UnclusteredMETCorr::~UnclusteredMETCorr()
{
  delete output_;
}

bool
UnclusteredMETCorr::initialize(char const* _outputName)
{
  TString outputName(_outputName);
  if(!outputName.Contains(".root")) outputName += "/metCorr.root";

  TFile* outputFile(TFile::Open(outputName, "recreate"));

  outputFile->cd();
  output_ = new TTree("metCorr", "MET correction");

  output_->Branch("metx", 0, "metx/F");
  output_->Branch("mety", 0, "mety/F");
  output_->Branch("unclx", 0, "unclx/F");
  output_->Branch("uncly", 0, "uncly/F");

  return true;
}

void
UnclusteredMETCorr::addInput(char const* _source)
{
  input_.Add(_source);
}

bool
UnclusteredMETCorr::run()
{
  input_.SetBranchStatus("*", 0);
  input_.SetBranchStatus("met_pfType1CorrectedMet*", 1);
  input_.SetBranchStatus("pfParticles*", 1);
  input_.SetBranchStatus("pfJets_ak5*", 1);
  susy::PhotonVars::setBranchStatus(input_);
  susy::ElectronVars::setBranchStatus(input_);
  susy::MuonVars::setBranchStatus(input_);

  float metx;
  float mety;
  float unclx;
  float uncly;

  output_->SetBranchAddress("metx", &metx);
  output_->SetBranchAddress("mety", &mety);
  output_->SetBranchAddress("unclx", &unclx);
  output_->SetBranchAddress("uncly", &uncly);

  susy::Event event;
  event.setInput(input_);

  long iEntry(0);
  while(event.getEntry(iEntry++)){
    if(iEntry % 10000 == 1) std::cout << iEntry << std::endl;

    susy::PFParticleCollection const& pfs(event.pfParticles);
    unsigned nPF(pfs.size());    

    std::set<const susy::PFParticle*> clustered;

    bool hasLepton(false);

    susy::ElectronCollection const& electrons(event.electrons["gsfElectrons"]);
    unsigned nE(electrons.size());
    for(unsigned iE(0); iE != nE; ++iE){
      susy::ElectronVars vars(electrons[iE], event);
      if(!vars.isVeto) continue;
      
      for(unsigned iPF(0); iPF != nPF; ++iPF){
        if(std::abs(pfs[iPF].pdgId) != 11) continue;
        if(pfs[iPF].momentum.DeltaR(electrons[iE].momentum) < 0.1){
          clustered.insert(&pfs[iPF]);
          break;
        }
      }

      if(vars.pt > 25. && vars.isMedium) hasLepton = true;
    }

    susy::MuonCollection const& muons(event.muons["muons"]);
    unsigned nM(muons.size());
    for(unsigned iM(0); iM != nM; ++iM){
      susy::MuonVars vars(muons[iM], event);
      if(!vars.isLoose) continue;
      
      for(unsigned iPF(0); iPF != nPF; ++iPF){
        if(std::abs(pfs[iPF].pdgId) != 13) continue;
        if(pfs[iPF].momentum.DeltaR(muons[iM].momentum) < 0.1){
          clustered.insert(&pfs[iPF]);
          break;
        }
      }

      if(vars.pt > 25. && vars.isTight) hasLepton = true;
    }

    if(!hasLepton) continue;

    bool hasPhoton(false);

    susy::PhotonCollection const& photons(event.photons["photons"]);
    unsigned nP(photons.size());
    for(unsigned iP(0); iP != nP; ++iP){
      susy::PhotonVars vars(photons[iP], event);
      if(!vars.isLoose) continue;

      for(unsigned iPF(0); iPF != nPF; ++iPF){
        if(pfs[iPF].pdgId != 22) continue;
        if(pfs[iPF].momentum.DeltaR(photons[iP].momentum) < 0.1){
          clustered.insert(&pfs[iPF]);
          break;
        }
      }

      if(vars.isLoose && vars.nPixelSeeds == 0) hasPhoton = true;
    }

    if(!hasPhoton) continue;

    susy::PFJetCollection const& jets(event.pfJets["ak5"]);
    unsigned nJ(jets.size());
    
    for(unsigned iJ(0); iJ != nJ; ++iJ)
      clustered.insert(jets[iJ].pfParticles.begin(), jets[iJ].pfParticles.end());

    std::set<const susy::PFParticle*>::const_iterator cEnd(clustered.end());

    unclx = 0.;
    uncly = 0.;

    for(unsigned iPF(0); iPF != nPF; ++iPF){
      if(clustered.find(&pfs[iPF]) == cEnd){
        unclx += pfs[iPF].momentum.X();
        uncly += pfs[iPF].momentum.Y();
      }
    }

    susy::MET const& met(event.metMap["pfType1CorrectedMet"]);
    metx = met.mEt.X();
    mety = met.mEt.Y();

    output_->Fill();
  }

  event.releaseTrees();

  return true;
}

void
UnclusteredMETCorr::clearInput()
{
  input_.Reset();
}

bool
UnclusteredMETCorr::finalize()
{
  TFile* outputFile(output_->GetCurrentFile());
  outputFile->cd();
  output_->Write();
  delete outputFile;

  output_ = 0;

  return true;
}

void
unclustered(TString const& _sourceName, TString const& _outputName)
{
  UnclusteredMETCorr corr;
  if(!corr.initialize(_outputName)) return;
  corr.addInput(_sourceName);
  if(!corr.run()) return;
  corr.clearInput();
  corr.finalize();
}
