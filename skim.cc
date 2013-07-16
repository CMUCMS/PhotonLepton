#include "TObjArray.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include "SusyEvent.h"
#include "SusyTriggerEvent.h"
#include "../Common/Utilities.h"
#include "../Common/ObjectTree.h"

#include <iostream>
#include <vector>

void
skim(TString const& _skimTypes, TString const& _triggers, TString const& _goodLumis, TObjArray* _urls, TObjArray* _outputNames)
{
  bool skimPhoton(false);
  bool skimElectron(false);
  bool skimMuon(false);
  bool skimGenLepton(false);

  TFile* outputFile(TFile::Open(_outputNames->At(0)->GetName(), "recreate"));
  if(!outputFile || outputFile->IsZombie()){
    std::cerr << "Error: Failed I/O" << std::endl;
    delete outputFile;
    return;
  }

  TObjArray* types(_skimTypes.Tokenize(","));
  for(int iT(0); iT != types->GetEntries(); ++iT){
    TString type(types->At(iT)->GetName());
    if(type == "Photon") skimPhoton = true;
    else if(type == "Electron") skimElectron = true;
    else if(type == "Muon") skimMuon = true;
    else if(type == "GenLepton") skimGenLepton = true;
  }
  delete types;

  TObjArray* trigArr(_triggers.Tokenize(","));
  trigArr->SetOwner(true);
  std::vector<TString> triggerNames;
  for(int iT(0); iT != trigArr->GetEntries(); ++iT)
    triggerNames.push_back(TString(trigArr->At(iT)->GetName()) + "_v*");
  delete trigArr;

  int nTrig(triggerNames.size());

  susy::Event* event(new susy::Event);

  TChain input("susyTree");
  input.AddFileInfoList(_urls);
  input.SetBranchStatus("*", 0);
  input.SetBranchStatus("isRealData", 1);
  input.SetBranchStatus("runNumber", 1);
  input.SetBranchStatus("luminosityBlockNumber", 1);
  input.SetBranchStatus("hlt*", 1);
  if(skimGenLepton) input.SetBranchStatus("genParticles*", 1);
  susy::ObjectTree::setBranchStatus(input, skimPhoton, skimElectron, skimMuon, false);

  event->setInput(input);

  susy::Event* fullEvent(new susy::Event);

  TChain fullInput("susyTree");
  fullInput.AddFileInfoList(_urls);
  fullInput.SetBranchStatus("caloJets*", 0);
  fullInput.SetBranchStatus("jptJets*", 0);
  fullInput.SetBranchStatus("photons*", 0);
  fullInput.SetBranchStatus("photons_photons*", 1);
  fullInput.SetBranchStatus("electrons*", 0);
  fullInput.SetBranchStatus("electrons_gsfElectrons*", 1);
  fullInput.SetBranchStatus("muons*", 0);
  fullInput.SetBranchStatus("muons_muons*", 1);

  fullEvent->setInput(fullInput);

  susy::TriggerEvent* triggerInput(new susy::TriggerEvent);
  susy::TriggerEvent* triggerOutput(0);
  if(triggerInput->bindTree(&fullInput, "susyEvents", "susyTriggers")){
    triggerOutput = new susy::TriggerEvent;
    triggerOutput->bookTrees(_outputNames->At(1)->GetName());
  }
  else{
    delete triggerInput;
    triggerInput = 0;
  }

  outputFile->cd();
  TTree* outTree(new TTree("susyTree", "SUSY Event"));
  outTree->SetAutoSave(10000000);
  fullEvent->addOutput(*outTree);

  susy::GoodLumis lumiFilter;
  lumiFilter.parseJSON(_goodLumis);

  long iEntry(0);
  while(event->getEntry(iEntry++) != 0){
    if(iEntry % 1000 == 1) std::cout << "Processing event " << iEntry << std::endl;

    if(event->isRealData && !lumiFilter.isGoodLumi(event->runNumber, event->luminosityBlockNumber)) continue;

    if(nTrig != 0){
      int iTrig(0);
      for(; iTrig != nTrig; ++iTrig)
        if(event->hltMap.pass(triggerNames.at(iTrig))) break;
      if(iTrig == nTrig) continue;
    }

    bool pass(false);

    if(skimPhoton){
      susy::PhotonCollection& photons(event->photons["photons"]);

      std::vector<bool> idResults;

      unsigned nPh(photons.size());
      unsigned iPh(0);
      for(; iPh != nPh; ++iPh){
        susy::Photon& photon(photons[iPh]);

        if(photon.momentum.Pt() < 15.) continue;
        if(photon.hadTowOverEm > 0.1) continue;

        float absEta(std::abs(photon.caloPosition.Eta()));

        if(absEta < susy::etaGapBegin && photon.sigmaIetaIeta < 0.014) break;
        else if(absEta > susy::etaGapEnd && absEta < susy::etaMax && photon.sigmaIetaIeta < 0.035) break;
      }
      pass = (iPh != nPh);
    }

    if(!pass && skimElectron){
      susy::ElectronCollection& electrons(event->electrons["gsfElectrons"]);

      unsigned nEl(electrons.size());
      unsigned iEl(0);
      for(; iEl != nEl; ++iEl){
        susy::Electron& electron(electrons[iEl]);

        if(electron.momentum.Pt() < 15.) continue;

        susy::ElectronVars vars(electron, *event);

        if(vars.isLoose) break;
      }
      pass = (iEl != nEl);
    }

    if(!pass && skimMuon){
      susy::MuonCollection& muons(event->muons["muons"]);

      unsigned nMu(muons.size());
      unsigned iMu(0);
      for(; iMu != nMu; ++iMu){
        susy::Muon& muon(muons[iMu]);

        if(muon.momentum.Pt() < 15.) continue;

        susy::MuonVars vars(muon, *event);

        if(vars.isLoose) break;
      }
      pass = (iMu != nMu);
    }

    if(!pass && skimGenLepton){
      std::vector<susy::Particle>& genParticles(event->genParticles);

      unsigned nG(genParticles.size());
      unsigned iG(0);
      for(; iG != nG; ++iG){
        susy::Particle& particle(genParticles[iG]);

        if(particle.status != 1 || particle.momentum.Pt() < 2.) continue;
        switch(std::abs(particle.pdgId)){
        case 11:
        case 13:
        case 15:
          break;
        default:
          continue;
        }
        break;
      }
      pass = (iG != nG);
    }

    if(!pass) continue;

    fullInput.GetEntry(iEntry - 1);

    if(triggerOutput) triggerOutput->copyEvent(*triggerInput);
    outTree->Fill();
  }

  delete triggerInput;

  if(triggerOutput) triggerOutput->write();
  delete triggerOutput;

  delete event;
  delete fullEvent;

  outputFile->cd();
  outTree->Write();
  delete outputFile;
}
