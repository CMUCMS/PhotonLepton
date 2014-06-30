#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"

#include "../CommonCode/ObjectTree.h"
#include "../CommonCode/Utilities.h"
#include "../CommonCode/SimpleEventProducer.h"

#include <iostream>
#include <cmath>

enum InputType {
  kDataE,
  kDataM,
  kMC,
  nInputTypes
};

enum EventType {
  kEEG,
  kMMG,
  kEEF,
  kMMF,
  nEventTypes
};

void
llgcomp(char const* _sourceName, char const* _outputName, unsigned _inputType, double _xsWeight = 1.)
{
  TChain eventTree("eventVars");
  TChain objTree("allObjects");
  TChain genTree("eventVars");

  eventTree.Add(_sourceName);
  objTree.Add(_sourceName);
  genTree.Add(_sourceName);

  eventTree.SetBranchStatus("*", 0);
  eventTree.SetBranchStatus("PhotonAndElectron", 1);
  eventTree.SetBranchStatus("PhotonAndMuon", 1);
  eventTree.SetBranchStatus("FakePhotonAndElectron", 1);
  eventTree.SetBranchStatus("FakePhotonAndMuon", 1);
  if(eventTree.GetBranch("puWeight")) eventTree.SetBranchStatus("puWeight", 1);

  objTree.SetBranchStatus("*", 0);
  objTree.SetBranchStatus("photon*", 1);
  if(_inputType == kMC || _inputType == kDataE) objTree.SetBranchStatus("electron*", 1);
  if(_inputType == kMC || _inputType == kDataM) objTree.SetBranchStatus("muon*", 1);

  genTree.SetBranchStatus("*", 0);
  if(genTree.GetBranch("gen.size")) genTree.SetBranchStatus("gen.*", 1);

  bool PhotonAndElectron;
  bool PhotonAndMuon;
  bool FakePhotonAndElectron;
  bool FakePhotonAndMuon;
  float puWeight(1.);

  susy::PhotonVarsArray photons;
  susy::ElectronVarsArray electrons;
  susy::MuonVarsArray muons;
  bool photon_isCand[susy::NMAX];
  bool photon_isFake[susy::NMAX];
  bool electron_isCand[susy::NMAX];
  bool muon_isCand[susy::NMAX];

  susy::SimpleEventProducer::EventVars eventVars;
  unsigned char photon_genMatch[susy::NMAX];
  bool electron_genMatch[susy::NMAX];
  bool muon_genMatch[susy::NMAX];

  eventTree.SetBranchAddress("PhotonAndElectron", &PhotonAndElectron);
  eventTree.SetBranchAddress("PhotonAndMuon", &PhotonAndMuon);
  eventTree.SetBranchAddress("FakePhotonAndElectron", &FakePhotonAndElectron);
  eventTree.SetBranchAddress("FakePhotonAndMuon", &FakePhotonAndMuon);
  if(eventTree.GetBranchStatus("puWeight")) eventTree.SetBranchAddress("puWeight", &puWeight);

  photons.setAddress(objTree);
  objTree.SetBranchAddress("photon.isCand", photon_isCand);
  objTree.SetBranchAddress("photon.isFake", photon_isFake);
  electrons.setAddress(objTree);
  objTree.SetBranchAddress("electron.isCand", electron_isCand);
  muons.setAddress(objTree);
  objTree.SetBranchAddress("muon.isCand", muon_isCand);

  if(_inputType == kMC) eventVars.setAddress(genTree);

  TH2* photonIDSF(0);
  TH2* electronIDSF(0);
  TH2* muonIDSF(0);
  TH2* photonHLTESF(0);
  TH2* photonHLTMSF(0);
  TH2* electronHLTSF(0);
  TH2* muonHLTSF(0);
  double muegHLTSF(0);

  if(_inputType == kMC){
    TFile* idSFSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiency/scalefactors.root"));
    photonIDSF = static_cast<TH2*>(idSFSource->Get("photon"));
    electronIDSF = static_cast<TH2*>(idSFSource->Get("electron"));
    muonIDSF = static_cast<TH2*>(idSFSource->Get("muon"));

    TFile* hltSFSource(TFile::Open("/afs/cern.ch/user/y/yiiyama/output/GammaL/main/hltEfficiency/scalefactors.root"));
    photonHLTESF = static_cast<TH2*>(hltSFSource->Get("photon_e"));
    photonHLTMSF = static_cast<TH2*>(hltSFSource->Get("photon_mu"));
    electronHLTSF = static_cast<TH2*>(hltSFSource->Get("electron"));
    muonHLTSF = static_cast<TH2*>(hltSFSource->Get("muon"));

    TGraph* mueg(static_cast<TGraph*>(hltSFSource->Get("mueg")));
    muegHLTSF = mueg->GetY()[0];
  }

  TFile* outputFile(TFile::Open(_outputName, "recreate"));
  TTree* output(new TTree("llgTree", "llgamma"));

  //mkbranch
  unsigned char eventType;
  float eventWeight;
  float gPt;
  float gEta;
  float gPhi;
  float gW(1.);
  bool gGenMatch(0);
  float l1Pt;
  float l1Eta;
  float l1Phi;
  float l1DR;
  float l1W(1.);
  bool l1GenMatch(0);
  float l2Pt;
  float l2Eta;
  float l2Phi;
  float l2DR;
  float l2W(1.);
  bool l2GenMatch(0);
  float llMass;
  float llgMass;
  //mkbranch

  output->Branch("eventType", &eventType, "eventType/b");
  output->Branch("eventWeight", &eventWeight, "eventWeight/F");
  output->Branch("puWeight", &puWeight, "puWeight/F");
  output->Branch("gPt", &gPt, "gPt/F");
  output->Branch("gEta", &gEta, "gEta/F");
  output->Branch("gPhi", &gPhi, "gPhi/F");
  output->Branch("gW", &gW, "gW/F");
  output->Branch("gGenMatch", &gGenMatch, "gGenMatch/O");
  output->Branch("l1Pt", &l1Pt, "l1Pt/F");
  output->Branch("l1Eta", &l1Eta, "l1Eta/F");
  output->Branch("l1Phi", &l1Phi, "l1Phi/F");
  output->Branch("l1DR", &l1DR, "l1DR/F");
  output->Branch("l1W", &l1W, "l1W/F");
  output->Branch("l1GenMatch", &l1GenMatch, "l1GenMatch/O");
  output->Branch("l2Pt", &l2Pt, "l2Pt/F");
  output->Branch("l2Eta", &l2Eta, "l2Eta/F");
  output->Branch("l2Phi", &l2Phi, "l2Phi/F");
  output->Branch("l2DR", &l2DR, "l2DR/F");
  output->Branch("l2W", &l2W, "l2W/F");
  output->Branch("l2GenMatch", &l2GenMatch, "l2GenMatch/O");
  output->Branch("llMass", &llMass, "llMass/F");
  output->Branch("llgMass", &llgMass, "llgMass/F");

  long iEntry(0);
  while(eventTree.GetEntry(iEntry++) > 0){
    if(iEntry % 10000 == 1) (std::cout << "\r" << iEntry).flush();

    if(PhotonAndElectron || FakePhotonAndElectron){
      if(_inputType == kDataM) continue;
    }
    else if(PhotonAndMuon || FakePhotonAndMuon){
      if(_inputType == kDataE) continue;
    }
    else continue;

    objTree.GetEntry(iEntry - 1);

    if(_inputType == kMC){
      genTree.GetEntry(iEntry - 1);
      std::fill_n(photon_genMatch, susy::NMAX, 0);
      std::fill_n(electron_genMatch, susy::NMAX, false);
      std::fill_n(muon_genMatch, susy::NMAX, false);
      susy::genMatch(eventVars, &photons, &electrons, &muons, photon_genMatch, electron_genMatch, muon_genMatch);
    }

    for(unsigned iP(0); iP != photons.size; ++iP){
      if(!photon_isCand[iP] && !photon_isFake[iP]) continue;

      gPt = photons.pt[iP];
      gEta = photons.eta[iP];
      gPhi = photons.phi[iP];

      if(_inputType == kMC){
	gGenMatch = photon_genMatch[iP] == 1;

        TVector3 caloPosition(photons.caloX[iP], photons.caloY[iP], photons.caloZ[iP]);
        double absEta(std::abs(caloPosition.Eta()));
        gW = photonIDSF->GetBinContent(photonIDSF->FindBin(gPt, absEta));
        if(PhotonAndElectron)
          gW *= photonHLTESF->GetBinContent(photonHLTESF->FindBin(gPt, absEta));
        else{
          gW *= photonHLTMSF->GetBinContent(photonHLTMSF->FindBin(gPt, absEta));
          gW *= muegHLTSF;
        }
      }

      TLorentzVector pG(photons.px[iP], photons.py[iP], photons.pz[iP], photons.energy[iP]);

      if(PhotonAndElectron || FakePhotonAndElectron){

        for(unsigned iL1(0); iL1 != electrons.size; ++iL1){
          if(!electron_isCand[iL1]) continue;

          TLorentzVector pL1(electrons.px[iL1], electrons.py[iL1], electrons.pz[iL1], electrons.energy[iL1]);
	  l1DR = pL1.DeltaR(pG);

          l1Pt = electrons.pt[iL1];
          l1Eta = electrons.eta[iL1];
          l1Phi = electrons.phi[iL1];

          if(_inputType == kMC){
	    l1GenMatch = electron_genMatch[iL1];

            TVector3 caloPosition(electrons.caloX[iL1], electrons.caloY[iL1], electrons.caloZ[iL1]);
            double absEta(std::abs(caloPosition.Eta()));
            l1W = electronIDSF->GetBinContent(electronIDSF->FindBin(l1Pt, absEta));
            l1W *= electronHLTSF->GetBinContent(electronHLTSF->FindBin(l1Pt, absEta));
          }

          for(unsigned iL2(iL1 + 1); iL2 != electrons.size; ++iL2){
            if(!electrons.isMedium[iL2]) continue;

            TLorentzVector pL2(electrons.px[iL2], electrons.py[iL2], electrons.pz[iL2], electrons.energy[iL2]);
	    l2DR = pL2.DeltaR(pG);

            eventType = photon_isCand[iP] ? kEEG : kEEF;

            l2Pt = electrons.pt[iL2];
            l2Eta = electrons.eta[iL2];
            l2Phi = electrons.phi[iL2];

            if(_inputType == kMC){
	      l2GenMatch = electron_genMatch[iL2];

              TVector3 caloPosition(electrons.caloX[iL2], electrons.caloY[iL2], electrons.caloZ[iL2]);
              double absEta(std::abs(caloPosition.Eta()));
              
              int binX(electronIDSF->GetXaxis()->FindFixBin(l2Pt));
              if(binX == 0){
                if(l2Pt < 15.){
                  if(absEta < 0.8) l2W = 0.865;
                  else if(absEta < 1.442) l2W = 0.967;
                  else if(absEta < 1.556) l2W = 1.064;
                  else if(absEta < 2.) l2W = 0.939;
                  else l2W = 1.05;
                }
                else{
                  if(absEta < 0.8) l2W = 0.958;
                  else if(absEta < 1.442) l2W = 0.971;
                  else if(absEta < 1.556) l2W = 0.902;
                  else if(absEta < 2.) l2W = 0.897;
                  else l2W = 0.941;
                }
              }
              else
                l2W = electronIDSF->GetBinContent(binX, electronIDSF->GetYaxis()->FindFixBin(absEta));
            }
     
            llMass = (pL1 + pL2).M();
            llgMass = (pG + pL1 + pL2).M();

	    eventWeight = _xsWeight * puWeight * gW * l1W * l2W;

            output->Fill();
          }
	}
      }
      else{
        for(unsigned iL1(0); iL1 != muons.size; ++iL1){
          if(!muon_isCand[iL1]) continue;

          TLorentzVector pL1(muons.px[iL1], muons.py[iL1], muons.pz[iL1], muons.energy[iL1]);
	  l1DR = pL1.DeltaR(pG);

          l1Pt = muons.pt[iL1];
          l1Eta = muons.eta[iL1];
          l1Phi = muons.phi[iL1];

          if(_inputType == kMC){
	    l1GenMatch = muon_genMatch[iL1];

            l1W = muonIDSF->GetBinContent(muonIDSF->FindBin(l1Pt, std::abs(l1Eta)));
            l1W *= muonHLTSF->GetBinContent(muonHLTSF->FindBin(l1Pt, std::abs(l1Eta)));
          }

          for(unsigned iL2(iL1 + 1); iL2 != muons.size; ++iL2){
            if(!muons.isTight[iL2]) continue;

            TLorentzVector pL2(muons.px[iL2], muons.py[iL2], muons.pz[iL2], muons.energy[iL2]);
	    l2DR = pL2.DeltaR(pG);

            eventType = photon_isCand[iP] ? kMMG : kMMF;

            l2Pt = muons.pt[iL2];
            l2Eta = muons.eta[iL2];
            l2Phi = muons.phi[iL2];

            if(_inputType == kMC){
	      l2GenMatch = muon_genMatch[iL2];

              l2W = muonIDSF->GetBinContent(muonIDSF->FindBin(25., std::abs(l2Eta)));
	    }

            llMass = (pL1 + pL2).M();
            llgMass = (pG + pL1 + pL2).M();

	    eventWeight = _xsWeight * puWeight * gW * l1W * l2W;

            output->Fill();
          }
        }
      }
    }
  }

  std::cout << std::endl;

  outputFile->cd();
  output->Write();
  delete outputFile;
}
