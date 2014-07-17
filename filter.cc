#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TChainElement.h"
#include "TMap.h"

#include "SusyEvent.h"
#include "SusyTriggerEvent.h"

#include "../CommonCode/Utilities.h"
#include "../CommonCode/ObjectSelector.h"
#include "../CommonCode/ObjectTree.h"
#include "../CommonCode/SimpleEventProducer.h"
#include "../CommonCode/PFParticleBugFix.h"
#include "../CommonCode/produceSimpleTree.cc"

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <set>
#include <map>
#include <bitset>

bool const USEBEAMSPOTIP(false);

enum FilterTypes {
  kPhotonAndElectron,
  kPhotonAndMuon,
  kElePhotonAndElectron,
  kElePhotonAndMuon,
  kFakePhotonAndElectron,
  kFakePhotonAndMuon,
  kPhotonAndFakeElectron,
  kPhotonAndFakeMuon,
  kFakePhotonAndFakeElectron,
  kFakePhotonAndFakeMuon,
  nFilterTypes
};

TString filterNames[] = {
  "PhotonAndElectron",
  "PhotonAndMuon",
  "ElePhotonAndElectron",
  "ElePhotonAndMuon",
  "FakePhotonAndElectron",
  "FakePhotonAndMuon",
  "PhotonAndFakeElectron",
  "PhotonAndFakeMuon",
  "FakePhotonAndFakeElectron",
  "FakePhotonAndFakeMuon"
};

enum CountPoints {
  kAllEvents,
  kGoodLumi,
  kMetFilter,
  kHLT,
  kGoodVertex,
  kRadiationVeto,
  kGoodElectron,
  kFakeElectron,
  kGoodMuon,
  kFakeMuon,
  kFlavorConflict,
  kChargedHadronVeto,
  kGoodPhoton,
  kFakePhoton,
  kElePhoton,
  kPassOneFilter,
  nCountPoints
};

TString countPoints[] = {
  "AllEvents",
  "GoodLumi",
  "MetFilter",
  "HLT",
  "GoodVertex",
  "RadiationVeto",
  "GoodElectron",
  "FakeElectron",
  "GoodMuon",
  "FakeMuon",
  "FlavorConflict",
  "ChargedHadronVeto",
  "GoodPhoton",
  "FakePhoton",
  "ElePhoton",
  "PassOneFilter",
  "PassPhotonAndElectron",
  "PassPhotonAndMuon",
  "PassElePhotonAndElectron",
  "PassElePhotonAndMuon",
  "PassFakePhotonAndElectron",
  "PassFakePhotonAndMuon",
  "PassPhotonAndFakeElectron",
  "PassPhotonAndFakeMuon"
};


class PhotonLeptonFilter : public SimpleTreeProducer {
public:
  PhotonLeptonFilter();
  ~PhotonLeptonFilter();
  bool initialize(char const*, char const*, double = 0.);
  void addInput(char const*, char const* = "");
  bool run();
  void clearInput();
  bool finalize();
private:
  bool readConfig_(char const*, std::map<TString, TString>&);

  TChain input_;
  TChain triggerInput_;
  TChain fullInput_;
  TChain fullTriggerInput_;

  unsigned eventCounter_[nCountPoints + nFilterTypes];
  std::bitset<nFilterTypes> useEvents_;

  susy::GoodLumis goodLumis_;

  bool filterResults_[nFilterTypes];

  bool ph_isCand_[susy::NMAX];
  bool ph_isFake_[susy::NMAX];
  bool ph_isEle_[susy::NMAX];
  bool el_isCand_[susy::NMAX];
  bool el_isFake_[susy::NMAX];
  bool mu_isCand_[susy::NMAX];
  bool mu_isFake_[susy::NMAX];
  bool jt_isCand_[susy::NMAX];

  bool applyHLTCut_;
  double radiationVetoThreshold_;
};

PhotonLeptonFilter::PhotonLeptonFilter() :
  SimpleTreeProducer(),
  input_("susyTree"),
  triggerInput_("triggerEventTree"),
  fullInput_("susyTree"),
  fullTriggerInput_("triggerEventTree"),
  useEvents_(),
  goodLumis_(),
  applyHLTCut_(true),
  radiationVetoThreshold_(0.)
{
}

PhotonLeptonFilter::~PhotonLeptonFilter()
{
}

bool
PhotonLeptonFilter::initialize(char const* _outputDir, char const* _configFileNames, double _radiationVetoThreshold/* = 0.*/)
{
  /* CONFIGURE */

  std::map<TString, TString> configRecords;

  TPRegexp configPat("^[ ]*([A-Z0-9_]+)[ ]*=[ ]*([^ ].*)$");
  std::string buf;
  std::stringstream bufs;
  TString line;

  TObjArray* configFileNames(TString(TString(_configFileNames).Strip(TString::kBoth, ' ')).Tokenize(" "));

  for(int iF(0); iF != configFileNames->GetEntries(); ++iF){
    std::ifstream configFile(configFileNames->At(iF)->GetName());
    if(!configFile.is_open()){
      delete configFileNames;
      std::cerr << "Cannot open config file" << std::endl;
      if(throw_) throw std::invalid_argument(_configFileNames);
      else return false;
    }

    while(true){
      std::getline(configFile, buf);
      if(!configFile.good()) break; // newline at the end necessary
      line = buf;

      TObjArray* matches(configPat.MatchS(line));
      if(matches->GetEntries() == 3){
        std::cout << matches->At(1)->GetName() << " = " << matches->At(2)->GetName() << std::endl;
        configRecords[matches->At(1)->GetName()] = matches->At(2)->GetName();
      }
      delete matches;
    }

    configFile.close();
  }

  delete configFileNames;

  if(!goodLumis_.parseJSON(configRecords["JSON"])){
    std::cerr << "Failed to parse JSON" << std::endl;
    if(throw_) throw std::invalid_argument(configRecords["JSON"].Data());
    else return false;
  }

  TObjArray* gridParams(configRecords["GRIDPARAMS"].Tokenize(" "));
  for(int iP(0); iP != gridParams->GetEntries(); ++iP)
    eventProducer_.addGridParam(gridParams->At(iP)->GetName());
  delete gridParams;

  TObjArray* filters(configRecords["FILTERS"].Tokenize(" "));
  for(int iF(0); iF != filters->GetEntries(); ++iF){
    TString filter(filters->At(iF)->GetName());
    TString* pF(std::find(filterNames, filterNames + nFilterTypes, filter));
    if(pF == filterNames + nFilterTypes){
      std::cerr << ("Filter type " + buf + " not defined") << std::endl;
      if(throw_) throw std::invalid_argument(filter.Data());
      else return false;
    }
    useEvents_.set(pF - filterNames);
  }
  delete filters;

  if(useEvents_.none()) useEvents_.set();

  if(configRecords.find("HLTCUT") != configRecords.end() && configRecords["HLTCUT"] == "No")
    applyHLTCut_ = false;

  /* INITIALIZE OUTPUT */

  TString outputName(_outputDir);
  if(!outputName.Contains(".root")) outputName += "/skim.root";

  if(!SimpleTreeProducer::initialize(outputName, configRecords["PUSCENARIO"], false, true, true, true)) return false;
  eventProducer_.setPhotonId(susy::PhLoose12Pix);

  allObjTree_->Branch("photon.isCand", ph_isCand_, "isCand[photon.size]/O");
  allObjTree_->Branch("photon.isFake", ph_isFake_, "isFake[photon.size]/O");
  allObjTree_->Branch("photon.isEle", ph_isEle_, "isEle[photon.size]/O");
  allObjTree_->Branch("electron.isCand", el_isCand_, "isCand[electron.size]/O");
  allObjTree_->Branch("electron.isFake", el_isFake_, "isFake[electron.size]/O");
  allObjTree_->Branch("muon.isCand", mu_isCand_, "isCand[muon.size]/O");
  allObjTree_->Branch("muon.isFake", mu_isFake_, "isFake[muon.size]/O");
  allObjTree_->Branch("jet.isCand", jt_isCand_, "isCand[jet.size]/O");

  for(unsigned iF(0); iF != nFilterTypes; ++iF)
    evtTree_->Branch(filterNames[iF], filterResults_ + iF, filterNames[iF] + "/O");

  outputFile_->cd();

  for(std::map<TString, TString>::iterator cItr(configRecords.begin()); cItr != configRecords.end(); ++cItr){
    if(cItr->second != ""){
      TObjString cLine(cItr->first + " = " + cItr->second);
      cLine.Write();
    }
  }

  /* SETUP RADIATION VETO */
  // veto event if an isolated photon with Pt > |threshold| is found
  // non-isolated photon has a sibling (shares motherIndex) particle that is not a neutrino and has Pt > 2 GeV within dR < 0.3
  // when threshold < 0, veto ISR only
  radiationVetoThreshold_ = _radiationVetoThreshold;

  /* INITIALIZE COUNTER */

  std::fill_n(eventCounter_, nCountPoints + nFilterTypes, 0);

  return true;
}

void
PhotonLeptonFilter::addInput(char const* _source, char const* _triggerSource/* = ""*/)
{
  input_.Add(_source);
  fullInput_.Add(_source);

  if(!_triggerSource || std::strlen(_triggerSource) == 0){
    TString inputName(_source);
    inputName.ReplaceAll("susyEvents", "susyTriggers");
    triggerInput_.Add(inputName);
    fullTriggerInput_.Add(inputName);
  }
  else{
    triggerInput_.Add(_triggerSource);
    fullTriggerInput_.Add(_triggerSource);
  }
}

bool
PhotonLeptonFilter::run()
{
  /* INITIALIZE INPUT FOR FILTERING */

  input_.SetBranchStatus("*", 0);
  input_.SetBranchStatus("eventNumber", 1);
  input_.SetBranchStatus("runNumber", 1);
  input_.SetBranchStatus("luminosityBlockNumber", 1);
  input_.SetBranchStatus("metFilter*", 1);
  input_.SetBranchStatus("hlt*", 1);
  input_.SetBranchStatus("pfParticles*", 1);
  if(radiationVetoThreshold_ != 0. && input_.GetBranch("genParticles")) input_.SetBranchStatus("genParticles*", 1);
  susy::ObjectTree::setBranchStatus(input_, true, true, true, false, true); // Photon, Electron, Muon, Jet, Vertex
  if(USEBEAMSPOTIP) input_.SetBranchStatus("beamSpot*", 1);

  susy::Event event;
  susy::TriggerEvent triggerEvent;
  if(!triggerEvent.bindTree(&input_, &triggerInput_)){
    if(throw_) throw std::runtime_error("bindTree");
    else return false;
  }

  event.setInput(input_);

  /* INITIALIZE INPUT TO COPY */

  susy::Event fullEvent;
  susy::TriggerEvent fullTriggerEvent;
  if(!initInput_(fullInput_, fullEvent, fullTriggerInput_, fullTriggerEvent)){
    std::cerr << "Input incompatible" << std::endl;
    if(throw_) throw std::runtime_error("");
    else return false;
  }

  /* TRIGGERS */

  TString electronHLT[] = {
    "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50"
  };
  unsigned const nElectronHLT(sizeof(electronHLT) / sizeof(TString));

  TString muonHLT[] = {
    "HLT_Mu22_Photon22_CaloIdL"
  };
  unsigned const nMuonHLT(sizeof(muonHLT) / sizeof(TString));

  std::vector<TString> electronHLTFilters[nElectronHLT];
  electronHLTFilters[0].push_back("hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded");

  std::vector<TString> muonHLTFilters[nMuonHLT];
  muonHLTFilters[0].push_back("hltL1Mu3p5EG12L3Filtered22");

  std::vector<TString> photonHLTFiltersE[nElectronHLT];
  photonHLTFiltersE[0].push_back("hltEG36CaloId10Iso50HcalIsoLastFilter");
  photonHLTFiltersE[0].push_back("hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded");

  std::vector<TString> photonHLTFiltersM[nMuonHLT];
  photonHLTFiltersM[0].push_back("hltMu22Photon22CaloIdLHEFilter");

  /* OBJECT ID MASKS */

  std::bitset<susy::nMuonCriteria> muIdResults;
  std::bitset<susy::nMuonCriteria> muNoIP(susy::ObjectSelector::muReferences[susy::MuTight12]);
  muNoIP.reset(susy::MuDxy);
  muNoIP.reset(susy::MuDz);
  std::bitset<susy::nElectronCriteria> elIdResults;
  std::bitset<susy::nElectronCriteria> elNoIP(susy::ObjectSelector::elReferences[susy::ElMedium12]);
  elNoIP.reset(susy::ElD0);
  elNoIP.reset(susy::ElDZ);
  std::bitset<susy::nPhotonCriteria> phIdResults;
  std::bitset<susy::nPhotonCriteria> phNoVeto(susy::ObjectSelector::phReferences[susy::PhLoose12Pix]);
  phNoVeto.reset(susy::PhElectronVeto);

  /* START LOOP */
  
  long iEntry(0);
  int nRead(0);
  
  while((nRead = event.getEntry(iEntry++)) != 0){
    try{
      if(nRead < 0)
        throw susy::Exception(susy::Exception::kIOError, "Corrupt input");
    
      if(iEntry % 1000 == 0) std::cout << "Analyzing event " << iEntry << std::endl;

      ++eventCounter_[kAllEvents];

      if(!goodLumis_.isGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;

      ++eventCounter_[kGoodLumi];

      if(!event.passMetFilters()) continue;

      ++eventCounter_[kMetFilter];

      std::bitset<nElectronHLT> passElectronHLT;
      std::bitset<nMuonHLT> passMuonHLT;
      passElectronHLT.set();
      passMuonHLT.set();

      if(applyHLTCut_){
        for(unsigned iHLT(0); iHLT != nElectronHLT; ++iHLT)
          if(!event.hltMap.pass(electronHLT[iHLT] + "_v*")) passElectronHLT.reset(iHLT);
        for(unsigned iHLT(0); iHLT != nMuonHLT; ++iHLT)
          if(!event.hltMap.pass(muonHLT[iHLT] + "_v*")) passMuonHLT.reset(iHLT);

        if(passElectronHLT.none() && passMuonHLT.none()) continue;
      }

      ++eventCounter_[kHLT];

      unsigned nV(event.vertices.size());
      unsigned iV(0);
      for(; iV != nV; ++iV){
        susy::VertexVars vars(event.vertices[iV]);
        if(vars.isGood) break;
      }
      if(iV == nV) continue;

      ++eventCounter_[kGoodVertex];

      if(radiationVetoThreshold_ != 0.){
        double ptThreshold(std::abs(radiationVetoThreshold_));

        susy::ParticleCollection const& genParticles(event.genParticles);
        unsigned nG(genParticles.size());

        unsigned iG(0);
        for(; iG != nG; ++iG){
          susy::Particle const& particle(genParticles[iG]);
          if(particle.status != 1 || particle.pdgId != 22 || particle.momentum.Pt() < ptThreshold) continue; // not a photon

          if(particle.mother && std::abs(particle.mother->pdgId) > 99) continue; // is fragmentation
          
          // unsigned iI(0);
          // for(; iI != nG; ++iI){
          //   if(iI == iG) continue;
          //   susy::Particle const& isoP(genParticles[iI]);
          //   if(isoP.status == 1 && !(isoP.charge == 0 && std::abs(isoP.pdgId) < 20) && isoP.momentum.Pt() > 2. && isoP.momentum.DeltaR(particle.momentum) < 0.3) break;
          // }
          // if(iI != nG) continue; // is not isolated

          double iso(0.);
          for(unsigned iI(0); iI != nG; ++iI){
            if(iI == iG) continue;
            susy::Particle const& isoP(genParticles[iI]);
            unsigned isoId(std::abs(isoP.pdgId));
            if(isoP.status == 1 && (isoId < 11 || isoId > 16) && isoP.momentum.Pt() > 2. && isoP.momentum.DeltaR(particle.momentum) < 0.3){
              iso += isoP.momentum.Pt();
              if(iso > 10.) break;
            }
          }
          if(iso > 10.) continue; // is not isolated

          if(radiationVetoThreshold_ < 0.){
            short idx(particle.motherIndex);
            while(idx != -1 && genParticles[idx].pdgId != 23 && std::abs(genParticles[idx].pdgId) != 24) idx = genParticles[idx].motherIndex;
            if(idx != -1) continue; // is FSR
          }

          break;
        }

        if(iG != nG) continue;
      }

      if(eventCounter_[kRadiationVeto] < 10) std::cout << std::endl << iEntry - 1 << std::endl;
      ++eventCounter_[kRadiationVeto];

      std::vector<const susy::Electron*> electrons(eventProducer_.sortElectrons(event.electrons["gsfElectrons"]));
      unsigned nEl(electrons.size());

      unsigned nCandElectron(0);
      unsigned nFakeElectron(0);

      /* SELECT ELECTRONS */

      std::vector<susy::TriggerObjectCollection> electronHLTObjects[nElectronHLT];
      if(applyHLTCut_){
        for(unsigned iHLT(0); iHLT != nElectronHLT; ++iHLT){
          if(!passElectronHLT[iHLT]) continue; // to save time
          for(unsigned iF(0); iF != electronHLTFilters[iHLT].size(); ++iF)
            electronHLTObjects[iHLT].push_back(triggerEvent.getFilterObjects(electronHLTFilters[iHLT][iF]));
        }
      }

      //      elLooseNoIso.assign(nEl, false);

      for(unsigned iEl(0); iEl < nEl; ++iEl){
        el_isCand_[iEl] = false;
        el_isFake_[iEl] = false;

        susy::Electron const& el(*electrons[iEl]);

        if(el.momentum.Pt() < 25.) continue;

        susy::ElectronVars vars(el, event);

        if(vars.iSubdet == -1) continue;

        if(applyHLTCut_){
          TLorentzVector dSC(el.superCluster->position, 0.);
          unsigned iHLT(0);
          for(; iHLT != nElectronHLT; ++iHLT){
            if(!passElectronHLT[iHLT]) continue;
            unsigned iF(0);
            for(; iF != electronHLTFilters[iHLT].size(); ++iF){
              susy::TriggerObjectCollection& objects(electronHLTObjects[iHLT][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(dSC) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == electronHLTFilters[iHLT].size()) break; // all filters of the path had a match
          }
          if(iHLT == nElectronHLT) continue;
        }

	bool isMedium(false);
	if(USEBEAMSPOTIP){
	  susy::ObjectSelector::isGoodElectron(vars, susy::ElMedium12, &elIdResults);
	  if((elIdResults & elNoIP) == elNoIP)
	    isMedium = std::abs(el.gsfTrack->dxy(event.beamSpot)) < 0.2;
	}
	else
	  isMedium = vars.isMedium;

        if(isMedium){
          el_isCand_[iEl] = true;
	  ++nCandElectron;
	}
        else{
          el_isFake_[iEl] = true;
	  ++nFakeElectron;
	}
      }

      if(passMuonHLT.none() && nCandElectron == 0 && nFakeElectron == 0) continue;

      if(nCandElectron != 0) ++eventCounter_[kGoodElectron];
      if(nFakeElectron != 0) ++eventCounter_[kFakeElectron];

      std::vector<const susy::Muon*> muons(eventProducer_.sortMuons(event.muons["muons"]));
      unsigned nMu(muons.size());

      unsigned nCandMuon(0);
      unsigned nFakeMuon(0);

      /* SELECT MUONS */

      std::vector<susy::TriggerObjectCollection> muonHLTObjects[nMuonHLT];
      if(applyHLTCut_){
        for(unsigned iHLT(0); iHLT != nMuonHLT; ++iHLT){
          if(!passMuonHLT[iHLT]) continue; // to save time
          for(unsigned iF(0); iF != muonHLTFilters[iHLT].size(); ++iF)
            muonHLTObjects[iHLT].push_back(triggerEvent.getFilterObjects(muonHLTFilters[iHLT][iF]));
        }
      }
   
      for(unsigned iMu(0); iMu < nMu; ++iMu){
        mu_isCand_[iMu] = false;
        mu_isFake_[iMu] = false;

        susy::Muon const& mu(*muons[iMu]);

        if(mu.momentum.Pt() < 25.) continue;

        susy::MuonVars vars(mu, event);

        if(vars.iSubdet == -1 || !vars.isLoose) continue;

        if(applyHLTCut_){
          unsigned iHLT(0);
          for(; iHLT != nMuonHLT; ++iHLT){
            if(!passMuonHLT[iHLT]) continue;
            unsigned iF(0);
            for(; iF != muonHLTFilters[iHLT].size(); ++iF){
              susy::TriggerObjectCollection& objects(muonHLTObjects[iHLT][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(mu.momentum) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == muonHLTFilters[iHLT].size()) break; // all filters of the path had a match
          }
          if(iHLT == nMuonHLT) continue;
        }

	bool isTight(false);
	if(USEBEAMSPOTIP){
	  susy::ObjectSelector::isGoodMuon(vars, susy::MuTight12, &muIdResults);
	  if((muIdResults & muNoIP) == muNoIP){
	    susy::Track const* track(vars.pt > 200. ? mu.highPtBestTrack : mu.bestTrack);
	    if(!track) track = mu.innerTrack;
	    if(track)
	      isTight = std::abs(track->dxy(event.beamSpot)) < 0.2;
	  }
	}
	else
	  isTight = vars.isTight;

        if(isTight){
          mu_isCand_[iMu] = true;
	  ++nCandMuon;
	}
        else{
          mu_isFake_[iMu] = true;
	  ++nFakeMuon;
	}
      }

      if(passElectronHLT.none() && nCandMuon == 0 && nFakeMuon == 0) continue;

      if(nCandMuon != 0) ++eventCounter_[kGoodMuon];
      if(nFakeMuon != 0) ++eventCounter_[kFakeMuon];

      if(nCandElectron != 0 && nCandMuon != 0) ++eventCounter_[kFlavorConflict];


      // bug fix for tag cms533v0 / cms538v1
      std::vector<susy::PFParticle const*> pfParticles(susy::cleanPFParticles(event.pfParticles));
      unsigned nPF(pfParticles.size());


      /* SELECT PHOTONS */

      std::vector<const susy::Photon*> photons(eventProducer_.sortPhotons(event.photons["photons"]));
      unsigned nPh(photons.size());

      std::vector<susy::TriggerObjectCollection> photonHLTObjectsE[nElectronHLT];
      std::vector<susy::TriggerObjectCollection> photonHLTObjectsM[nMuonHLT];
      if(applyHLTCut_){
        for(unsigned iHLT(0); iHLT != nElectronHLT; ++iHLT){
          if(!passElectronHLT[iHLT]) continue; // to save time
          for(unsigned iF(0); iF != photonHLTFiltersE[iHLT].size(); ++iF)
            photonHLTObjectsE[iHLT].push_back(triggerEvent.getFilterObjects(photonHLTFiltersE[iHLT][iF]));
        }
        for(unsigned iHLT(0); iHLT != nMuonHLT; ++iHLT){
          if(!passMuonHLT[iHLT]) continue; // to save time
          for(unsigned iF(0); iF != photonHLTFiltersM[iHLT].size(); ++iF)
            photonHLTObjectsM[iHLT].push_back(triggerEvent.getFilterObjects(photonHLTFiltersM[iHLT][iF]));
        }
      }

      bool passCHVeto(false);

      unsigned nCandPhoton(0);
      unsigned nFakePhoton(0);
      unsigned nElePhoton(0);

      for(unsigned iPh(0); iPh < nPh; ++iPh){
        ph_isCand_[iPh] = false;
        ph_isFake_[iPh] = false;
        ph_isEle_[iPh] = false;

        susy::Photon const& ph(*photons[iPh]);

        if(ph.momentum.Pt() < 25.) continue;
        if(std::abs(ph.caloPosition.Eta()) > susy::etaGapBegin) continue;

        /* require photons to match a PF object, but veto if the match is to a CH */
        /* better synchronizes PF MET and calo MET */

        unsigned iPF(0);
        for(; iPF != nPF; ++iPF){
          susy::PFParticle const& part(*pfParticles[iPF]);
          if(part.momentum.Pt() < 3.) continue;
          TVector3 dir(ph.caloPosition - part.vertex);
          if(std::abs(part.pdgId) == 211 &&
             std::abs(dir.Eta() - part.momentum.Eta()) < 0.005 &&
             std::abs(TVector2::Phi_mpi_pi(dir.Phi() - part.momentum.Phi())) < 0.02) break;
        }
        if(iPF != nPF) continue;

        passCHVeto = true;

        if(applyHLTCut_){
          TLorentzVector dSC(ph.superCluster->position, 0.);

          unsigned iHLTE(0);
          for(; iHLTE != nElectronHLT; ++iHLTE){
            if(!passElectronHLT[iHLTE]) continue;
            unsigned iF(0);
            for(; iF != photonHLTFiltersE[iHLTE].size(); ++iF){
              susy::TriggerObjectCollection& objects(photonHLTObjectsE[iHLTE][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(dSC) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == photonHLTFiltersE[iHLTE].size()) break; // all filters of the path matched
          }

          unsigned iHLTM(0);
          for(; iHLTM != nMuonHLT; ++iHLTM){
            if(!passMuonHLT[iHLTM]) continue;
            unsigned iF(0);
            for(; iF != photonHLTFiltersM[iHLTM].size(); ++iF){
              susy::TriggerObjectCollection& objects(photonHLTObjectsM[iHLTM][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(dSC) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == photonHLTFiltersM[iHLTM].size()) break; // all filters of the path matched
          }

          if(iHLTE == nElectronHLT && iHLTM == nMuonHLT) continue;
        }

        susy::PhotonVars vars(ph, event);
        bool isGood(susy::ObjectSelector::isGoodPhoton(vars, susy::PhLoose12Pix, &phIdResults));

        bool matchGSF(false);

        unsigned iL(0);
        for(; iL != nEl; ++iL){
          if(electrons[iL]->momentum.Pt() < 2.) continue;
          if(electrons[iL]->superClusterIndex == ph.superClusterIndex ||
             electrons[iL]->superCluster->position.DeltaR(ph.caloPosition) < 0.02) break;
        }
        matchGSF = iL != nEl;
        isGood = isGood && !matchGSF;

        if(isGood){
          ph_isCand_[iPh] = true;
	  ++nCandPhoton;
	}
        else if((phIdResults & phNoVeto) == phNoVeto){
          ph_isEle_[iPh] = true;
	  ++nElePhoton;
	}
        else if(!matchGSF && vars.nPixelSeeds == 0){
          ph_isFake_[iPh] = true;
	  ++nFakePhoton;
	}
      }

      if(passCHVeto) ++eventCounter_[kChargedHadronVeto];

      if(nCandPhoton != 0) ++eventCounter_[kGoodPhoton];
      if(nFakePhoton != 0) ++eventCounter_[kFakePhoton];
      if(nElePhoton != 0) ++eventCounter_[kElePhoton];

      /* DETERMINE RESULT OF EACH FILTER */

      filterResults_[kPhotonAndElectron] = passElectronHLT.any() && nCandPhoton != 0 && nCandElectron != 0;
      filterResults_[kPhotonAndMuon] = passMuonHLT.any() && nCandPhoton != 0 && nCandMuon != 0;
      filterResults_[kElePhotonAndElectron] = passElectronHLT.any() && nElePhoton != 0 && nCandElectron != 0;
      filterResults_[kElePhotonAndMuon] = passMuonHLT.any() && nElePhoton != 0 && nCandMuon != 0;
      filterResults_[kFakePhotonAndElectron] = passElectronHLT.any() && nFakePhoton != 0 && nCandElectron != 0;
      filterResults_[kFakePhotonAndMuon] = passMuonHLT.any() && nFakePhoton != 0 && nCandMuon != 0;
      filterResults_[kPhotonAndFakeElectron] = passElectronHLT.any() && nCandPhoton != 0 && nFakeElectron != 0;
      filterResults_[kPhotonAndFakeMuon] = passMuonHLT.any() && nCandPhoton != 0 && nFakeMuon != 0;

      if(filterResults_[kElePhotonAndElectron] && nElePhoton == 1 && nCandElectron == 1){
	// is the only elePhoton actually the candidate electron?
	unsigned iElePh(0);
	for(; iElePh != nPh; ++iElePh)
	  if(ph_isEle_[iElePh]) break;
	unsigned iCandEl(0);
	for(; iCandEl != nEl; ++iCandEl)
	  if(el_isCand_[iCandEl]) break;
	if(electrons[iCandEl]->superClusterIndex == photons[iElePh]->superClusterIndex ||
	   electrons[iCandEl]->superCluster->position.DeltaR(photons[iElePh]->caloPosition) < 0.02)
	  filterResults_[kElePhotonAndElectron] = false;
      }

      bool passOneFilter(false);
      for(unsigned iF(0); iF != nFilterTypes; ++iF){
        if(filterResults_[iF]){
          ++eventCounter_[nCountPoints + iF];
          if(useEvents_[iF]) passOneFilter = true;
        }
      }

      if(!passOneFilter) continue;

      ++eventCounter_[kPassOneFilter];

      /* EVENT PASSES AT LEAST ONE FILTER - LOAD FULL EVENT AND FILL OUTPUT */

      fullEvent.getEntry(iEntry - 1);

      std::vector<const susy::PFJet*> jets(eventProducer_.sortJets(fullEvent.pfJets["ak5"]));
      unsigned nJ(jets.size());

      for(unsigned iJ(0); iJ != nJ; ++iJ){
        jt_isCand_[iJ] = false;

        susy::PFJet const& jet(*jets[iJ]);

        susy::JetVars vars(jet, event);
        if(!vars.isLoose || !vars.passPUJetIdLoose) continue;
        if(vars.pt < 30. || std::abs(vars.eta) > 3.) continue;

        unsigned iC;

        iC = 0;
        for(; iC != nEl; ++iC)
          if(el_isCand_[iC] && electrons[iC]->momentum.DeltaR(jet.momentum) < 0.5) break;
        if(iC != nEl) continue;

        iC = 0;
        for(; iC != nMu; ++iC)
          if(mu_isCand_[iC] && muons[iC]->momentum.DeltaR(jet.momentum) < 0.5) break;
        if(iC != nMu) continue;

        iC = 0;
        for(; iC != nPh; ++iC)
          if((ph_isCand_[iC] || ph_isEle_[iC]) && photons[iC]->momentum.DeltaR(jet.momentum) < 0.5) break;
        if(iC != nPh) continue;

        jt_isCand_[iJ] = true;
      }

      eventProducer_.extractTriggerObjects(fullTriggerEvent);
      eventProducer_.produce(fullEvent);

      if(evtTree_->Fill() < 0)
        throw susy::Exception(susy::Exception::kIOError, "eventVars");
      if(allObjTree_->Fill() < 0)
        throw susy::Exception(susy::Exception::kIOError, "allObjects");
    }
    catch(std::exception& ex){
      if(processError_(ex, fullEvent, fullInput_)) continue;
      if(throw_) throw;
      else return false;
    }
  }

  return true;
}

void
PhotonLeptonFilter::clearInput()
{
  input_.Reset();
  triggerInput_.Reset();
  fullInput_.Reset();
  fullTriggerInput_.Reset();
}

bool
PhotonLeptonFilter::finalize()
{
  std::cout << "Cut flow: " << std::endl;
  for(unsigned iC(0); iC != nCountPoints + nFilterTypes; ++iC)
    std::cout << "[" << countPoints[iC] << "]: " << eventCounter_[iC] << std::endl;

  outputFile_->cd();
  evtTree_->Write();
  allObjTree_->Write();

  delete outputFile_;

  outputFile_ = 0;
  evtTree_ = 0;
  allObjTree_ = 0;

  return true;
}

void
filter(TString const& _configFileName, TString const& _inputPaths, TString const& _outputName, double _radiationVetoThreshold = 0.)
{
  PhotonLeptonFilter filter;
  filter.setThrow(true);
  filter.initialize(_outputName, _configFileName, _radiationVetoThreshold);
  filter.addInput(_inputPaths);
  filter.run();
  filter.clearInput();
  filter.finalize();
}
