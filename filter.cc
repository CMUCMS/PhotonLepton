#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TVector2.h"

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
#include <map>
#include <bitset>

#include "Toolset/GenTreeViewer/test/GenDecayFilterRA3.cc"

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
  kElePhotonAndFakeElectron,
  kElePhotonAndFakeMuon,
  kFakePhotonAndFakeElectron,
  kFakePhotonAndFakeMuon,
  nFilterTypes
};

TString filterNames[nFilterTypes] = {
  "PhotonAndElectron",
  "PhotonAndMuon",
  "ElePhotonAndElectron",
  "ElePhotonAndMuon",
  "FakePhotonAndElectron",
  "FakePhotonAndMuon",
  "PhotonAndFakeElectron",
  "PhotonAndFakeMuon",
  "ElePhotonAndFakeElectron",
  "ElePhotonAndFakeMuon",
  "FakePhotonAndFakeElectron",
  "FakePhotonAndFakeMuon"
};

enum Cuts {
  kAllEvents,
  kHLT,
  kGoodLumi,
  kMetFilter,
  kGoodVertex,
  kRadiationVeto,
  kGoodPhoton,
  kGoodLepton,
  nCuts
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

  TTree* effTree_;
  TTree* cutTree_;
  
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

  struct EffTreeVariables {
    //mkbranch
    char pdgId;
    bool pass;
    float pt;
    float eta;
    float phi;
    float genPt;
    float genIso;
    unsigned char nVtx;
    bool primVtx;
    //mkbranch
  } effVars_;

  double radiationVetoThreshold_;

  GenDecayFilterRA3* genFilter_;
};

PhotonLeptonFilter::PhotonLeptonFilter() :
  SimpleTreeProducer(),
  input_("susyTree"),
  triggerInput_("triggerEventTree"),
  fullInput_("susyTree"),
  fullTriggerInput_("triggerEventTree"),
  effTree_(0),
  cutTree_(0),
  goodLumis_(),
  radiationVetoThreshold_(0.),
  genFilter_(0)
{
}

PhotonLeptonFilter::~PhotonLeptonFilter()
{
  delete effTree_;
  delete cutTree_;
  delete genFilter_;
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

  if(configRecords["GENFILTER"].Length() != 0)
    genFilter_ = new GenDecayFilterRA3(configRecords["GENFILTER"]);

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

  if(configRecords["PUSCENARIO"] != ""){
    effTree_ = new TTree("effTree", "Efficiency Tree");
    effTree_->Branch("pdgId", &effVars_.pdgId, "pdgId/B");
    effTree_->Branch("pass", &effVars_.pass, "pass/O");
    effTree_->Branch("pt", &effVars_.pt, "pt/F");
    effTree_->Branch("eta", &effVars_.eta, "eta/F");
    effTree_->Branch("phi", &effVars_.phi, "phi/F");
    effTree_->Branch("genPt", &effVars_.genPt, "genPt/F");
    effTree_->Branch("genIso", &effVars_.genIso, "genIso/F");
    effTree_->Branch("nVtx", &effVars_.nVtx, "nVtx/b");
    effTree_->Branch("primVtx", &effVars_.primVtx, "primVtx/O");
  }

  cutTree_ = new TTree("cutTree", "Cutflow");
  cutTree_->Branch("run", 0, "run/i");
  cutTree_->Branch("lumi", 0, "lumi/i");
  cutTree_->Branch("event", 0, "event/i");
  cutTree_->Branch("cutflowE", 0, "cutflowE/b");
  cutTree_->Branch("cutflowM", 0, "cutflowM/b");

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
  if(input_.GetBranch("genParticles")) input_.SetBranchStatus("genParticles*", 1);
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

  /* SETUP CUTFLOW TREE */
  
  cutTree_->SetBranchAddress("run", &event.runNumber);
  cutTree_->SetBranchAddress("lumi", &event.luminosityBlockNumber);
  cutTree_->SetBranchAddress("event", &event.eventNumber);
  unsigned char cutflowE;
  unsigned char cutflowM;
  cutTree_->SetBranchAddress("cutflowE", &cutflowE);
  cutTree_->SetBranchAddress("cutflowM", &cutflowM);

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
  
  int nRead(0);
  for(long iEntry(0); (nRead = event.getEntry(iEntry)) != 0; cutTree_->Fill(), ++iEntry){
    try{
      if(nRead < 0)
        throw susy::Exception(susy::Exception::kIOError, "Corrupt input");
    
      if(iEntry % 1000 == 0) std::cout << "Analyzing event " << iEntry << std::endl;

      if(genFilter_ && !genFilter_->pass(event)) continue;

      cutflowE = 0;
      cutflowM = 0;

      std::bitset<nElectronHLT> passElectronHLT;
      std::bitset<nMuonHLT> passMuonHLT;

      for(unsigned iHLT(0); iHLT != nElectronHLT; ++iHLT)
        if(event.hltMap.pass(electronHLT[iHLT] + "_v*")) passElectronHLT.set(iHLT);
      for(unsigned iHLT(0); iHLT != nMuonHLT; ++iHLT)
        if(event.hltMap.pass(muonHLT[iHLT] + "_v*")) passMuonHLT.set(iHLT);

      if(passElectronHLT.any() && cutflowE == kHLT - 1) ++cutflowE;
      if(passMuonHLT.any() && cutflowM == kHLT - 1) ++cutflowM;

      if(!goodLumis_.isGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;

      if(cutflowE == kGoodLumi - 1) ++cutflowE;
      if(cutflowM == kGoodLumi - 1) ++cutflowM;

      if(!event.passMetFilters()) continue;

      if(cutflowE == kMetFilter - 1) ++cutflowE;
      if(cutflowM == kMetFilter - 1) ++cutflowM;

      unsigned nV(event.vertices.size());
      unsigned iV(0);
      for(; iV != nV; ++iV){
        susy::VertexVars vars(event.vertices[iV]);
        if(vars.isGood) break;
      }
      if(iV == nV) continue;

      if(cutflowE == kGoodVertex - 1) ++cutflowE;
      if(cutflowM == kGoodVertex - 1) ++cutflowM;

      if(radiationVetoThreshold_ != 0.){
        double ptThreshold(std::abs(radiationVetoThreshold_));

        susy::ParticleCollection const& genParticles(event.genParticles);
        unsigned nG(genParticles.size());

        unsigned iG(0);
        for(; iG != nG; ++iG){
          susy::Particle const& particle(genParticles[iG]);
          if(particle.status != 1 || particle.pdgId != 22 || particle.momentum.Pt() < ptThreshold) continue; // not a photon

          if(particle.mother && std::abs(particle.mother->pdgId) > 99) continue; // is fragmentation
          
          double iso(0.);
          for(unsigned iI(0); iI != nG; ++iI){
            if(iI == iG) continue;
            susy::Particle const& isoP(genParticles[iI]);
            if(isoP.status != 1) continue;
            unsigned isoId(std::abs(isoP.pdgId));
            if(isoId == 12 || isoId == 14 || isoId == 16 || isoId == 1000022 || isoId == 1000039) continue;
            if(isoP.momentum.DeltaR(particle.momentum) > 0.3) continue;
            iso += isoP.momentum.Pt();
            if(iso > 10.) break;
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
      if(cutflowE == kRadiationVeto - 1) ++cutflowE;
      if(cutflowM == kRadiationVeto - 1) ++cutflowM;

      std::vector<std::pair<susy::Particle const*, float> > genPhotons;
      std::vector<std::pair<susy::Particle const*, float> > genElectrons;
      std::vector<std::pair<susy::Particle const*, float> > genMuons;
      std::vector<susy::Photon const*> genMatchedPhotons;
      std::vector<susy::Electron const*> genMatchedElectrons;
      std::vector<susy::Muon const*> genMatchedMuons;

      if(effTree_){
        susy::ParticleCollection const& genParticles(event.genParticles);
        unsigned nG(genParticles.size());

        TVector3 const* genVertex(0);
        for(unsigned iG(0); iG != nG; ++iG){
          if(genParticles[iG].mother){
            genVertex = &genParticles[iG].vertex;
            break;
          }
        }

        double dZPrim(0.);
        effVars_.nVtx = 0;
        effVars_.primVtx = true;
        for(iV = 0; iV != nV; ++iV){
          susy::VertexVars vars(event.vertices[iV]);
          if(!vars.isGood) continue;
          double dZ(std::abs(vars.z - genVertex->Z()));
          if(effVars_.nVtx == 0) dZPrim = dZ;
          else if(dZ < dZPrim) effVars_.primVtx = false;
          ++effVars_.nVtx;
        }

        for(unsigned iG(0); iG != nG; ++iG){
          susy::Particle const& part(genParticles[iG]);
          if(part.status != 1 || !part.mother || std::abs(part.mother->pdgId) > 99) continue;
          if(part.momentum.Pt() < 20.) continue;
          unsigned absId(std::abs(part.pdgId));
          if(absId != 22 && absId != 11 && absId != 13) continue;
          if(absId == 11 || absId == 13){
            susy::Particle const* mother(part.mother);
            while(mother && mother->pdgId != 23 && std::abs(mother->pdgId) != 24) mother = mother->mother;
            if(!mother) continue;
          }

          double iso(0.);
          for(unsigned iI(0); iI != nG; ++iI){
            if(iI == iG) continue;
            susy::Particle const& isoP(genParticles[iI]);
            if(isoP.status != 1) continue;
            unsigned isoId(std::abs(isoP.pdgId));
            if(isoId == 12 || isoId == 14 || isoId == 16 || isoId == 1000022 || isoId == 1000039) continue;
            if(isoP.momentum.DeltaR(part.momentum) > 0.3) continue;
            iso += isoP.momentum.Pt();
          }

          switch(absId){
          case 22:
            genPhotons.push_back(std::make_pair(&part, iso));
            genMatchedPhotons.push_back(0);
            break;
          case 11:
            genElectrons.push_back(std::make_pair(&part, iso));
            genMatchedElectrons.push_back(0);
            break;
          case 13:
            genMuons.push_back(std::make_pair(&part, iso));
            genMatchedMuons.push_back(0);
            break;
          }
        }
      }

      // bug fix for tag cms533v0 / cms538v1
      std::vector<susy::PFParticle const*> pfParticles(susy::cleanPFParticles(event.pfParticles));
      unsigned nPF(pfParticles.size());

      /* SORT OBJECTS (ELECTRON SIZE WILL CHANGE) */

      std::vector<const susy::Photon*> photons(eventProducer_.sortPhotons(event.photons["photons"]));
      std::vector<const susy::Electron*> electrons(eventProducer_.sortElectrons(event.electrons["gsfElectrons"]));
      std::vector<const susy::Muon*> muons(eventProducer_.sortMuons(event.muons["muons"]));

      unsigned nPh(photons.size());
      unsigned nEl(electrons.size());
      unsigned nMu(muons.size());

      /* SELECT PHOTONS */

      unsigned nCandPhoton(0);
      unsigned nFakePhoton(0);
      unsigned nElePhoton(0);

      std::fill_n(ph_isCand_, susy::NMAX, false);
      std::fill_n(ph_isFake_, susy::NMAX, false);
      std::fill_n(ph_isEle_, susy::NMAX, false);

      int leadGoodPhoton(-1);

      for(unsigned iPh(0); iPh != nPh; ++iPh){
        susy::Photon const& ph(*photons[iPh]);

        if(ph.momentum.Pt() < 25.) break;
        if(std::abs(ph.caloPosition.Eta()) > susy::etaGapBegin) continue;

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

        unsigned iL(0);
        for(; iL != nMu; ++iL){
          if(muons[iL]->momentum.Pt() < 2.){
            iL = nMu;
            break;
          }
          if(muons[iL]->momentum.DeltaR(ph.momentum) < 0.3) break;
        }
        if(iL != nMu) continue;

        bool gsfVeto(true);
        for(iL = 0; iL != nEl; ++iL){
          if(electrons[iL]->momentum.Pt() < 2.){
            iL = nEl;
            break;
          }
          if(electrons[iL]->superClusterIndex == ph.superClusterIndex){
            gsfVeto = false;
            continue;
          }
          double dR(electrons[iL]->superCluster->position.DeltaR(ph.caloPosition));
          if(dR < 0.02) gsfVeto = false;
          else if(dR < 0.3) break;
        }
        if(iL != nEl) continue;

        susy::PhotonVars vars(ph, event);
        bool isGood(susy::ObjectSelector::isGoodPhoton(vars, susy::PhLoose12Pix, &phIdResults) && gsfVeto); // order matters! isGoodPhoton must execute

        if(isGood){
          if(leadGoodPhoton < 0) leadGoodPhoton = iPh;
          
          ph_isCand_[iPh] = true;
          ++nCandPhoton;

          if(effTree_){
            for(unsigned iG(0); iG != genPhotons.size(); ++iG){
              susy::Particle const& genPhoton(*genPhotons[iG].first);
              if(genPhoton.momentum.Vect().DeltaR(ph.caloPosition - genPhoton.vertex) < 0.1)
                genMatchedPhotons[iG] = &ph;
            }
          }
        }
        else if((phIdResults & phNoVeto) == phNoVeto){
          ph_isEle_[iPh] = true;
          ++nElePhoton;
        }
        else if(gsfVeto && phIdResults[susy::PhElectronVeto]){
          ph_isFake_[iPh] = true;
          ++nFakePhoton;
        }
      }

      if(nCandPhoton == 0 && nElePhoton == 0 && nFakePhoton == 0) continue;

      if(leadGoodPhoton >= 0 && photons[leadGoodPhoton]->momentum.Pt() > 40.){
        if(cutflowE == kGoodPhoton - 1){
          std::vector<susy::TriggerObjectCollection> photonHLTObjects[nElectronHLT];
          for(unsigned iHLT(0); iHLT != nElectronHLT; ++iHLT){
            if(!passElectronHLT[iHLT]) continue; // to save time
            for(unsigned iF(0); iF != photonHLTFiltersE[iHLT].size(); ++iF)
              photonHLTObjects[iHLT].push_back(triggerEvent.getFilterObjects(photonHLTFiltersE[iHLT][iF]));
          }

          unsigned iPh(0);
          for(; iPh != nPh; ++iPh){
            if(!ph_isCand_[iPh]) continue;
            if(photons[iPh]->momentum.Pt() < 40.){
              iPh = nPh;
              break;
            }

            TLorentzVector dSC(photons[iPh]->superCluster->position, 0.);

            unsigned iHLT(0);
            for(; iHLT != nElectronHLT; ++iHLT){
              if(!passElectronHLT[iHLT]) continue;
              unsigned iF(0);
              for(; iF != photonHLTFiltersE[iHLT].size(); ++iF){
                susy::TriggerObjectCollection& objects(photonHLTObjects[iHLT][iF]);
                unsigned iO(0);
                for(; iO != objects.size(); ++iO)
                  if(objects[iO].deltaR(dSC) < 0.15) break; // object matched
                if(iO == objects.size()) break; // no matching object for the filter
              }
              if(iF == photonHLTFiltersE[iHLT].size()) break; // all filters of the path matched
            }
            if(iHLT != nElectronHLT) break; // at least one path had its objects matching
          }

          if(iPh != nPh) ++cutflowE;
        }

        if(cutflowM == kGoodPhoton - 1){
          std::vector<susy::TriggerObjectCollection> photonHLTObjects[nMuonHLT];
          for(unsigned iHLT(0); iHLT != nMuonHLT; ++iHLT){
            if(!passMuonHLT[iHLT]) continue; // to save time
            for(unsigned iF(0); iF != photonHLTFiltersM[iHLT].size(); ++iF)
              photonHLTObjects[iHLT].push_back(triggerEvent.getFilterObjects(photonHLTFiltersM[iHLT][iF]));
          }

          unsigned iPh(0);
          for(; iPh != nPh; ++iPh){
            if(!ph_isCand_[iPh]) continue;
            if(photons[iPh]->momentum.Pt() < 40.){
              iPh = nPh;
              break;
            }

            TLorentzVector dSC(photons[iPh]->superCluster->position, 0.);

            unsigned iHLT(0);
            for(; iHLT != nMuonHLT; ++iHLT){
              if(!passMuonHLT[iHLT]) continue;
              unsigned iF(0);
              for(; iF != photonHLTFiltersM[iHLT].size(); ++iF){
                susy::TriggerObjectCollection& objects(photonHLTObjects[iHLT][iF]);
                unsigned iO(0);
                for(; iO != objects.size(); ++iO)
                  if(objects[iO].deltaR(dSC) < 0.15) break; // object matched
                if(iO == objects.size()) break; // no matching object for the filter
              }
              if(iF == photonHLTFiltersM[iHLT].size()) break; // all filters of the path matched
            }
            if(iHLT != nMuonHLT) break; // at least one path had its objects matching
          }

          if(iPh != nPh) ++cutflowM;
        }
      }

      /* SELECT ELECTRONS */

      unsigned nCandElectron(0);
      unsigned nFakeElectron(0);

      std::fill_n(el_isCand_, susy::NMAX, false);
      std::fill_n(el_isFake_, susy::NMAX, false);

      for(unsigned iEl(0); iEl != nEl; ++iEl){
        susy::Electron const& el(*electrons[iEl]);

        if(el.momentum.Pt() < 25.) break;

        susy::ElectronVars vars(el, event);

        if(vars.iSubdet == -1) continue;

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

          if(effTree_){
            for(unsigned iG(0); iG != genElectrons.size(); ++iG){
              susy::Particle const& genElectron(*genElectrons[iG].first);
              if(genElectron.momentum.DeltaR(el.momentum) < 0.05)
                genMatchedElectrons[iG] = &el;
            }
          }
        }
        else{
          el_isFake_[iEl] = true;
          ++nFakeElectron;
        }
      }

      if(nCandElectron != 0 && cutflowE == kGoodLepton - 1){
        std::vector<susy::TriggerObjectCollection> electronHLTObjects[nElectronHLT];
        for(unsigned iHLT(0); iHLT != nElectronHLT; ++iHLT){
          if(!passElectronHLT[iHLT]) continue; // to save time
          for(unsigned iF(0); iF != electronHLTFilters[iHLT].size(); ++iF)
            electronHLTObjects[iHLT].push_back(triggerEvent.getFilterObjects(electronHLTFilters[iHLT][iF]));
        }

        unsigned iEl(0);
        for(; iEl != nEl; ++iEl){
          if(!el_isCand_[iEl]) continue;
          if(electrons[iEl]->momentum.Pt() < 25.){
            iEl = nEl;
            break;
          }

          TLorentzVector dSC(electrons[iEl]->superCluster->position, 0.);
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
          if(iHLT != nElectronHLT) break; // at least one path had its objects matching
        }

        if(iEl != nEl) ++cutflowE;
      }

      /* SELECT MUONS */

      unsigned nCandMuon(0);
      unsigned nFakeMuon(0);

      std::fill_n(mu_isCand_, susy::NMAX, false);
      std::fill_n(mu_isFake_, susy::NMAX, false);

      for(unsigned iMu(0); iMu != nMu; ++iMu){
        susy::Muon const& mu(*muons[iMu]);

        if(mu.momentum.Pt() < 25.) break;

        susy::MuonVars vars(mu, event);

        if(vars.iSubdet == -1 || !vars.isLoose) continue;

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

          if(effTree_){
            for(unsigned iG(0); iG != genMuons.size(); ++iG){
              susy::Particle const& genMuon(*genMuons[iG].first);
              if(genMuon.momentum.DeltaR(mu.momentum) < 0.05)
                genMatchedMuons[iG] = &mu;
            }
          }
        }
        else{
          mu_isFake_[iMu] = true;
          ++nFakeMuon;
        }
      }

      if(nCandMuon != 0 && cutflowM == kGoodLepton - 1){
        std::vector<susy::TriggerObjectCollection> muonHLTObjects[nMuonHLT];
        for(unsigned iHLT(0); iHLT != nMuonHLT; ++iHLT){
          if(!passMuonHLT[iHLT]) continue; // to save time
          for(unsigned iF(0); iF != muonHLTFilters[iHLT].size(); ++iF)
            muonHLTObjects[iHLT].push_back(triggerEvent.getFilterObjects(muonHLTFilters[iHLT][iF]));
        }

        unsigned iMu(0);
        for(; iMu != nMu; ++iMu){
          if(!mu_isCand_[iMu]) continue;
          if(muons[iMu]->momentum.Pt() < 25.){
            iMu = nMu;
            break;
          }

          unsigned iHLT(0);
          for(; iHLT != nMuonHLT; ++iHLT){
            if(!passMuonHLT[iHLT]) continue;
            unsigned iF(0);
            for(; iF != muonHLTFilters[iHLT].size(); ++iF){
              susy::TriggerObjectCollection& objects(muonHLTObjects[iHLT][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(muons[iMu]->momentum) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == muonHLTFilters[iHLT].size()) break; // all filters of the path had a match
          }
          if(iHLT != nMuonHLT) break; // at least one path had its objects matching
        }

        if(iMu != nMu) ++cutflowM;
      }

      /* DETERMINE THE RESULT OF EACH FILTER */

      filterResults_[kPhotonAndElectron] = nCandPhoton != 0 && nCandElectron != 0;
      filterResults_[kPhotonAndMuon] = nCandPhoton != 0 && nCandMuon != 0;
      filterResults_[kElePhotonAndElectron] = nElePhoton != 0 && nCandElectron != 0;
      filterResults_[kElePhotonAndMuon] = nElePhoton != 0 && nCandMuon != 0;
      filterResults_[kFakePhotonAndElectron] = nFakePhoton != 0 && nCandElectron != 0;
      filterResults_[kFakePhotonAndMuon] = nFakePhoton != 0 && nCandMuon != 0;
      filterResults_[kPhotonAndFakeElectron] = nCandPhoton != 0 && nFakeElectron != 0;
      filterResults_[kPhotonAndFakeMuon] = nCandPhoton != 0 && nFakeMuon != 0;
      filterResults_[kElePhotonAndFakeElectron] = nElePhoton != 0 && nFakeElectron != 0;
      filterResults_[kElePhotonAndFakeMuon] = nElePhoton != 0 && nFakeMuon != 0;
      filterResults_[kFakePhotonAndFakeElectron] = nFakePhoton != 0 && nFakeElectron != 0;
      filterResults_[kFakePhotonAndFakeMuon] = nFakePhoton != 0 && nFakeMuon != 0;

      if(nElePhoton == 1 && nCandElectron == 1){
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

      unsigned iF(0);
      for(; iF != nFilterTypes; ++iF)
        if(filterResults_[iF] && useEvents_[iF]) break;

      if(iF == nFilterTypes) continue;

      /* EVENT PASSES AT LEAST ONE FILTER - LOAD FULL EVENT AND FILL OUTPUT */

      fullEvent.getEntry(iEntry);

      std::vector<const susy::PFJet*> jets(eventProducer_.sortJets(fullEvent.pfJets["ak5"]));
      unsigned nJ(jets.size());

      for(unsigned iJ(0); iJ != nJ; ++iJ){
        jt_isCand_[iJ] = false;

        susy::PFJet const& jet(*jets[iJ]);

        susy::JetVars vars(jet, event);
        if(!vars.isLoose || !vars.passPUJetIdLoose) continue;

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

      if(effTree_){
        for(unsigned iG(0); iG != genPhotons.size(); ++iG){
          effVars_.pdgId = 22;
          effVars_.genPt = genPhotons[iG].first->momentum.Pt();
          effVars_.genIso = genPhotons[iG].second;
          susy::Photon const* reco(genMatchedPhotons[iG]);
          if(reco){
            effVars_.pass = true;
            effVars_.pt = reco->momentum.Pt();
            effVars_.eta = reco->caloPosition.Eta();
            effVars_.phi = reco->caloPosition.Phi();
          }
          else{
            effVars_.pass = false;
            effVars_.pt = 0.;
            effVars_.eta = 0.;
            effVars_.phi = 0.;
          }
          
          effTree_->Fill();
        }

        for(unsigned iG(0); iG != genElectrons.size(); ++iG){
          effVars_.pdgId = genElectrons[iG].first->pdgId;
          effVars_.genPt = genElectrons[iG].first->momentum.Pt();
          effVars_.genIso = genElectrons[iG].second;
          susy::Electron const* reco(genMatchedElectrons[iG]);
          if(reco){
            effVars_.pass = true;
            effVars_.pt = reco->momentum.Pt();
            effVars_.eta = reco->superCluster->position.Eta();
            effVars_.phi = reco->superCluster->position.Phi();
          }
          else{
            effVars_.pass = false;
            effVars_.pt = 0.;
            effVars_.eta = 0.;
            effVars_.phi = 0.;
          }

          effTree_->Fill();
        }

        for(unsigned iG(0); iG != genMuons.size(); ++iG){
          effVars_.pdgId = genMuons[iG].first->pdgId;
          effVars_.genPt = genMuons[iG].first->momentum.Pt();
          effVars_.genIso = genMuons[iG].second;
          susy::Muon const* reco(genMatchedMuons[iG]);
          if(reco){
            effVars_.pass = true;
            effVars_.pt = reco->momentum.Pt();
            effVars_.eta = reco->momentum.Eta();
            effVars_.phi = reco->momentum.Phi();
          }
          else{
            effVars_.pass = false;
            effVars_.pt = 0.;
            effVars_.eta = 0.;
            effVars_.phi = 0.;
          }

          effTree_->Fill();
        }
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
  outputFile_->cd();
  evtTree_->Write();
  allObjTree_->Write();
  cutTree_->Write();
  if(effTree_) effTree_->Write();

  delete outputFile_;

  outputFile_ = 0;
  evtTree_ = 0;
  allObjTree_ = 0;
  cutTree_ = 0;
  effTree_ = 0;

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
