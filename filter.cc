#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TEntryList.h"

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
  kSoftPhotonAndElectron,
  kSoftPhotonAndMuon,
  kSoftElePhotonAndElectron,
  kSoftElePhotonAndMuon,
  kSoftFakePhotonAndElectron,
  kSoftFakePhotonAndMuon,
  kSoftPhotonAndFakeElectron,
  kSoftPhotonAndFakeMuon,
  nFilterTypes,
  nHardPhotonFilters = kSoftPhotonAndElectron
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
  "FakePhotonAndFakeMuon",
  "SoftPhotonAndElectron",
  "SoftPhotonAndMuon",
  "SoftElePhotonAndElectron",
  "SoftElePhotonAndMuon",
  "SoftFakePhotonAndElectron",
  "SoftFakePhotonAndMuon",
  "SoftPhotonAndFakeElectron",
  "SoftPhotonAndFakeMuon"
};

enum Cut {
  kAllEvents,
  kGoodLumi,
  kHLTE,
  kHLTM,
  kHLTESingle,
  kHLTMSingle,
  kMetFilter,
  kGoodVertex,
  kGoodPhotonE,
  kGoodPhotonM,
  kGoodSoftPhoton,
  kGoodElectron,
  kGoodHardElectron,
  kGoodMuon,
  kGoodHardMuon,
  nCuts
};

TString cuts[nCuts] = {
  "AllEvents",
  "GoodLumi",
  "HLTE",
  "HLTM",
  "HLTESingle",
  "HLTMSingle",
  "MetFilter",
  "GoodVertex",
  "GoodPhotonE",
  "GoodPhotonM",
  "GoodSoftPhoton",
  "GoodElectron",
  "GoodHardElectron",
  "GoodMuon",
  "GoodHardMuon"
};

class PhotonLeptonFilter : public SimpleTreeProducer {
public:
  PhotonLeptonFilter();
  ~PhotonLeptonFilter();
  bool initialize(char const*, char const*, bool = false);
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

  TEntryList* hardPhotonList_;
  TEntryList* softPhotonList_;
  TTree* effTree_;
  TTree* cutTree_;

  std::bitset<nFilterTypes> useEvents_;

  susy::GoodLumis goodLumis_;

  bool filterResults_[nFilterTypes];

  bool ph_isCand_[susy::NMAX];
  bool ph_isFake_[susy::NMAX];
  bool ph_isEle_[susy::NMAX];
  bool ph_isSoftCand_[susy::NMAX];
  bool ph_isSoftFake_[susy::NMAX];
  bool ph_isSoftEle_[susy::NMAX];
  bool el_isCand_[susy::NMAX];
  bool el_isFake_[susy::NMAX];
  bool mu_isCand_[susy::NMAX];
  bool mu_isFake_[susy::NMAX];
  bool jt_isCand_[susy::NMAX];

  struct EffTreeVariables {
    //mkbranch
    char pdgId;
    bool reco;
    bool pass;
    float pt;
    float eta;
    float phi;
    float genPt;
    float genEta;
    float genPhi;
    float genIso;
    unsigned char nVtx;
    float puWeight;
    bool primVtx;
    //mkbranch
  } effVars_;

  struct GenMatch {
    susy::Particle const* particle;
    susy::Photon const* photon;
    susy::Electron const* electron;
    susy::Muon const* muon;
    float genIso;
    bool pass;
    GenMatch(susy::Particle const& _part, float _iso) : particle(&_part), photon(0), electron(0), muon(0), genIso(_iso), pass(false) {}
    GenMatch(GenMatch const& _orig) : particle(_orig.particle), photon(_orig.photon), electron(_orig.electron), muon(_orig.muon), genIso(_orig.genIso), pass(_orig.pass) {}
    GenMatch& operator=(GenMatch const& _rhs)
    {
      particle = _rhs.particle; photon = _rhs.photon; electron = _rhs.electron; muon = _rhs.muon; genIso = _rhs.genIso; pass = _rhs.pass;
      return *this;
    }
  };

  bool storeRadiation_;
  float radPt_;
  bool isFSR_;

  GenDecayFilterRA3* genFilter_;
};

PhotonLeptonFilter::PhotonLeptonFilter() :
  SimpleTreeProducer(),
  input_("susyTree"),
  triggerInput_("triggerEventTree"),
  fullInput_("susyTree"),
  fullTriggerInput_("triggerEventTree"),
  hardPhotonList_(0),
  softPhotonList_(0),
  effTree_(0),
  cutTree_(0),
  goodLumis_(),
  storeRadiation_(false),
  radPt_(0.),
  isFSR_(false),
  genFilter_(0)
{
}

PhotonLeptonFilter::~PhotonLeptonFilter()
{
  delete hardPhotonList_;
  delete softPhotonList_;
  delete effTree_;
  delete cutTree_;
  delete genFilter_;
}

bool
PhotonLeptonFilter::initialize(char const* _outputDir, char const* _configFileNames, bool _storeRadiation/* = false*/)
{
  /* CONFIGURE */

  std::map<TString, TString> configRecords;

  TPRegexp configPat("^[ ]*([A-Z0-9_]+)[ ]*([+]?=)[ ]*([^ ].*)$");
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
      if(matches->GetEntries() == 4){
        std::cout << matches->At(1)->GetName() << " " << matches->At(2)->GetName() << " " << matches->At(3)->GetName() << std::endl;
        if(TString(matches->At(2)->GetName()) == "=")
          configRecords[matches->At(1)->GetName()] = matches->At(3)->GetName();
        else{
          configRecords[matches->At(1)->GetName()] += " ";
          configRecords[matches->At(1)->GetName()] += matches->At(3)->GetName();
        }
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
      std::cerr << ("Filter type " + filter + " not defined") << std::endl;
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
  allObjTree_->Branch("photon.isSoftCand", ph_isSoftCand_, "isSoftCand[photon.size]/O");
  allObjTree_->Branch("photon.isSoftFake", ph_isSoftFake_, "isSoftFake[photon.size]/O");
  allObjTree_->Branch("photon.isSoftEle", ph_isSoftEle_, "isSoftEle[photon.size]/O");
  allObjTree_->Branch("electron.isCand", el_isCand_, "isCand[electron.size]/O");
  allObjTree_->Branch("electron.isFake", el_isFake_, "isFake[electron.size]/O");
  allObjTree_->Branch("muon.isCand", mu_isCand_, "isCand[muon.size]/O");
  allObjTree_->Branch("muon.isFake", mu_isFake_, "isFake[muon.size]/O");
  allObjTree_->Branch("jet.isCand", jt_isCand_, "isCand[jet.size]/O");

  for(unsigned iF(0); iF != nFilterTypes; ++iF)
    evtTree_->Branch(filterNames[iF], filterResults_ + iF, filterNames[iF] + "/O");

  /* SETUP RADIATION VETO */
  // Store highest Pt of gen-isolated photon
  // Non-isolated photon has a sibling (shares motherIndex) particle that is not a neutrino and has Pt > 2 GeV within dR < 0.3
  storeRadiation_ = _storeRadiation;

  if(storeRadiation_){
    evtTree_->Branch("radPt", &radPt_, "radPt/F");
    evtTree_->Branch("isFSR", &isFSR_, "isFSR/O");
  }

  outputFile_->cd();

  hardPhotonList_ = new TEntryList("hardPhotonList", "Hard photon events", evtTree_);
  softPhotonList_ = new TEntryList("softPhotonList", "Soft photon events", evtTree_);

  if(configRecords["PUSCENARIO"] != ""){
    effTree_ = new TTree("effTree", "Efficiency Tree");
    effTree_->Branch("pdgId", &effVars_.pdgId, "pdgId/B");
    effTree_->Branch("reco", &effVars_.reco, "reco/O");
    effTree_->Branch("pass", &effVars_.pass, "pass/O");
    effTree_->Branch("pt", &effVars_.pt, "pt/F");
    effTree_->Branch("eta", &effVars_.eta, "eta/F");
    effTree_->Branch("phi", &effVars_.phi, "phi/F");
    effTree_->Branch("genPt", &effVars_.genPt, "genPt/F");
    effTree_->Branch("genEta", &effVars_.genEta, "genEta/F");
    effTree_->Branch("genPhi", &effVars_.genPhi, "genPhi/F");
    effTree_->Branch("genIso", &effVars_.genIso, "genIso/F");
    effTree_->Branch("nVtx", &effVars_.nVtx, "nVtx/b");
    effTree_->Branch("puWeight", &effVars_.puWeight, "puWeight/F");
    effTree_->Branch("primVtx", &effVars_.primVtx, "primVtx/O");
  }

  cutTree_ = new TTree("cutTree", "Cutflow");
  cutTree_->Branch("run", 0, "run/i");
  cutTree_->Branch("lumi", 0, "lumi/i");
  cutTree_->Branch("event", 0, "event/i");
  for(unsigned iC(0); iC != nCuts; ++iC)
    cutTree_->Branch(cuts[iC], 0, cuts[iC] + "/O");

  TTree* cfgTree(new TTree("cfgTree", "Configuration parameters"));
  TString* cfgName(new TString);
  TString* cfgVal(new TString);
  cfgTree->Branch("name", "TString", &cfgName);
  cfgTree->Branch("value", "TString", &cfgVal);
  for(std::map<TString, TString>::iterator cItr(configRecords.begin()); cItr != configRecords.end(); ++cItr){
    if(cItr->second == "") continue;
    *cfgName = cItr->first;
    *cfgVal = cItr->second;
    cfgTree->Fill();
  }
  cfgTree->Write();

  delete cfgTree;
  delete cfgName;
  delete cfgVal;

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
  susy::ObjectTree::setBranchStatus(input_, true, true, true, false, true); // Photon, Electron, Muon, Jet, Vertex
  if(USEBEAMSPOTIP) input_.SetBranchStatus("beamSpot*", 1);
  if(eventProducer_.isMC()){
    input_.SetBranchStatus("genParticles*", 1);
    input_.SetBranchStatus("pu*", 1);
  }

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
  bool cutflow[nCuts] = {true};
  for(unsigned iC(0); iC != nCuts; ++iC)
    cutTree_->SetBranchAddress(cuts[iC], cutflow + iC);

  /* TRIGGERS (ONLY FOR CUTFLOW - NO FILTERING DONE) */

  TString electronHLT[] = {
    "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50"
  };
  unsigned const nElectronHLT(sizeof(electronHLT) / sizeof(TString));

  TString singleElectronHLT[] = {
    "HLT_Ele27_WP80"
  };
  unsigned const nSingleElectronHLT(sizeof(singleElectronHLT) / sizeof(TString));

  TString muonHLT[] = {
    "HLT_Mu22_Photon22_CaloIdL"
  };
  unsigned const nMuonHLT(sizeof(muonHLT) / sizeof(TString));

  TString singleMuonHLT[] = {
    "HLT_IsoMu24",
    "HLT_IsoMu24_eta2p1",
    "HLT_IsoMu24_eta2p1" // filter implementation changed; repeat to get the OR of the different implementations
  };
  unsigned const nSingleMuonHLT(sizeof(singleMuonHLT) / sizeof(TString));

  std::vector<TString> electronHLTFilters[nElectronHLT];
  electronHLTFilters[0].push_back("hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded");

  std::vector<TString> singleElectronHLTFilters[nSingleElectronHLT];
  singleElectronHLTFilters[0].push_back("hltEle27WP80TrackIsoFilter");

  std::vector<TString> muonHLTFilters[nMuonHLT];
  muonHLTFilters[0].push_back("hltL1Mu3p5EG12L3Filtered22");

  std::vector<TString> singleMuonHLTFilters[nSingleMuonHLT];
  singleMuonHLTFilters[0].push_back("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15");
  singleMuonHLTFilters[1].push_back("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10");
  singleMuonHLTFilters[2].push_back("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15");

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
  std::bitset<susy::nMuonCriteria> muBaseline(susy::ObjectSelector::muReferences[susy::MuTight12]);
  muBaseline.reset(susy::MuCombIso);
  muBaseline.reset(susy::MuDxy);
  muBaseline.reset(susy::MuDz);
  std::bitset<susy::nElectronCriteria> elIdResults;
  std::bitset<susy::nElectronCriteria> elNoIP(susy::ObjectSelector::elReferences[susy::ElMedium12]);
  elNoIP.reset(susy::ElD0);
  elNoIP.reset(susy::ElDZ);
  std::bitset<susy::nElectronCriteria> elBaseline(susy::ObjectSelector::elReferences[susy::ElMedium12]);
  elBaseline.reset(susy::ElCombIso);
  elBaseline.reset(susy::ElDeltaEta);
  elBaseline.reset(susy::ElDeltaPhi);
  elBaseline.reset(susy::ElD0);
  elBaseline.reset(susy::ElDZ);
  std::bitset<susy::nPhotonCriteria> phIdResults;
  std::bitset<susy::nPhotonCriteria> phNoVeto(susy::ObjectSelector::phReferences[susy::PhLoose12Pix]);
  phNoVeto.reset(susy::PhElectronVeto);

  /* EVENT TYPE MASKS */

  bool useHardPhoton(false);
  bool useSoftPhoton(false);
  for(unsigned iF(0); iF != nHardPhotonFilters; ++iF)
    if(useEvents_[iF]) useHardPhoton = true;
  for(unsigned iF(nHardPhotonFilters); iF != nFilterTypes; ++iF)
    if(useEvents_[iF]) useSoftPhoton = true;

  /* START LOOP */
  
  int nRead(0);
  for(long iEntry(0); (nRead = event.getEntry(iEntry)) != 0; cutTree_->Fill(), ++iEntry){
    try{
      if(nRead < 0)
        throw susy::Exception(susy::Exception::kIOError, "Corrupt input");
    
      if(iEntry % 1000 == 0) std::cout << "Analyzing event " << iEntry << std::endl;

      if(genFilter_ && !genFilter_->pass(event)) continue;

      std::fill_n(cutflow + 1, nCuts - 1, false);

      if(!goodLumis_.isGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;

      cutflow[kGoodLumi] = true;

      std::bitset<nElectronHLT> passElectronHLT;
      std::bitset<nMuonHLT> passMuonHLT;
      std::bitset<nSingleElectronHLT> passSingleElectronHLT;
      std::bitset<nSingleMuonHLT> passSingleMuonHLT;

      for(unsigned iHLT(0); iHLT != nElectronHLT; ++iHLT)
        passElectronHLT[iHLT] = event.hltMap.pass(electronHLT[iHLT] + "_v*");
      for(unsigned iHLT(0); iHLT != nSingleElectronHLT; ++iHLT)
        passSingleElectronHLT[iHLT] = event.hltMap.pass(singleElectronHLT[iHLT] + "_v*");
      for(unsigned iHLT(0); iHLT != nMuonHLT; ++iHLT)
        passMuonHLT[iHLT] = event.hltMap.pass(muonHLT[iHLT] + "_v*");
      for(unsigned iHLT(0); iHLT != nSingleMuonHLT; ++iHLT)
        passSingleMuonHLT[iHLT] = event.hltMap.pass(singleMuonHLT[iHLT] + "_v*");

      if(passElectronHLT.any()) cutflow[kHLTE] = true;
      if(passMuonHLT.any()) cutflow[kHLTM] = true;
      if(passSingleElectronHLT.any()) cutflow[kHLTESingle] = true;
      if(passSingleMuonHLT.any()) cutflow[kHLTMSingle] = true;

      if(!event.passMetFilters()) continue;

      cutflow[kMetFilter] = true;

      unsigned nV(event.vertices.size());
      unsigned iV(0);
      for(; iV != nV; ++iV){
        susy::VertexVars vars(event.vertices[iV]);
        if(vars.isGood) break;
      }
      if(iV == nV) continue;

      cutflow[kGoodVertex] = true;

      if(storeRadiation_){
        susy::ParticleCollection const& genParticles(event.genParticles);
        unsigned nG(genParticles.size());

        radPt_ = 0.;
        isFSR_ = false;

        for(unsigned iG(0); iG != nG; ++iG){
          susy::Particle const& particle(genParticles[iG]);
          if(particle.status != 1 || particle.pdgId != 22) continue;
          double pt(particle.momentum.Pt());
          if(pt < 2. || pt < radPt_) continue;

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

          radPt_ = pt;

          short idx(particle.motherIndex);
          while(idx != -1 && genParticles[idx].pdgId != 23 && std::abs(genParticles[idx].pdgId) != 24) idx = genParticles[idx].motherIndex;
          isFSR_ = (idx != -1);
        }
      }

      std::vector<GenMatch> genMatches;

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
        
        effVars_.puWeight = 1.;
        for(unsigned iPU(0); iPU != event.pu.size(); ++iPU){
          if(event.pu[iPU].BX != 0) continue;
          effVars_.puWeight = eventProducer_.getPUWeight(event.pu[iPU].trueNumInteractions);
          break;
        }

        for(unsigned iG(0); iG != nG; ++iG){
          susy::Particle const& part(genParticles[iG]);
          if(part.status != 1 || !part.mother || std::abs(part.mother->pdgId) > 99) continue;
          double absEta(std::abs(part.momentum.Eta()));
          if(part.momentum.Pt() < 20. || absEta > 3.) continue;
          unsigned absId(std::abs(part.pdgId));
          if(absId != 22 && absId != 11 && absId != 13) continue;
          if(absId == 11 || absId == 13){
            susy::Particle const* mother(part.mother);
            while(mother && mother->pdgId != 23 && std::abs(mother->pdgId) != 24) mother = mother->mother;
            if(!mother) continue;
          }
          else if(absEta > 1.6) continue;

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

          genMatches.push_back(GenMatch(part, iso));
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
      unsigned nSoftCandPhoton(0);
      unsigned nSoftFakePhoton(0);
      unsigned nSoftElePhoton(0);

      std::fill_n(ph_isCand_, susy::NMAX, false);
      std::fill_n(ph_isFake_, susy::NMAX, false);
      std::fill_n(ph_isEle_, susy::NMAX, false);
      std::fill_n(ph_isSoftCand_, susy::NMAX, false);
      std::fill_n(ph_isSoftFake_, susy::NMAX, false);
      std::fill_n(ph_isSoftEle_, susy::NMAX, false);

      int leadGoodPhoton(-1);

      for(unsigned iPh(0); iPh != nPh; ++iPh){
        susy::Photon const& ph(*photons[iPh]);

        unsigned iGen(0);
        for(; iGen != genMatches.size(); ++iGen){
          GenMatch& gen(genMatches[iGen]);
          if(gen.particle->pdgId != 22) continue;
          if(gen.particle->momentum.Vect().DeltaR(ph.caloPosition - gen.particle->vertex) < 0.1){
            gen.photon = &ph;
            break;
          }
        }

        if(ph.hadTowOverEm > 0.05) continue;

        susy::PhotonVars vars(ph, event);

        if(useSoftPhoton){
          if(((vars.iSubdet == 0 && ph.sigmaIetaIeta > 0.005) ||
              (vars.iSubdet == 1 && std::abs(ph.caloPosition.Eta()) < 2.5 && ph.sigmaIetaIeta > 0.019))){
            bool isPhoton(susy::ObjectSelector::isGoodPhoton(vars, susy::PhMedium12, &phIdResults));

            if(isPhoton){
              ph_isSoftCand_[iPh] = true;
              ++nSoftCandPhoton;
            }
            else if((phIdResults & phNoVeto) == phNoVeto){
              ph_isSoftEle_[iPh] = true;
              ++nSoftElePhoton;
            }
            else if(phIdResults[susy::PhElectronVeto] && vars.chargedHadronIso < 15. && vars.neutralHadronIso < 10. && vars.photonIso < 10.){
              ph_isSoftFake_[iPh] = true;
              ++nSoftFakePhoton;
            }
          }
        }

        if(useHardPhoton){
          if(ph.momentum.Pt() > 25. && std::abs(ph.caloPosition.Eta()) < susy::etaGapBegin){

            unsigned iPF(0);
            for(; iPF != nPF; ++iPF){
              susy::PFParticle const& part(*pfParticles[iPF]);
              if(part.momentum.Pt() < 3.) continue;
              TVector3 dir(ph.caloPosition - part.vertex);
              if(std::abs(part.pdgId) == 211 &&
                 std::abs(dir.Eta() - part.momentum.Eta()) < 0.005 &&
                 std::abs(TVector2::Phi_mpi_pi(dir.Phi() - part.momentum.Phi())) < 0.02) break;
            }

            if(iPF == nPF){

              unsigned iL(0);
              for(; iL != nMu; ++iL){
                if(muons[iL]->momentum.Pt() < 2.){
                  iL = nMu;
                  break;
                }
                if(muons[iL]->momentum.DeltaR(ph.momentum) < 0.3) break;
              }

              if(iL == nMu){

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

                if(iL == nEl){

                  bool isPhoton(susy::ObjectSelector::isGoodPhoton(vars, susy::PhLoose12Pix, &phIdResults) && gsfVeto); // order matters! isGoodPhoton must execute                  

                  if(isPhoton){
                    if(leadGoodPhoton < 0) leadGoodPhoton = iPh;
          
                    ph_isCand_[iPh] = true;
                    ++nCandPhoton;

                    if(iGen != genMatches.size()) genMatches[iGen].pass = true;
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
              }
            }
          }
        }
      }

      if(leadGoodPhoton >= 0 && photons[leadGoodPhoton]->momentum.Pt() > 40.){
        std::vector<susy::TriggerObjectCollection> photonHLTObjectsE[nElectronHLT];
        for(unsigned iHLT(0); iHLT != nElectronHLT; ++iHLT){
          if(!passElectronHLT[iHLT]) continue; // to save time
          for(unsigned iF(0); iF != photonHLTFiltersE[iHLT].size(); ++iF)
            photonHLTObjectsE[iHLT].push_back(triggerEvent.getFilterObjects(photonHLTFiltersE[iHLT][iF]));
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
              susy::TriggerObjectCollection& objects(photonHLTObjectsE[iHLT][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(dSC) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == photonHLTFiltersE[iHLT].size()) break; // all filters of the path matched
          }
          if(iHLT != nElectronHLT) break; // at least one path had its objects matching
        }

        if(iPh != nPh) cutflow[kGoodPhotonE] = true;

        std::vector<susy::TriggerObjectCollection> photonHLTObjectsM[nMuonHLT];
        for(unsigned iHLT(0); iHLT != nMuonHLT; ++iHLT){
          if(!passMuonHLT[iHLT]) continue; // to save time
          for(unsigned iF(0); iF != photonHLTFiltersM[iHLT].size(); ++iF)
            photonHLTObjectsM[iHLT].push_back(triggerEvent.getFilterObjects(photonHLTFiltersM[iHLT][iF]));
        }

        iPh = 0;
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
              susy::TriggerObjectCollection& objects(photonHLTObjectsM[iHLT][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(dSC) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == photonHLTFiltersM[iHLT].size()) break; // all filters of the path matched
          }
          if(iHLT != nMuonHLT) break; // at least one path had its objects matching
        }

        if(iPh != nPh) cutflow[kGoodPhotonM] = true;
      }

      if(nSoftCandPhoton != 0) cutflow[kGoodSoftPhoton] = true;

      /* SELECT ELECTRONS */

      unsigned nCandElectron(0);
      unsigned nFakeElectron(0);
      unsigned nVetoElectron(0);

      std::fill_n(el_isCand_, susy::NMAX, false);
      std::fill_n(el_isFake_, susy::NMAX, false);

      int leadGoodElectron(-1);

      for(unsigned iEl(0); iEl != nEl; ++iEl){
        susy::Electron const& el(*electrons[iEl]);

        if(el.momentum.Pt() < 10.) break;

        unsigned iGen(0);
        for(; iGen != genMatches.size(); ++iGen){
          GenMatch& gen(genMatches[iGen]);
          if(std::abs(gen.particle->pdgId) != 11) continue;
          if(gen.particle->momentum.DeltaR(el.momentum) < 0.05){
            gen.electron = &el;
            break;
          }
        }

        susy::ElectronVars vars(el, event);

        if(vars.isVeto) ++nVetoElectron;

        if(vars.pt < 25. || vars.iSubdet == -1) continue;

        susy::ObjectSelector::isGoodElectron(vars, susy::ElMedium12, &elIdResults);
        if((elIdResults & elBaseline) != elBaseline) continue;

        bool isMedium(false);
        if(USEBEAMSPOTIP){
          if((elIdResults & elNoIP) == elNoIP)
            isMedium = std::abs(el.gsfTrack->dxy(event.beamSpot)) < 0.2;
        }
        else
          isMedium = vars.isMedium;

        if(isMedium){
          if(leadGoodElectron < 0) leadGoodElectron = iEl;

          el_isCand_[iEl] = true;
          ++nCandElectron;

          if(iGen != genMatches.size()) genMatches[iGen].pass = true;
        }
        else{
          el_isFake_[iEl] = true;
          ++nFakeElectron;
        }
      }

      if(nCandElectron != 0){
        if(useHardPhoton){
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

          if(iEl != nEl) cutflow[kGoodElectron] = true;
        }

        if(useSoftPhoton && nVetoElectron == 1 && electrons[leadGoodElectron]->momentum.Pt() > 30.){
          std::vector<susy::TriggerObjectCollection> singleElectronHLTObjects[nSingleElectronHLT];
          for(unsigned iHLT(0); iHLT != nSingleElectronHLT; ++iHLT){
            if(!passSingleElectronHLT[iHLT]) continue; // to save time
            for(unsigned iF(0); iF != singleElectronHLTFilters[iHLT].size(); ++iF)
              singleElectronHLTObjects[iHLT].push_back(triggerEvent.getFilterObjects(singleElectronHLTFilters[iHLT][iF]));
          }

          TLorentzVector dSC(electrons[leadGoodElectron]->superCluster->position, 0.);
          unsigned iHLT(0);
          for(; iHLT != nSingleElectronHLT; ++iHLT){
            if(!passSingleElectronHLT[iHLT]) continue;
            unsigned iF(0);
            for(; iF != singleElectronHLTFilters[iHLT].size(); ++iF){
              susy::TriggerObjectCollection& objects(singleElectronHLTObjects[iHLT][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(dSC) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == singleElectronHLTFilters[iHLT].size()) break; // all filters of the path had a match
          }

          if(iHLT != nSingleElectronHLT) cutflow[kGoodHardElectron] = true; // at least one path had its objects matching
        }
      }

      /* SELECT MUONS */

      unsigned nCandMuon(0);
      unsigned nFakeMuon(0);
      unsigned nVetoMuon(0);

      std::fill_n(mu_isCand_, susy::NMAX, false);
      std::fill_n(mu_isFake_, susy::NMAX, false);

      int leadGoodMuon(-1);

      for(unsigned iMu(0); iMu != nMu; ++iMu){
        susy::Muon const& mu(*muons[iMu]);

        double pt(mu.momentum.Pt());

        if(pt < 10.) break;

        unsigned iGen(0);
        for(; iGen != genMatches.size(); ++iGen){
          GenMatch& gen(genMatches[iGen]);
          if(std::abs(gen.particle->pdgId) != 13) continue;
          if(gen.particle->momentum.DeltaR(mu.momentum) < 0.05){
            gen.muon = &mu;
            break;
          }
        }

        if(std::abs(mu.momentum.Eta()) < 2.4) ++nVetoMuon;

        if(pt < 25.) continue;

        susy::MuonVars vars(mu, event);

        if(vars.iSubdet == -1 || !vars.isLoose || vars.combRelIso > 1.) continue;

        susy::ObjectSelector::isGoodMuon(vars, susy::MuTight12, &muIdResults);
        if((muIdResults & muBaseline) != muBaseline) continue;

        bool isTight(false);
        if(USEBEAMSPOTIP){
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
          if(leadGoodMuon < 0) leadGoodMuon = iMu;

          mu_isCand_[iMu] = true;
          ++nCandMuon;

          if(iGen != genMatches.size()) genMatches[iGen].pass = true;
        }
        else{
          mu_isFake_[iMu] = true;
          ++nFakeMuon;
        }
      }

      if(nCandMuon != 0){
        if(useHardPhoton){
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

          if(iMu != nMu) cutflow[kGoodMuon] = true;
        }

        if(useSoftPhoton && nVetoMuon == 1 && muons[leadGoodMuon]->momentum.Pt() > 26. && std::abs(muons[leadGoodMuon]->momentum.Eta()) < 2.1){
          std::vector<susy::TriggerObjectCollection> singleMuonHLTObjects[nSingleMuonHLT];
          for(unsigned iHLT(0); iHLT != nSingleMuonHLT; ++iHLT){
            if(!passSingleMuonHLT[iHLT]) continue; // to save time
            for(unsigned iF(0); iF != singleMuonHLTFilters[iHLT].size(); ++iF)
              singleMuonHLTObjects[iHLT].push_back(triggerEvent.getFilterObjects(singleMuonHLTFilters[iHLT][iF]));
          }

          unsigned iHLT(0);
          for(; iHLT != nSingleMuonHLT; ++iHLT){
            if(!passSingleMuonHLT[iHLT]) continue;
            unsigned iF(0);
            for(; iF != singleMuonHLTFilters[iHLT].size(); ++iF){
              susy::TriggerObjectCollection& objects(singleMuonHLTObjects[iHLT][iF]);
              unsigned iO(0);
              for(; iO != objects.size(); ++iO)
                if(objects[iO].deltaR(muons[leadGoodMuon]->momentum) < 0.15) break; // object matched
              if(iO == objects.size()) break; // no matching object for the filter
            }
            if(iF == singleMuonHLTFilters[iHLT].size()) break; // all filters of the path had a match
          }

          if(iHLT != nSingleMuonHLT) cutflow[kGoodHardMuon] = true; // at least one path had its objects matching
        }
      }

      /* FILL GEN EFFICIENCY TREE */

      if(effTree_){
        for(unsigned iG(0); iG != genMatches.size(); ++iG){
          GenMatch& gen(genMatches[iG]);
          susy::Particle const& part(*gen.particle);
          effVars_.pdgId = part.pdgId;
          effVars_.genPt = part.momentum.Pt();
          effVars_.genEta = part.momentum.Eta();
          effVars_.genPhi = part.momentum.Phi();
          effVars_.genIso = gen.genIso;

          if(gen.photon){
            susy::Photon const& photon(*gen.photon);
            effVars_.reco = true;
            effVars_.pt = photon.momentum.Pt();
            effVars_.eta = photon.caloPosition.Eta();
            effVars_.phi = photon.caloPosition.Phi();
          }
          else if(gen.electron){
            susy::Electron const& electron(*gen.electron);
            effVars_.reco = true;
            effVars_.pt = electron.momentum.Pt();
            effVars_.eta = electron.superCluster->position.Eta();
            effVars_.phi = electron.superCluster->position.Phi();
          }
          else if(gen.muon){
            susy::Muon const& muon(*gen.muon);
            effVars_.reco = true;
            effVars_.pt = muon.momentum.Pt();
            effVars_.eta = muon.momentum.Eta();
            effVars_.phi = muon.momentum.Phi();
          }
          else{
            effVars_.reco = false;
            effVars_.pt = effVars_.eta = effVars_.phi = 0.;
          }

          effVars_.pass = gen.pass;

          effTree_->Fill();
        }
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

      filterResults_[kSoftPhotonAndElectron] = nSoftCandPhoton != 0 && nCandElectron != 0;
      filterResults_[kSoftPhotonAndMuon] = nSoftCandPhoton != 0 && nCandMuon != 0;
      filterResults_[kSoftElePhotonAndElectron] = nSoftElePhoton != 0 && nCandElectron != 0;
      filterResults_[kSoftElePhotonAndMuon] = nSoftElePhoton != 0 && nCandMuon != 0;
      filterResults_[kSoftFakePhotonAndElectron] = nSoftFakePhoton != 0 && nCandElectron != 0;
      filterResults_[kSoftFakePhotonAndMuon] = nSoftFakePhoton != 0 && nCandMuon != 0;
      filterResults_[kSoftPhotonAndFakeElectron] = nSoftCandPhoton != 0 && nFakeElectron != 0;
      filterResults_[kSoftPhotonAndFakeMuon] = nSoftCandPhoton != 0 && nFakeMuon != 0;

      if(nCandElectron == 1){
        unsigned iCandEl(0);
        for(; iCandEl != nEl; ++iCandEl)
          if(el_isCand_[iCandEl]) break;

        if(nElePhoton == 1){
          // is the only elePhoton actually the candidate electron?
          unsigned iElePh(0);
          for(; iElePh != nPh; ++iElePh)
            if(ph_isEle_[iElePh]) break;

          if(electrons[iCandEl]->superClusterIndex == photons[iElePh]->superClusterIndex ||
             electrons[iCandEl]->superCluster->position.DeltaR(photons[iElePh]->caloPosition) < 0.02)
            filterResults_[kElePhotonAndElectron] = false;
        }

        if(nSoftElePhoton == 1){
          // is the only elePhoton actually the candidate electron?
          unsigned iElePh(0);
          for(; iElePh != nPh; ++iElePh)
            if(ph_isSoftEle_[iElePh]) break;

          if(electrons[iCandEl]->superClusterIndex == photons[iElePh]->superClusterIndex ||
             electrons[iCandEl]->superCluster->position.DeltaR(photons[iElePh]->caloPosition) < 0.02)
            filterResults_[kSoftElePhotonAndElectron] = false;
        }
      }

      /* DID THE EVENT PASS AT LEAST ONE FILTER? */

      unsigned iF(0);
      for(; iF != nFilterTypes; ++iF)
        if(filterResults_[iF] && useEvents_[iF]) break;

      if(iF == nFilterTypes) continue;

      /* EVENT PASSES AT LEAST ONE FILTER - LOAD FULL EVENT AND FILL OUTPUT */

      for(iF = 0; iF != nHardPhotonFilters; ++iF){
        if(filterResults_[iF]){
          hardPhotonList_->Enter(evtTree_->GetEntries());
          break;
        }
      }
      for(iF = nHardPhotonFilters; iF != nFilterTypes; ++iF){
        if(filterResults_[iF]){
          softPhotonList_->Enter(evtTree_->GetEntries());
          break;
        }
      }

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
  hardPhotonList_->Write();
  softPhotonList_->Write();
  evtTree_->Write();
  allObjTree_->Write();
  cutTree_->Write();
  if(effTree_) effTree_->Write();

  delete outputFile_;

  outputFile_ = 0;
  hardPhotonList_ = 0;
  softPhotonList_ = 0;
  evtTree_ = 0;
  allObjTree_ = 0;
  cutTree_ = 0;
  effTree_ = 0;

  return true;
}

void
filter(TString const& _configFileName, TString const& _inputPaths, TString const& _outputName, bool _storeRadiation = false)
{
  PhotonLeptonFilter filter;
  filter.setThrow(true);
  filter.initialize(_outputName, _configFileName, _storeRadiation);
  filter.addInput(_inputPaths);
  filter.run();
  filter.clearInput();
  filter.finalize();
}
