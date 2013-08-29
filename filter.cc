#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "TObjArray.h"

#include "SusyEvent.h"

#include "../CommonCode/Utilities.h"
#include "../CommonCode/ObjectSelector.h"
#include "../CommonCode/ObjectTree.h"
#include "../CommonCode/SimpleEventProducer.h"

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <set>
#include <bitset>

template<class C>
void
sortIndicesByPt(std::vector<unsigned>& _indices, std::vector<C> const& _original)
{
  std::vector<C const*> objects;
  unsigned nI(_indices.size());
  for(unsigned iI(0); iI != nI; ++iI) objects.push_back(&_original[_indices[iI]]);

  sortByPt(objects, &_indices);
}

void
filter(TString const& _configFileName, TObjArray* _urls, TObjArray* _outputName)
{
  using namespace susy;

  enum FilterTypes {
    kPhotonAndElectron,
    kPhotonAndMuon,
    kDielectron,
    kElectronAndMuon,
    kFakePhotonAndElectron,
    kFakePhotonAndMuon,
    kPhotonAndFakeElectron,
    kPhotonAndFakeMuon,
    nFilterTypes
  };

  TString filterNames[nFilterTypes];
  filterNames[kPhotonAndElectron] = "PhotonAndElectron";
  filterNames[kPhotonAndMuon] = "PhotonAndMuon";
  filterNames[kDielectron] = "Dielectron";
  filterNames[kElectronAndMuon] = "ElectronAndMuon";
  filterNames[kFakePhotonAndElectron] = "FakePhotonAndElectron";
  filterNames[kFakePhotonAndMuon] = "FakePhotonAndMuon";
  filterNames[kPhotonAndFakeElectron] = "PhotonAndFakeElectron";
  filterNames[kPhotonAndFakeMuon] = "PhotonAndFakeMuon";

  /* SETUP FROM THE CONFIGURATION FILE */

  TString jsonPath;
  std::bitset<nFilterTypes> filterEnabled;
  std::vector<TString> hltPaths;
  std::vector<TString> gridParamNames;

  std::ifstream configFile(_configFileName);
  if(!configFile.is_open())
    throw std::invalid_argument("Cannot open config file");

  TObjArray* configRecords(new TObjArray);
  configRecords->SetOwner(true);
  TPRegexp configPat("^[ ]*([A-Z0-9_]+)[ ]*=[ ]*([^ ].*)$");
  std::string buf;
  std::stringstream bufs;
  TString line;
  while(true){
    std::getline(configFile, buf);
    if(!configFile.good()) break;
    line = buf;

    if(!configPat.MatchB(line)) continue;
    TObjArray* matches(configPat.MatchS(line));
    matches->SetOwner(true);
    TString confName(matches->At(1)->GetName());
    TString confVal(matches->At(2)->GetName());
    delete matches;

    if(confName == "JSON"){
      jsonPath = confVal;
      std::cout << line << std::endl;
      configRecords->Add(new TObjString(line));
    }
    else if(confName == "TRIGGERPATHS"){
      std::cout << "TRIGGERPATHS = ";
      bufs.clear();
      bufs.str("");
      bufs << confVal;
      while(true){
        bufs >> buf;
        std::cout << "\"" << TString(buf) << "\" ";
        hltPaths.push_back(TString(buf) + "_v*");
        if(!bufs.good()) break;
      }
      std::cout << std::endl;
      configRecords->Add(new TObjString(line));
    }
    else if(confName == "GRIDPARAMS"){
      std::cout << "GRIDPARAMS = ";
      bufs.clear();
      bufs.str("");
      bufs << confVal;
      while(true){
        bufs >> buf;
        std::cout << "\"" << TString(buf) << "\" ";
        gridParamNames.push_back(TString(buf) + "_v*");
        if(!bufs.good()) break;
      }
      std::cout << std::endl;
      configRecords->Add(new TObjString(line));
    }
    else if(confName == "FILTERS"){
      std::cout << "FILTERS = ";
      bufs.clear();
      bufs.str("");
      bufs << confVal;
      while(true){
        bufs >> buf;
        std::cout << "\"" << buf << "\" ";

	TString* pF(std::find(filterNames, filterNames + nFilterTypes, buf));
        if(pF == filterNames + nFilterTypes)
          throw std::invalid_argument(("Filter type " + buf + " not defined").c_str());

        filterEnabled[pF - filterNames] = true;

        if(!bufs.good()) break;
      }
      std::cout << std::endl;
      configRecords->Add(new TObjString(line));
    }
  }

  configFile.close();

  unsigned const nHLT(hltPaths.size());

  if(filterEnabled.none()) filterEnabled.set();

  GoodLumis goodLumis;
  if(jsonPath != ""){
    if(!goodLumis.parseJSON(jsonPath))
      throw std::invalid_argument("Failed to parse JSON");
  }

  /* INITIALIZE INPUT */

  TChain input("susyTree");
  if(input.AddFileInfoList(_urls) == 0)
    throw std::invalid_argument("Input list invalid");

  input.SetBranchStatus("*", 0);
  input.SetBranchStatus("runNumber", 1);
  input.SetBranchStatus("luminosityBlockNumber", 1);
  input.SetBranchStatus("eventNumber", 1);
  input.SetBranchStatus("isRealData", 1);
  input.SetBranchStatus("metFilter*", 1);
  input.SetBranchStatus("hlt*", 1);
  ObjectTree::setBranchStatus(input, true, true, true, false, true); // Photon, Electron, Muon, Jet, Vertex

  Event* event(new Event);
  event->setInput(input);

  // use first event to determine data / MC
  if(event->getEntry(0) <= 0)
    throw std::runtime_error("Input empty or corrupted");

  TChain fullInput("susyTree");
  fullInput.AddFileInfoList(_urls);

  fullInput.SetBranchStatus("*", 0);
  fullInput.SetBranchStatus("runNumber", 1);
  fullInput.SetBranchStatus("luminosityBlockNumber", 1);
  fullInput.SetBranchStatus("eventNumber", 1);
  fullInput.SetBranchStatus("metFilter*", 1);
  fullInput.SetBranchStatus("hlt*", 1);
  fullInput.SetBranchStatus("genParticles*", 1);
  fullInput.SetBranchStatus("pfParticles*", 1);
  fullInput.SetBranchStatus("met_pfType01CorrectedMet*", 1);
  if(!event->isRealData) fullInput.SetBranchStatus("gridParams*", 1);
  susy::ObjectTree::setBranchStatus(fullInput);

  Event* fullEvent(new Event);
  fullEvent->setInput(fullInput);

  /* SETUP OUTPUT */

  TFile* outputFile(TFile::Open(_outputName->At(0)->GetName(), "recreate"));
  if(!outputFile || outputFile->IsZombie()){
    delete outputFile;
    throw std::runtime_error("IOError");
  }

  outputFile->cd();

  for(int iR(0); iR != configRecords->GetEntries(); ++iR)
    configRecords->At(iR)->Write();
  delete configRecords;

  TTree* evtTree(new TTree("eventVars", "Event variables"));
  TTree* allObjTree(new TTree("allObjects", "All objects"));
  TTree* selectedObjTree(new TTree("selectedObjects", "Signal-quality objects"));
  TTree* fakeObjTree(new TTree("fakeObjects", "Fake proxy objects"));
  TTree* eleObjTree(new TTree("eleObjects", "Electron-fake proxy objects"));

  evtTree->SetAutoSave(10000000);
  allObjTree->SetAutoSave(10000000);
  selectedObjTree->SetAutoSave(10000000);
  fakeObjTree->SetAutoSave(10000000);
  eleObjTree->SetAutoSave(10000000);

  std::vector<unsigned> goodPhotons;
  std::vector<unsigned> fakePhotons;
  std::vector<unsigned> elePhotons;
  std::vector<unsigned> goodElectrons;
  std::vector<unsigned> fakeElectrons;
  std::vector<unsigned> goodMuons;
  std::vector<unsigned> fakeMuons;
  std::vector<unsigned> goodJets;
  std::vector<unsigned> goodVertices;

  SimpleEventProducer output;

  TString hltPathList[] = {
    "HLT_IsoMu24_eta2p1",
    "HLT_IsoMu24",
    "HLT_Ele27_WP80",
    "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50",
    "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50",
    "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60",
    "HLT_Photon26_R9Id85_IR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60",
    "HLT_Mu22_Photon22_CaloIdL",
    "HLT_Photon70_CaloIdXL_PFHT400",
    "HLT_Photon70_CaloIdXL_PFHT500",
    "HLT_Photon70_CaloIdXL_PFNoPUHT400",
    "HLT_Photon70_CaloIdXL_PFNoPUHT500",
    "HLT_Photon70_CaloIdXL_PFMET100",
    "HLT_Photon60_CaloIdL_HT300",
    "HLT_Photon60_CaloIdL_MHT70",
    "HLT_Photon135",
    "HLT_Photon150"
  };
  output.setHLTPaths(std::vector<TString>(hltPathList, hltPathList + sizeof(hltPathList) / sizeof(TString)));
  if(!event->isRealData) output.setGridParams(gridParamNames);

  output.setPhotonId(susy::PhLoose12LV);

  output.setSavePF(true);

  output.initialize(evtTree, 0, allObjTree, event->isRealData);
  output.addPreselected(*fakeObjTree, event->isRealData, &fakePhotons, &fakeElectrons, &fakeMuons, 0, 0);
  output.addPreselected(*selectedObjTree, event->isRealData, &goodPhotons, &goodElectrons, &goodMuons, &goodJets, &goodVertices);
  output.addPreselected(*eleObjTree, event->isRealData, &elePhotons, 0, 0, 0, 0);

  bool filterResults[nFilterTypes];
  for(unsigned iF(0); iF != nFilterTypes; ++iF)
    evtTree->Branch(filterNames[iF], filterResults + iF, filterNames[iF] + "/O");

  /* START LOOP */

  enum CountPoints {
    kAllEvents, //0
    kGoodLumi, //1
    kMetFilter, //2
    kHLT, //3
    kGoodVertex, //4
    kGoodElectron, //5
    kFakeElectron, //6
    kGoodMuon, //7
    kFakeMuon, //8
    kGoodPhoton, //9
    kFakePhoton, //10
    kElePhoton, //11
    kFinal, //12-(11+nFilterTypes)
    nCountPoints = kFinal + nFilterTypes
  };

  unsigned eventCounter[nCountPoints];
  for(unsigned i(0); i != nCountPoints; ++i) eventCounter[i] = 0;

  long iEvent(0);
  int iTree(-1);
  while(event->getEntry(iEvent++) > 0){
    if(iEvent % 1000 == 0) std::cout << "Analyzing event " << iEvent << std::endl;

    if(iTree != input.GetTreeNumber()){
      std::cout << "Current file: " << input.GetCurrentFile()->GetName() << std::endl;
      iTree = input.GetTreeNumber();
    }

    ++eventCounter[kAllEvents];

    if(event->isRealData && !goodLumis.isGoodLumi(event->runNumber, event->luminosityBlockNumber)) continue;

    ++eventCounter[kGoodLumi];

    if(!event->passMetFilters()) continue;

    ++eventCounter[kMetFilter];

    if(nHLT != 0){
      unsigned iHLT(0);
      for(; iHLT != nHLT; ++iHLT)
        if(event->hltMap.pass(hltPaths[iHLT])) break;
      if(iHLT == nHLT) continue;
    }

    ++eventCounter[kHLT];

    goodPhotons.clear();
    fakePhotons.clear();
    elePhotons.clear();
    goodElectrons.clear();
    fakeElectrons.clear();
    goodMuons.clear();
    fakeMuons.clear();
    goodJets.clear();
    goodVertices.clear();

    if(countVertices(event->vertices, &goodVertices) == 0) continue;

    ++eventCounter[kGoodVertex];

    /* SELECT ELECTRONS */

    ElectronCollection const& electrons(event->electrons["gsfElectrons"]);

    // for electron duplicate problem in 52X fastsim
    std::set<short> superClusterIndices;
    std::set<short> trackIndices;

    std::vector<Electron const*> electronsNoIso;

    std::bitset<nElectronCriteria> elIdResults;
    std::bitset<nElectronCriteria> elBaseline;
    elBaseline.set(ElFiducial);
    elBaseline.set(ElCombIso);
    elBaseline.set(ElSigmaIetaIeta);
    elBaseline.set(ElHOverE);
    std::bitset<nElectronCriteria> elNoIso(ObjectSelector::elReferences[ElMedium12]);
    elNoIso.reset(ElCombIso);

    float minDeltaEta[2] = {0.007, 0.01};
    float minDeltaPhi[2] = {0.8, 0.7};

    unsigned nEl(electrons.size());
    for(unsigned iEl(0); iEl < nEl; ++iEl){
      Electron const& el(electrons[iEl]);

      // FIX FOR 52X FASTSIM BUG
      if(!event->isRealData){
        if(superClusterIndices.find(el.superClusterIndex) != superClusterIndices.end() ||
           trackIndices.find(el.gsfTrackIndex) != trackIndices.end()) continue;
        superClusterIndices.insert(el.superClusterIndex);
        trackIndices.insert(el.gsfTrackIndex);
      }

      // FIX FOR NAN-MOMENTUM PROBLEM
      if(el.momentum.X() != el.momentum.X()) continue;

      ElectronVars vars(el, *event);

      if(vars.iSubdet == -1) continue;

      bool isGood(ObjectSelector::isGoodElectron(vars, ElMedium12, &elIdResults));

      if((elIdResults & elNoIso) == elNoIso) electronsNoIso.push_back(&el);

      if(vars.pt < 25.) continue;

      if(isGood)
        goodElectrons.push_back(iEl);
      else if((elIdResults & elBaseline) == elBaseline &&
              vars.deltaEta > minDeltaEta[vars.iSubdet] && vars.deltaPhi > minDeltaPhi[vars.iSubdet] &&
              vars.d0 < 0.2 && vars.dz < 1.)
        fakeElectrons.push_back(iEl);
    }

    sortIndicesByPt(goodElectrons, electrons);
    sortIndicesByPt(fakeElectrons, electrons);

    unsigned nGoodElectrons(goodElectrons.size());
    unsigned nFakeElectrons(fakeElectrons.size());
    unsigned nNoIsoElectrons(electronsNoIso.size());

    bool hasGoodElectron(nGoodElectrons != 0);
    bool hasFakeElectron(nFakeElectrons != 0);

    if(hasGoodElectron) ++eventCounter[kGoodElectron];
    if(hasFakeElectron) ++eventCounter[kFakeElectron];

    /* SELECT MUONS */

    MuonCollection const& muons(event->muons["muons"]);

    std::vector<Muon const*> muonsNoIso;

    std::bitset<nMuonCriteria> muIdResults;
    std::bitset<nMuonCriteria> muBaseline(ObjectSelector::muReferences[MuTight12]);
    muBaseline.reset(MuCombIso);

    unsigned nMu(muons.size());
    for(unsigned iMu(0); iMu < nMu; ++iMu){
      Muon const& mu(muons[iMu]);

      MuonVars vars(mu, *event);

     if(vars.iSubdet == -1) continue;

      bool isGood(ObjectSelector::isGoodMuon(vars, MuTight12, &muIdResults));

      if((muIdResults & muBaseline) == muBaseline)
        muonsNoIso.push_back(&mu);

      if(vars.pt < 25.) continue;

      if(isGood)
        goodMuons.push_back(iMu);
      else if((muIdResults & muBaseline) == muBaseline)
        fakeMuons.push_back(iMu);
    }

    sortIndicesByPt(goodMuons, muons);
    sortIndicesByPt(fakeMuons, muons);

    unsigned nGoodMuons(goodMuons.size());
    unsigned nNoIsoMuons(muonsNoIso.size());

    bool hasGoodMuon(nGoodMuons != 0);
    bool hasFakeMuon(fakeMuons.size() != 0);

    if(hasGoodMuon) ++eventCounter[kGoodMuon];
    if(hasFakeMuon) ++eventCounter[kFakeMuon];

    /* SELECT PHOTONS */

    PhotonCollection const& photons(event->photons["photons"]);

    std::bitset<nPhotonCriteria> phIdResults;
    std::bitset<nPhotonCriteria> phIsolation;
    phIsolation.set(PhChargedHadronIso);
    phIsolation.set(PhNeutralHadronIso);
    phIsolation.set(PhPhotonIso);
    std::bitset<nPhotonCriteria> phNoVeto(ObjectSelector::phReferences[PhLoose12LV]);
    phNoVeto.reset(PhElectronVeto);

    unsigned nPh(photons.size());
    for(unsigned iPh(0); iPh < nPh; ++iPh){
      Photon const& ph(photons[iPh]);

      if(ph.momentum.Pt() < 25.) continue;
      if(std::abs(ph.caloPosition.Eta()) > etaGapBegin) continue;

      /* veto bremming electrons */
      /* must allow the ecalDriven electron that shares the superCluster - this will be vetoed in the electronVeto if missingHits <= 1 */
      /* superCluster matching will miss the trackerDriven electrons. */
      /* so if a trackerDriven electron of medium sans isolation quality is pointing at this photon, we will drop this, which should happen because such case will
         also pass the c-safe electron veto and lead to double counting. */
      /* after this cut, any selected photon/fake candidates will have no >= medium\iso electrons within dR < 0.3 except for the one that shares the supercluster with the photon/fake */
      unsigned iL(0);
      for(; iL != nNoIsoElectrons; ++iL){
        if(electronsNoIso[iL]->superClusterIndex == ph.superClusterIndex) continue;
        if(electronsNoIso[iL]->momentum.DeltaR(ph.momentum) < 0.3) break;
      }
      if(iL != nNoIsoElectrons) continue;

      iL = 0;
      for(; iL != nNoIsoMuons; ++iL)
        if(muonsNoIso[iL]->momentum.DeltaR(ph.momentum) < 0.3) break;
      if(iL != nNoIsoMuons) continue;

      PhotonVars vars(ph, *event);

      bool isGood(ObjectSelector::isGoodPhoton(vars, PhLoose12LV, &phIdResults));

      if(isGood)
        goodPhotons.push_back(iPh); // c-safe electron ID will not allow medium\iso electron to share the supercluster with this photon
      else if((phIdResults & phNoVeto) == phNoVeto)
        elePhotons.push_back(iPh);
      else if(phIdResults[PhFiducial] && vars.hOverE > 0. && (phIdResults & phIsolation).count() == 1)
        fakePhotons.push_back(iPh);
    }

    sortIndicesByPt(goodPhotons, photons);
    sortIndicesByPt(fakePhotons, photons);
    sortIndicesByPt(elePhotons, photons);

    unsigned nGoodPhotons(goodPhotons.size());
    unsigned nFakePhotons(fakePhotons.size());
    unsigned nElePhotons(elePhotons.size());

    bool hasGoodPhoton(nGoodPhotons != 0);
    bool hasFakePhoton(nFakePhotons != 0);
    bool hasElePhoton(nElePhotons != 0);

    if(hasGoodPhoton) ++eventCounter[kGoodPhoton];
    if(hasFakePhoton) ++eventCounter[kFakePhoton];
    if(hasElePhoton) ++eventCounter[kElePhoton];

    /* DETERMINE RESULT OF EACH FILTER */

    filterResults[kPhotonAndElectron] = hasGoodPhoton && hasGoodElectron && photons[goodPhotons[0]].momentum.Pt() > 40.;
    filterResults[kPhotonAndMuon] = hasGoodPhoton && hasGoodMuon;
    filterResults[kDielectron] = false;
    if(hasElePhoton && hasGoodElectron && photons[elePhotons[0]].momentum.Pt() > 40.){
      if(nElePhotons > 1 || nGoodElectrons > 1)
        filterResults[kDielectron] = true;
      else if(electrons[goodElectrons[0]].superClusterIndex != photons[elePhotons[0]].superClusterIndex)
        filterResults[kDielectron] = true;
    }
    filterResults[kElectronAndMuon] = hasGoodElectron && hasGoodMuon; 
    filterResults[kFakePhotonAndElectron] = false;
    if(hasFakePhoton && hasGoodElectron && photons[fakePhotons[0]].momentum.Pt() > 40.){
      if(nFakePhotons > 1 || nGoodElectrons > 1)
        filterResults[kFakePhotonAndElectron] = true;
      else if(electrons[goodElectrons[0]].superClusterIndex != photons[fakePhotons[0]].superClusterIndex)
        filterResults[kFakePhotonAndElectron] = true;
    }
    filterResults[kFakePhotonAndMuon] = hasFakePhoton && hasGoodMuon;
    filterResults[kPhotonAndFakeElectron] = false;
    if(hasGoodPhoton && hasFakeElectron && photons[goodPhotons[0]].momentum.Pt() > 40.){
      if(nGoodPhotons > 1 || nFakeElectrons > 1)
        filterResults[kPhotonAndFakeElectron] = true;
      else if(electrons[fakeElectrons[0]].superClusterIndex != photons[goodPhotons[0]].superClusterIndex)
        filterResults[kPhotonAndFakeElectron] = true;
    }
    filterResults[kPhotonAndFakeMuon] = hasGoodPhoton && hasFakeMuon;

    std::bitset<nFilterTypes> filterBits;
    for(unsigned iF(0); iF != nFilterTypes; ++iF){
      if(filterResults[iF]){
        ++eventCounter[kFinal + iF];
        filterBits.set(iF);
      }
    }

    if((std::bitset<nFilterTypes>(filterBits) & filterEnabled).none()) continue;

    /* EVENT PASSES AT LEAST ONE FILTER - LOAD FULL EVENT AND FILL OUTPUT */

    fullEvent->getEntry(iEvent - 1);

    PFJetCollection const& jets(fullEvent->pfJets["ak5"]);
    unsigned nJ(jets.size());
    for(unsigned iJ(0); iJ != nJ; ++iJ){
      PFJet const& jet(jets[iJ]);

      JetVars vars(jet, *event);
      if(!vars.isLoose) continue;
      if(vars.pt < 25. || std::abs(vars.eta) > 3.) continue;

      unsigned iG;

      iG = 0;
      for(; iG != nGoodElectrons; ++iG){
        unsigned iEl(goodElectrons[iG]);
        if(electrons[iEl].momentum.DeltaR(jet.momentum) < 0.3) break;
      }
      if(iG != nGoodElectrons) continue;

      iG = 0;
      for(; iG != nGoodMuons; ++iG){
        unsigned iMu(goodMuons[iG]);
        if(muons[iMu].momentum.DeltaR(jet.momentum) < 0.3) break;
      }
      if(iG != nGoodMuons) continue;

      iG = 0;
      for(; iG != nGoodPhotons; ++iG){
        unsigned iPh(goodPhotons[iG]);
        if(photons[iPh].momentum.DeltaR(jet.momentum) < 0.3) break;
      }
      if(iG != nGoodPhotons) continue;

      iG = 0;
      for(; iG != nElePhotons; ++iG){
        unsigned iPh(elePhotons[iG]);
        if(photons[iPh].momentum.DeltaR(jet.momentum) < 0.3) break;
      }
      if(iG != nElePhotons) continue;

      goodJets.push_back(iJ);
    }

    output.produce(*fullEvent);

    evtTree->Fill();
    allObjTree->Fill();
    selectedObjTree->Fill();
    fakeObjTree->Fill();
    eleObjTree->Fill();
  }

  std::cout << "Cut flow: ";
  for(unsigned i(0); i < nCountPoints; ++i)
    std::cout << "[" << i << "]: " << eventCounter[i] << " ";
  std::cout << std::endl;

  /* FINALIZE OUTPUT */

  outputFile->cd();
  outputFile->Write();
  delete outputFile;

  delete event;
  delete fullEvent;
}

void
filter(TString const& _configFileName, TString const& _dataset, TObjArray* _fileNames, TObjArray* _outputDir)
{
  // For "externalList" mode of dcmu job scheduler
  // The list must have the grid point name in the (N+1)n-th rows (N=files per point, n=0,1,...) and the file names in the folloing N rows
  // example: TChiwg.list
  //   TChiwg_1000_1000
  //   TChiwg_1000_1000_10_4_2zc.root
  //   ...
  // Then the submit command is
  // submitDCMUJobs.py -n 51 -J TChiwg -d 20 -x "$PWD/TChiwg.list" -a '"'$PWD'/filter.cfg"' -s 'output_directory' -z /store/RA3Ntuples/... filter.cc

  TObjArray output;
  output.SetOwner();

  TString gridPoint(_fileNames->At(0)->GetName());

  TObjString* outstr(new TObjString(_outputDir->At(0)->GetName()));
  outstr->String() += "/skim_" + gridPoint + ".root";

  output.Add(outstr);

  TObjArray input;
  input.SetOwner();

  for(int iF(1); iF != _fileNames->GetEntries(); ++iF)
    input.Add(new TObjString(_dataset + "/" + _fileNames->At(iF)->GetName()));

  std::cout << "Running filter on " << gridPoint << std::endl;

  filter(_configFileName, &input, &output);
}
