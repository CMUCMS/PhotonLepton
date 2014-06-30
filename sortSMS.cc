#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <map>
#include <utility>

#include "SusyEvent.h"
#include "SusyTriggerEvent.h"

class SMSSorting {
public:
  SMSSorting();
  ~SMSSorting();

  bool initialize(char const*, char const*, bool);
  void addInput(char const*, char const*);
  bool run();
  void clearInput();
  bool finalize();

private:
  TString model_;
  bool singlePoint_;

  TChain input_;
  TChain triggerInput_;

  TString outputDir_;
  std::map<std::pair<int, int>, std::pair<TTree*, susy::TriggerEvent*> > output_;
};

SMSSorting::SMSSorting() :
  model_(""),
  input_("susyTree"),
  triggerInput_("triggerEventTree"),
  outputDir_(""),
  output_()
{
}

SMSSorting::~SMSSorting()
{
  for(std::map<std::pair<int, int>, std::pair<TTree*, susy::TriggerEvent*> >::iterator oItr(output_.begin()); oItr != output_.end(); ++oItr){
    delete oItr->second.first;
    delete oItr->second.second;
  }
}

bool
SMSSorting::initialize(char const* _outputDir, char const* _model, bool _singlePoint)
{
  outputDir_ = _outputDir;
  model_ = _model;
  singlePoint_ = _singlePoint;

  if(model_ != "T5wg") return false;

  return true;
}

void
SMSSorting::addInput(char const* _source, char const* _triggerSource)
{
  input_.Add(_source);
  triggerInput_.Add(_triggerSource);
}

bool
SMSSorting::run()
{
  susy::Event event;
  susy::TriggerEvent triggerEvent;

  triggerEvent.bindTree(&input_, &triggerInput_);
  event.setInput(input_);

  for(std::map<std::pair<int, int>, std::pair<TTree*, susy::TriggerEvent*> >::iterator oItr(output_.begin()); oItr != output_.end(); ++oItr)
    event.addOutput(*oItr->second.first);

  unsigned pdg1(0);
  unsigned pdg2(0);
  if(model_ == "T5wg"){
    pdg1 = 1000021;
    pdg2 = 1000024;
  }

  std::pair<std::pair<int, int>, std::pair<TTree*, susy::TriggerEvent*> > current;
  current.first = std::pair<int, int>(0, 0);
  current.second = std::pair<TTree*, susy::TriggerEvent*>(0, 0);

  long iEntry(0);
  while(event.getEntry(iEntry++) > 0){
    std::pair<int, int> masses;
    for(unsigned iG(0); iG != event.genParticles.size(); ++iG){
      susy::Particle const& part(event.genParticles[iG]);
      unsigned absId(std::abs(part.pdgId));
      if(absId == pdg1)
        masses.first = int(round(part.momentum.M() / 50.)) * 50;
      else if(absId == pdg2)
        masses.second = int(round((part.momentum.M() - 25.) / 50.)) * 50 + 25;
    }

    if(masses != current.first){
      std::map<std::pair<int, int>, std::pair<TTree*, susy::TriggerEvent*> >::iterator oItr(output_.find(masses));
      if(oItr != output_.end())
        current = *oItr;
      else{
        current.first = masses;

        TString pName;
        if(!singlePoint_){
          pName = "_" + model_;
          pName += TString::Format("_%d", masses.first);
          if(masses.second != 0) pName += TString::Format("_%d", masses.second);
        }
        pName += ".root";
        TFile::Open(outputDir_ + "/susyEvents" + pName, "recreate");
        current.second.first = new TTree("susyTree", "SUSY Event");
        current.second.first->SetAutoSave(10000000);
        current.second.second = new susy::TriggerEvent;
        current.second.second->bookTrees(outputDir_ + "/susyTriggers" + pName);
        event.addOutput(*current.second.first);

        output_.insert(current);
      }

    }

    if(!current.second.second->copyEvent(triggerEvent)) return false;
    if(current.second.first->Fill() < 0) return false;
  }

  return true;
}

void
SMSSorting::clearInput()
{
  input_.Reset();
  triggerInput_.Reset();
}

bool
SMSSorting::finalize()
{
  for(std::map<std::pair<int, int>, std::pair<TTree*, susy::TriggerEvent*> >::iterator oItr(output_.begin()); oItr != output_.end(); ++oItr){
    susy::TriggerEvent* trigEvent(oItr->second.second);
    trigEvent->write();
    delete trigEvent;

    TTree* eventTree(oItr->second.first);
    TFile* eventFile(eventTree->GetCurrentFile());
    eventFile->cd();
    eventTree->Write();
    delete eventFile;
  }

  output_.clear();

  return true;
}  

void
sortSMS(TString const& _eventSource, TString const& _outputDir, TString const& _model, bool _singlePoint)
{
  SMSSorting sort;
  if(!sort.initialize(_outputDir, _model, _singlePoint)) return;
  sort.addInput(_eventSource, TString(_eventSource).ReplaceAll("susyEvents", "susyTriggers"));
  sort.run();
  sort.clearInput();
  sort.finalize();
}
