// Sort full model by sparticle production process

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <iostream>
#include <map>
#include <utility>

#include "SusyEvent.h"
#include "SusyTriggerEvent.h"

#include "Toolset/GenTreeViewer/test/GenDecayFilterRA3.cc"

enum Process {
  NeuChip,
  NeuChim,
  ChipChim,
  NeuGlu,
  ChipGlu,
  ChimGlu,
  GluGlu,
  SqGlu,
  nProcesses
};

TString processNames[nProcesses] = {
  "ncp",
  "ncm",
  "cc",
  "ng",
  "cpg",
  "cmg",
  "gg",
  "sg"
};

class FullModelSorting {
public:
  FullModelSorting();
  ~FullModelSorting();

  bool initialize(char const*);
  void addInput(char const*, char const*);
  bool run();
  void clearInput();
  bool finalize();

private:
  TChain input_;
  TChain triggerInput_;

  TString outputDir_;
  std::map<int, std::pair<TTree*, susy::TriggerEvent*> > output_;

  GenDecaySorterRA3 sorter_;
};

FullModelSorting::FullModelSorting() :
  input_("susyTree"),
  triggerInput_("triggerEventTree"),
  outputDir_(""),
  output_(),
  sorter_()
{
}

FullModelSorting::~FullModelSorting()
{
  for(std::map<int, std::pair<TTree*, susy::TriggerEvent*> >::iterator oItr(output_.begin()); oItr != output_.end(); ++oItr){
    delete oItr->second.first;
    delete oItr->second.second;
  }
}

bool
FullModelSorting::initialize(char const* _outputDir)
{
  outputDir_ = _outputDir;

  sorter_.addCategory(NeuChip, "j>1000022 && j>+1000024");
  sorter_.addCategory(NeuChim, "j>1000022 && j>-1000024");
  sorter_.addCategory(ChipChim, "j>+1000024 && j>-1000024");
  sorter_.addCategory(NeuGlu, "j>1000022 && j>1000021");
  sorter_.addCategory(ChipGlu, "j>+1000024 && j>1000021");
  sorter_.addCategory(ChimGlu, "j>-1000024 && j>1000021");
  sorter_.addCategory(GluGlu, "2 * (j>1000021)");
  sorter_.addCategory(SqGlu, "(j>1000001 || j>1000002 || j>1000003 || j>1000004 || j>1000005 || j>1000006 || j>2000001 || j>2000002 || j>2000003 || j>2000004 || j>2000005 || j>2000006) && j>1000021");

  return true;
}

void
FullModelSorting::addInput(char const* _source, char const* _triggerSource)
{
  input_.Add(_source);
  triggerInput_.Add(_triggerSource);
}

bool
FullModelSorting::run()
{
  susy::Event event;
  susy::TriggerEvent triggerEvent;

  triggerEvent.bindTree(&input_, &triggerInput_);
  event.setInput(input_);

  for(std::map<int, std::pair<TTree*, susy::TriggerEvent*> >::iterator oItr(output_.begin()); oItr != output_.end(); ++oItr)
    event.addOutput(*oItr->second.first);

  std::pair<int, std::pair<TTree*, susy::TriggerEvent*> > current;
  current.first = -1;
  current.second = std::pair<TTree*, susy::TriggerEvent*>(0, 0);

  long iEntry(0);
  while(event.getEntry(iEntry++) > 0){
    int proc(sorter_.sort(event));
    if(proc < 0) continue;

    if(proc != current.first){
      std::map<int, std::pair<TTree*, susy::TriggerEvent*> >::iterator oItr(output_.find(proc));
      if(oItr != output_.end())
        current = *oItr;
      else{
        current.first = proc;

        TString pName(gSystem->BaseName(input_.GetCurrentFile()->GetName()));
        pName.ReplaceAll(".root", "").ReplaceAll("susyEvents", "");
        pName += "_" + processNames[proc] + ".root";

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
FullModelSorting::clearInput()
{
  input_.Reset();
  triggerInput_.Reset();
}

bool
FullModelSorting::finalize()
{
  for(std::map<int, std::pair<TTree*, susy::TriggerEvent*> >::iterator oItr(output_.begin()); oItr != output_.end(); ++oItr){
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
sortFullModel(TString const& _eventSource, TString const& _outputDir)
{
  FullModelSorting sort;
  if(!sort.initialize(_outputDir)) return;
  sort.addInput(_eventSource, TString(_eventSource).ReplaceAll("susyEvents", "susyTriggers"));
  sort.run();
  sort.clearInput();
  sort.finalize();
}
