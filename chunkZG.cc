#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"

void
chunkZG(TString const& _sourceName, TString const& _outputName)
{
  TFile* source(TFile::Open(_sourceName));
  TTree* eventIn(static_cast<TTree*>(source->Get("eventVars")));
  TTree* objIn(static_cast<TTree*>(source->Get("allObjects")));
  TTree* cutIn(static_cast<TTree*>(source->Get("cutTree")));

  TFile* outputFile(TFile::Open(_outputName, "recreate"));
  TTree* eventOut(eventIn->CloneTree(0));
  TTree* objOut(objIn->CloneTree(0));
  TTree* cutOut(cutIn->CloneTree(0));

  TChain skimInput("eventVars");
  skimInput.Add(_sourceName);

  skimInput.SetBranchStatus("*", 0);
  skimInput.SetBranchStatus("lhePhotonPt", 1);

  float lhePhotonPt;

  skimInput.SetBranchAddress("lhePhotonPt", &lhePhotonPt);

  long iEntry(0);
  while(skimInput.GetEntry(iEntry++) > 0){
    if(lhePhotonPt > 130.) continue;

    eventIn->GetEntry(iEntry - 1);
    objIn->GetEntry(iEntry - 1);
    cutIn->GetEntry(iEntry - 1);

    eventOut->Fill();
    objOut->Fill();
    cutOut->Fill();
  }

  delete source;

  outputFile->cd();
  outputFile->Write();
  delete outputFile;
}
