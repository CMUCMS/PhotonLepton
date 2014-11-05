import re
import array
import sys

import ROOT

ROOT.gROOT.SetBatch(True)

remoteDir = 'rooth://ncmu40//store/glskim/'
outputDir = '/afs/cern.ch/user/y/yiiyama/output/GammaL/main/idEfficiencies/'

if len(sys.argv) == 2:
    inputName = sys.argv[1]
    outputName = sys.argv[1]
else:
    remoteDir += sys.argv[1] + '/'
    outputDir += sys.argv[1] + '/'
    inputName = 'skim_' + sys.argv[2]
    outputName = sys.argv[2]

photonPtBins = array.array('d', [40., 50., 8000.])
photonEtaBins = array.array('d', [0., 0.8, 1.4442])
electronPtBins = array.array('d', [20., 30., 40., 50., 8000.])
electronEtaBins = array.array('d', [0., 0.8, 1.442, 1.556, 2., 2.5])
muonPtBins = array.array('d', [25., 30., 35., 40., 50., 60., 90., 140., 8000.])
muonEtaBins = array.array('d', [0., 0.9, 1.2, 2.1, 2.4])

tree = ROOT.TChain('effTree')

dp = ROOT.gSystem.OpenDirectory(remoteDir)
if not dp:
    raise RuntimeError('Remote directory not found')

while True:
    entry = ROOT.gSystem.GetDirEntry(dp)
    if not entry: break
    if entry == '.' or entry == '..': continue
    if re.match(inputName + '(?:|_[0-9]+)[.]root$', entry):
        tree.Add(remoteDir + entry)

ROOT.gSystem.FreeDirectory(dp)

outputFile = ROOT.TFile.Open(outputDir + outputName + '.root', 'recreate')

photon_eff = ROOT.TProfile2D('photon_eff', 'Photon ID efficiency', len(photonPtBins) - 1, photonPtBins, len(photonEtaBins) - 1, photonEtaBins)
electron_eff = ROOT.TProfile2D('electron_eff', 'Electron ID efficiency', len(electronPtBins) - 1, electronPtBins, len(electronEtaBins) - 1, electronEtaBins)
muon_eff = ROOT.TProfile2D('muon_eff', 'Muon ID efficiency', len(muonPtBins) - 1, muonPtBins, len(muonEtaBins) - 1, muonEtaBins)

photon_eff.Sumw2()
electron_eff.Sumw2()
muon_eff.Sumw2()

tree.Draw('pass:TMath::Abs(eta):pt>>+photon_eff', 'puWeight * (pdgId == 22 && genIso < 5. && reco && pt > 40. && TMath::Abs(eta) < 1.4442)', 'prof goff')
tree.Draw('pass:TMath::Abs(eta):pt>>+electron_eff', 'puWeight * (TMath::Abs(pdgId) == 11 && reco && pt > 25. && (TMath::Abs(eta) < 1.4442 || (TMath::Abs(eta) > 1.56 && TMath::Abs(eta) < 2.5)))', 'prof goff')
tree.Draw('pass:TMath::Abs(eta):pt>>+muon_eff', 'puWeight * (TMath::Abs(pdgId) == 13 && reco && pt > 25. && TMath::Abs(eta) < 2.4)', 'prof goff')

outputFile.cd()
outputFile.Write()
outputFile.Close()
