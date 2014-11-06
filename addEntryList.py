import sys
import os
import re
import ROOT

# TTJets SingleE SingleM WW

sourceDir = '/store/glskim/'

dataset = sys.argv[1]

for fileName in os.listdir(sourceDir):
    if not re.match(dataset + '(?:|_[0-9]+)[.]root$', fileName): continue

    print fileName

    skimFile = ROOT.TFile.Open(sourceDir + fileName, 'update')
    tree = skimFile.Get('eventVars')

    if skimFile.Get('entrylist'):
        tree.SetEntryList(skimFile.Get('entrylist'))

    tree.Draw('>>hardPhotonList', 'PhotonAndElectron || PhotonAndMuon || ElePhotonAndElectron || ElePhotonAndMuon || FakePhotonAndElectron || FakePhotonAndMuon || PhotonAndFakeElectron || PhotonAndFakeMuon || ElePhotonAndFakeElectron || ElePhotonAndFakeMuon || FakePhotonAndFakeElectron || FakePhotonAndFakeMuon', 'entrylist')
    elist = ROOT.gDirectory.Get('hardPhotonList')
    elist.SetTitle('Hard photon events')
    elist.SetFileName('rooth://ncmu40/' + sourceDir + fileName)
    elist.Write()

    tree.Draw('>>softPhotonList', 'SoftPhotonAndElectron || SoftPhotonAndMuon || SoftElePhotonAndElectron || SoftElePhotonAndMuon || SoftFakePhotonAndElectron || SoftFakePhotonAndMuon || SoftPhotonAndFakeElectron || SoftPhotonAndFakeMuon', 'entrylist')
    elist = ROOT.gDirectory.Get('softPhotonList')
    elist.SetTitle('Soft photon events')
    elist.SetFileName('rooth://ncmu40/' + sourceDir + fileName)
    elist.Write()

    skimFile.Close()
