import sys
import os
import re
import ROOT

sourceDir = '/store/glskim/'

dataset = sys.argv[1]

for fileName in os.listdir(sourceDir):
    if not re.match(dataset + '(?:|_[0-9]+)[.]root$', fileName): continue

    print fileName

    skimFile = ROOT.TFile.Open(sourceDir + fileName, 'update')

    cut = ''

    if dataset == 'ZGToLLG':
        oldName = fileName
        fileName = fileName.replace('ZGToLLG', 'ZGToLLG_PtG-5-130')

        cut = 'lhePhotonPt < 130.'

    elif dataset == 'TTJetsFullLept' or dataset == 'TTJetsSemiLept' or dataset == 'WW':
        cut = 'radPt < 20. || isFSR'

    if cut:
        tree = skimFile.Get('eventVars')

        tree.Draw('>>hardPhotonList', '(PhotonAndElectron || PhotonAndMuon || ElePhotonAndElectron || ElePhotonAndMuon || FakePhotonAndElectron || FakePhotonAndMuon || PhotonAndFakeElectron || PhotonAndFakeMuon || ElePhotonAndFakeElectron || ElePhotonAndFakeMuon || FakePhotonAndFakeElectron || FakePhotonAndFakeMuon) && (' + cut + ')', 'entrylist')
        hardPhotonList = ROOT.gDirectory.Get('hardPhotonList')
        hardPhotonList.SetTitle('Hard photon events')

        tree.Draw('>>softPhotonList', '(SoftPhotonAndElectron || SoftPhotonAndMuon || SoftElePhotonAndElectron || SoftElePhotonAndMuon || SoftFakePhotonAndElectron || SoftFakePhotonAndMuon || SoftPhotonAndFakeElectron || SoftPhotonAndFakeMuon) && (' + cut + ')', 'entrylist')
        softPhotonList = ROOT.gDirectory.Get('softPhotonList')
        softPhotonList.SetTitle('Soft photon events')

    else:
        hardPhotonList = skimFile.Get('hardPhotonList')
        softPhotonList = skimFile.Get('softPhotonList')

    hardPhotonList.SetFileName('rooth://ncmu40/' + sourceDir + fileName)
    hardPhotonList.Write()
    softPhotonList.SetFileName('rooth://ncmu40/' + sourceDir + fileName)
    softPhotonList.Write()

    skimFile.Close()
    
    if dataset == 'ZGToLLG':
        os.rename(sourceDir + oldName, sourceDir + fileName)
