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
    hardPhotons = skimFile.Get('hardPhotons')
    softPhotons = skimFile.Get('softPhotons')

    cut = ''

    if dataset == 'ZGToLLG':
        oldName = fileName
        fileName = fileName.replace('ZGToLLG', 'ZGToLLG_PtG-5-130')

        cut = 'lhePhotonPt < 130.'

    elif dataset == 'TTJetsFullLept' or dataset == 'TTJetsSemiLept' or dataset == 'WW':
        cut = 'radPt < 20. || isFSR'

    if cut:
        tree = skimFile.Get('eventVars')

        hardPhotons.SetFileName(sourceDir + fileName)
        hardPhotons.SetName('oldHardPhotons')
        tree.SetEntryList(hardPhotons)
        tree.GetEntryNumber(0) # needed to initialize the list
        tree.Draw('>>hardPhotons', cut, 'entrylist')
        hardPhotons = ROOT.gDirectory.Get('hardPhotons')
    
        softPhotons.SetFileName(sourceDir + fileName)
        softPhotons.SetName('oldSoftPhotons')
        tree.SetEntryList(softPhotons)
        tree.GetEntryNumber(0) # needed to initialize the list
        tree.Draw('>>softPhotons', cut, 'entrylist')        
        softPhotons = ROOT.gDirectory.Get('softPhotons')

        tree.SetEntryList(0)

    hardPhotons.SetFileName('rooth://ncmu40/' + sourceDir + fileName)
    hardPhotons.Write()
    softPhotons.SetFileName('rooth://ncmu40/' + sourceDir + fileName)
    softPhotons.Write()

    skimFile.Close()
    
    if dataset == 'ZGToLLG':
        os.rename(sourceDir + oldName, sourceDir + fileName)
