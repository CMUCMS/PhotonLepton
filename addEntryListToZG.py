import os
import re
import ROOT

sourceDir = '/store/glskim/'

for fileName in os.listdir(sourceDir):
    if not re.match('ZGToLLG(?:|_[0-9]+)[.]root$', fileName): continue

    newName = fileName.replace('ZGToLLG', 'ZGToLLG_PtG-5-130')
    
    skimFile = ROOT.TFile.Open(sourceDir + fileName, 'update')
    if skimFile.Get('entrylist'):
        skimFile.Close()
        continue

    tree = skimFile.Get('eventVars')
    tree.Draw('>>entrylist', 'lhePhotonPt < 130.', 'entrylist')
    
    elist = ROOT.gDirectory.Get('entrylist')
    elist.SetFileName('rooth://ncmu40/' + sourceDir + newName)

    elist.Write()
    skimFile.Close()

    os.rename(sourceDir + fileName, sourceDir + newName)
