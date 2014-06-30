import sys
import os
import re

dirs = sys.argv[1:]

points = [
    'AllEvents',
    'GoodLumi',
    'MetFilter',
    'HLT',
    'GoodVertex',
    'RadiationVeto',
    'GoodElectron',
    'FakeElectron',
    'GoodMuon',
    'FakeMuon',
    'FlavorConflict',
    'ChargedHadronVeto',
    'GoodPhoton',
    'FakePhoton',
    'ElePhoton',
    'PassOneFilter',
    'PassPhotonAndElectron',
    'PassPhotonAndMuon',
    'PassElePhotonAndElectron',
    'PassElePhotonAndMuon',
    'PassFakePhotonAndElectron',
    'PassFakePhotonAndMuon',
    'PassPhotonAndFakeElectron',
    'PassPhotonAndFakeMuon'
]

counts = dict([(key, 0) for key in points])

countRe = re.compile('\[?([^\]]+)\]?: ([0-9]+)')

for dir in dirs:
    dir = os.path.realpath(dir)
    for fileName in os.listdir(dir):
        localCounts = dict([(key, 0) for key in points])
        with open(dir + '/' + fileName) as logFile:
            for line in logFile:
                matches = countRe.match(line)
                if matches is None: continue
    
                localCounts[matches.group(1)] = int(matches.group(2))
    
        for key in points:
            counts[key] += localCounts[key]

for key in points:
    print key + ':', counts[key]
