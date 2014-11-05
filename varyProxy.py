import ROOT

for dataset, lepton in [('DataE', 'Electron'), ('DataM', 'Muon')]:
    for sample in ['FakePhotonAnd', 'PhotonAndFake']:

        nomSource = ROOT.TFile.Open('rooth://ncmu40//store/glweighted/' + dataset + '_' + sample + lepton + '.root')
        altSource = ROOT.TFile.Open('rooth://ncmu40//store/glweighted_altfake/' + dataset + '_' + sample + lepton + '.root')

        nom = nomSource.Get('eventList')
        alt = altSource.Get('eventList')

        if lepton == 'Electron':
            baseline = '(mass2 < 81. || mass2 > 101.)'
        else:
            baseline = '1'

        cut = 'met > 120. && mt > 100.'

        nomN = nom.GetEntries(baseline + ' && ' + cut)
        nomD = nom.GetEntries(baseline)
        altN = alt.GetEntries(baseline + ' && ' + cut)
        altD = alt.GetEntries(baseline)

        print sample + lepton, nomN, '/', nomD, '=', float(nomN) / nomD, '->', altN, '/', altD, '=', float(altN) / altD

        nomSource.Close()
        altSource.Close()
