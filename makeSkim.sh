#!/bin/bash

spec=_cflow

JOBS="$@"

get_arguments() {
    vetoThreshold="0."
    format='susyEvents*.root,susyTriggers*.root'
    reducer=Hadder
    target=''
    case $1 in
        PhotonA)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012A-22Jan2013-v1/Photon/Photon36Photon22_FT_53_V21_AN4
            cfg="Data.cfg Electron.cfg"
            ;;
        DoublePhotonB)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/Photon36Photon22_FT_53_V21_AN4
            cfg="Data.cfg Electron.cfg"
            ;;
        DoublePhotonC)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Run2012C-22Jan2013-v2/DoublePhoton/Photon36Photon22_FT_53_V21_AN4
            cfg="Data.cfg Electron.cfg"
            ;;
        DoublePhotonD)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Run2012D-22Jan2013-v1/DoublePhoton/Photon36Photon22_FT_53_V21_AN4
            cfg="Data.cfg Electron.cfg"
            ;;
        MuEGA)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Run2012A-22Jan2013-v1/MuEG/Mu22Photon22_FT_53_V21_AN4
            cfg="Data.cfg Muon.cfg"
            ;;
        MuEGB)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Run2012B-22Jan2013-v1/MuEG/Mu22Photon22_FT_53_V21_AN4
            cfg="Data.cfg Muon.cfg"
            ;;
        MuEGC)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Run2012C-22Jan2013-v1/MuEG/Mu22Photon22_FT_53_V21_AN4
            cfg="Data.cfg Muon.cfg"
            ;;
        MuEGD)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Run2012D-22Jan2013-v1/MuEG/Mu22Photon22_FT_53_V21_AN4
            cfg="Data.cfg Muon.cfg"
            ;;
        SingleMuA)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012A-22Jan2013-v1/SingleMu/IsoMu24_FT_53_V21_AN4
            cfg="Data.cfg FakePhoton.cfg"
            ;;
        SingleMuB)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/SingleMu/IsoMu24_FT_53_V21_AN4
            cfg="Data.cfg FakePhoton.cfg"
            ;;
        SingleMuC)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012C-22Jan2013-v1/SingleMu/IsoMu24_FT_53_V21_AN4
            cfg="Data.cfg FakePhoton.cfg"
            ;;
        SingleMuD)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012D-22Jan2013-v1/SingleMu/IsoMu24_FT_53_V21_AN4
            cfg="Data.cfg FakePhoton.cfg"
            ;;
        SingleElectronA)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012A-22Jan2013-v1/SingleElectron/Ele27WP80_FT_53_V21_AN4
            cfg="Data.cfg EleFakePhoton.cfg"
            ;;
        SingleElectronB)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/SingleElectron/Ele27WP80_FT_53_V21_AN4
            cfg="Data.cfg EleFakePhoton.cfg"
            ;;
        SingleElectronC)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012C-22Jan2013-v1/SingleElectron/Ele27WP80_FT_53_V21_AN4
            cfg="Data.cfg EleFakePhoton.cfg"
            ;;
        SingleElectronD)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Run2012D-22Jan2013-v1/SingleElectron/Ele27WP80_FT_53_V21_AN4
            cfg="Data.cfg EleFakePhoton.cfg"
            ;;
        WGToLNuG)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/WGToLNuG_TuneZ2star_8TeV-madgraph-tauola/PU_S10_START53_V7A-v1
            cfg="MC.cfg All.cfg"
            ;;
        WGToLNuG_PtG-30-50)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/WGToLNuG_PtG-30-50_8TeV-madgraph/PU_S10_START53_V7C-v1
            cfg="MC.cfg All.cfg"
            ;;
        WGToLNuG_PtG-50-130)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/WGToLNuG_PtG-50-130_8TeV-madgraph/PU_S10_START53_V7C-v1
            cfg="MC.cfg All.cfg"
            ;;
        WGToLNuG_PtG-130)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/WGToLNuG_PtG-130_8TeV-madgraph-pythia6_tauola/PU_S10_START53_V7A-v1
            cfg="MC.cfg All.cfg"
            ;;
        ZGToLLG)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/ZGToLLG_8TeV-madgraph/PU_S10_START53_V7A-v1
            cfg="MC.cfg All.cfg"
            ;;
        ZGToLLG_PtG-130)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/ZGToLLG_PtG-130_8TeV-madgraph-pythia6_tauola/PU_S10_START53_V7A-v1
            cfg="MC.cfg All.cfg"
            ;;
        WJetsToLNu)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X-PrivateSkim/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/PU_S10_START53_V7A-v2_PtW0To50_Photon25
            cfg="MC.cfg FakePhoton.cfg"
            ;;
        WJetsToLNu_PtW-50To70)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X-PrivateSkim/WJetsToLNu_PtW-50To70_TuneZ2star_8TeV-madgraph/PU_S10_START53_V7A-v1_Photon25
            cfg="MC.cfg FakePhoton.cfg"
            ;;
        WJetsToLNu_PtW-70To100)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X-PrivateSkim/WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph/PU_S10_START53_V7A-v1_Photon25
            cfg="MC.cfg FakePhoton.cfg"
            ;;
        WJetsToLNu_PtW-100)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X-PrivateSkim/WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph/PU_S10_START53_V7A-v1_Photon25
            cfg="MC.cfg FakePhoton.cfg"
            ;;
        DYJetsToLL)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/PU_S10_START53_V7A-v1
            cfg="MC.cfg EleFakePhoton.cfg"
            ;;
        TTJetsSemiLept)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X-PrivateSkim/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/PU_S10_START53_V7C-v1_Photon25
            cfg="MC.cfg All.cfg"
            [ $spec != "nohlt" ] && vetoThreshold="-20"
            ;;
        TTJetsFullLept)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X-PrivateSkim/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/PU_S10_START53_V7C-v2_Photon25
            cfg="MC.cfg All.cfg"
            [ $spec != "nohlt" ] && vetoThreshold="-20"
            ;;
        WW)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/WW_TuneZ2star_8TeV_pythia6_tauola/PU_S10_START53_V7A-v1
            cfg="MC.cfg All.cfg"
            [ $spec != "nohlt" ] && vetoThreshold="-20"
            ;;
        WZ)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/WZ_TuneZ2star_8TeV_pythia6_tauola/PU_S10_START53_V7A-v1
            cfg="MC.cfg All.cfg"
            ;;
        WWJetsTo2L2Nu)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/PU_S10_START53_V7A-v1
            cfg="MC.cfg EleFakePhoton.cfg"
            [ $spec != "nohlt" ] && vetoThreshold="-20"
            ;;
        WWGJets)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/WWGJets_8TeV-madgraph_v2/PU_S10_START53_V7A-v1
            cfg="MC.cfg All.cfg"
            ;;
        TTGJets)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_DR53X/TTGJets_8TeV-madgraph/PU_S10_START53_V19-v1
            cfg="MC.cfg All.cfg"
            ;;
        GJet_Pt-20to40)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Summer12_DR53X/GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/PU_S10_START53_V7A-v1
            cfg="MC.cfg FakeLepton.cfg"
            ;;
        GJet_Pt40)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1/Summer12_DR53X/GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/PU_S10_START53_V7A-v1
            cfg="MC.cfg FakeLepton.cfg"
            ;;
        QCD_Pt-30to40)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v0p1/Summer12_DR53X/QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/PU_S10_START53_V7A-v1
            cfg="MC.cfg FakeLepton.cfg"
            ;;
        QCD_Pt-40)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v0p1/Summer12_DR53X/QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/PU_S10_START53_V7A-v1
            cfg="MC.cfg FakeLepton.cfg"
            ;;
        TChiwg_300)
            sourcedir=/store/RA3Ntuples/TestData/TChiwg_300
            cfg="MC.cfg EleFakePhoton.cfg"
            ;;
    T5wg_400to550)
        sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_FS53/SMS-T5wg_2J_mGl-400to550_mLSP-25to525_TuneZ2star_8TeV-madgraph-tauola/START53_V19_FSIM_PU_S12-v1-sorted
        cfg="MC.cfg All.cfg"
        format='susyEvents_T5wg_{[0-9]+_[0-9]+}_*.root,susyTriggers_T5wg_{}_*.root'
        target='/T5wg'
        reducer='None'
        ;;
    T5wg_600to750)
        sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_FS53/SMS-T5wg_2J_mGl-600to750_mLSP-25to725_TuneZ2star_8TeV-madgraph-tauola/START53_V19_FSIM_PU_S12-v1-sorted
        cfg="MC.cfg All.cfg"
        format='susyEvents_T5wg_{[0-9]+_[0-9]+}_*.root,susyTriggers_T5wg_{}_*.root'
        target='/T5wg'
        reducer='None'
        ;;
        T5wg_800to950)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_FS53/SMS-T5wg_2J_mGl-800to950_mLSP-25to925_TuneZ2star_8TeV-madgraph-tauola/START53_V19_FSIM_PU_S12-v1-sorted
        cfg="MC.cfg All.cfg"
        format='susyEvents_T5wg_{[0-9]+_[0-9]+}_*.root,susyTriggers_T5wg_{}_*.root'
        target='/T5wg'
        reducer='None'
            ;;
        T5wg_1000to1150)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_FS53/SMS-T5wg_2J_mGl-1000to1150_mLSP-25to1125_TuneZ2star_8TeV-madgraph-tauola/START53_V19_FSIM_PU_S12-v1-sorted
            cfg="MC.cfg All.cfg"
            format='susyEvents_T5wg_{[0-9]+_[0-9]+}_*.root,susyTriggers_T5wg_{}_*.root'
            target='/T5wg'
            reducer='None'
            ;;
        T5wg_1200to1300)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_FS53/SMS-T5wg_2J_mGl-1200to1300_mLSP-25to1275_TuneZ2star_8TeV-madgraph-tauola/START53_V19_FSIM_PU_S12-v1-sorted
            cfg="MC.cfg All.cfg"
            format='susyEvents_T5wg_{[0-9]+_[0-9]+}_*.root,susyTriggers_T5wg_{}_*.root'
            target='/T5wg'
            reducer='None'
            ;;
        TChiwg)
            sourcedir=/store/RA3Ntuples/SusyNtuples/cms538v1p2/Summer12_FS53/SMS-TChiwg_2J_mChargino-100to800_TuneZ2star_8TeV-madgraph-tauola/START53_V19_FSIM_PU_S12-v1-sorted
            cfg="MC.cfg All.cfg"
            format='susyEvents_TChiwg_{[0-9]+}_*.root,susyTriggers_TChiwg_{}_*.root'
            target='/TChiwg'
            reducer='None'
            ;;
        *)
            ;;
    esac
}

cd $HOME/src/GammaL/PhotonLepton

for JOB in $JOBS; do
    get_arguments $JOB

    if [ $spec = "nohlt" ]; then
        TARGETDIR=/home/rootd/data/glskim_nohlt$target
        cfg="${cfg} noHLT.cfg"
    else
        TARGETDIR=/home/rootd/data/glskim${spec}${target}
    fi

    ssh ncmu40 "mkdir -p $TARGETDIR 2> /dev/null"

    filesPerJob=$(($(ls $sourcedir | grep susyEvents | wc -l) / 50))
    [ $filesPerJob -lt 10 ] && filesPerJob=10

    cfgstr=""
    for c in $cfg; do
        cfgstr=$cfgstr" $HOME/src/GammaL/PhotonLepton/cfg/$c"
    done

    rm -r skims/$JOB > /dev/null 2>&1
    dap.py -w "skims/$JOB" -a "\"$cfgstr\", $vetoThreshold" -d "ncmu40:$TARGETDIR" -o $JOB.root -u $reducer -f $format -n $filesPerJob -S $sourcedir filter.cc:PhotonLeptonFilter || exit 1
done
