#!/bin/bash
### adapted from https://github.com/rkansal47/HHbbVV/blob/main/src/HHbbVV/combine/run_blinded.sh
### Author: Raghav Kansal, Yuzhe Zhao, Farouk Mokhtar

####################################################################################################
# Script for fits
#
# 1) Combines cards and makes a workspace (--workspace / -w)
# 2) Background-only fit (--bfit / -b)
# 3) Expected asymptotic limits (--limits / -l)
# 4) Expected significance (--significance / -s)
# 5) Fit diagnostics (--dfit / -d)
# 6) GoF on data (--gofdata / -g)
# 7) GoF on toys (--goftoys / -t),
# 8) Impacts: initial fit (--impactsi / -i), per-nuisance fits (--impactsf $nuisance), collect (--impactsc $nuisances)
# 9) Bias test: run a bias test on toys (using post-fit nuisances) with expected signal strength
#    given by --bias X.
#
# Specify seed with --seed (default 42) and number of toys with --numtoys (default 100)
#
# Usage ./run_1l.sh [-wblsdgt] [--numtoys 100] [--seed 42]

#TO BE NOTICED: should only be run in the datacards directory
####################################################################################################


####################################################################################################
# Read options
####################################################################################################


workspace=0
bfit=0
limits=0
significance=0
dfit=0
dfit_asimov=0
gofdata=0
goftoys=0
impactsi=0
impactsf=0
impactsc=0
seed=444
numtoys=100
bias=-1
mintol=0.5 # --cminDefaultMinimizerTolerance
# maxcalls=1000000000  # --X-rtd MINIMIZER_MaxCalls

options=$(getopt -o "wblsdrgti" --long "workspace,bfit,limits,significance,dfit,dfitasimov,resonant,gofdata,goftoys,impactsi,impactsf:,impactsc:,bias:,seed:,numtoys:,mintol:" -- "$@")
eval set -- "$options"

while true; do
    case "$1" in
        -w|--workspace)
            workspace=1
            ;;
        -b|--bfit)
            bfit=1
            ;;
        -l|--limits)
            limits=1
            ;;
        -s|--significance)
            significance=1
            ;;
        -d|--dfit)
            dfit=1
            ;;
        --dfitasimov)
            dfit_asimov=1
            ;;
        -g|--gofdata)
            gofdata=1
            ;;
        -t|--goftoys)
            goftoys=1
            ;;
        -i|--impactsi)
            impactsi=1
            ;;
        --impactsf)
            shift
            impactsf=$1
            ;;
        --impactsc)
            shift
            impactsc=$1
            ;;
        --seed)
            shift
            seed=$1
            ;;
        --numtoys)
            shift
            numtoys=$1
            ;;
        --mintol)
            shift
            mintol=$1
            ;;
        --bias)
            shift
            bias=$1
            ;;
        --)
            shift
            break;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
    shift
done

echo "Arguments: workspace=$workspace bfit=$bfit limits=$limits \
significance=$significance dfit=$dfit gofdata=$gofdata goftoys=$goftoys \
seed=$seed numtoys=$numtoys"



####################################################################################################
# Set up fit arguments
#
# We use channel masking to "mask" the blinded and "unblinded" regions in the same workspace.
# (mask = 1 means the channel is turned off)
####################################################################################################

dataset=data_obs
cards_dir="."
ws=${cards_dir}/combined
wsm=${ws}_withmasks
wsm_snapshot=higgsCombineSnapshot.MultiDimFit.mH125

CMS_PARAMS_LABEL="CMS_HWW_boosted"

outsdir=${cards_dir}/outs
mkdir -p $outsdir

# ####################################################################################################
# # Combine cards, text2workspace, fit, limits, significances, fitdiagnositcs, GoFs
# ####################################################################################################
# # # need to run this for large # of nuisances
# # # https://cms-talk.web.cern.ch/t/segmentation-fault-in-combine/20735

# ADD REGIONS
VBF="VBF"
ggFpt250to350="ggFpt250to350"
ggFpt350to500="ggFpt350to500"
ggFpt500toInf="ggFpt500toInf"
TopCR="TopCR"
WJetsCR="WJetsCR"

########################### define SR/CR with datacards

ccargs_1l="VBF=${cards_dir}/${VBF}.txt ggFpt250to350=${cards_dir}/${ggFpt250to350}.txt ggFpt350to500=${cards_dir}/${ggFpt350to500}.txt ggFpt500toInf=${cards_dir}/${ggFpt500toInf}.txt TopCR=${cards_dir}/${TopCR}.txt WJetsCR=${cards_dir}/${WJetsCR}.txt"

if [ $workspace = 1 ]; then
    echo "Combining cards:"
    for file in $ccargs_1l; do
    echo "  ${file##*/}"
    done
    echo "-------------------------"
    combineCards.py $ccargs_1l > $ws.txt
    echo "Running text2workspace"
    text2workspace.py $ws.txt --channel-masks -o $wsm.root 2>&1 | tee $outsdir/text2workspace.txt
else
    if [ ! -f "$wsm.root" ]; then
        echo "Workspace doesn't exist! Use the -w|--workspace option to make workspace first"
        exit 1
    fi
fi


if [ $significance = 1 ]; then
    echo "Expected significance"

    # combine -M Significance -d $ws -m 125 -t -1 --expectSignal=1 --rMin -1 --rMax 5

    combine -M Significance -d ${wsm}.root -t -1 --expectSignal 1 2>&1 | tee $outsdir/ExpectedSignificance.txt
    echo "Observed significance"

    combine -M Significance -d ${wsm}.root  2>&1 | tee $outsdir/ObservedSignificance.txt

    echo "Expected NLL scan, split to stat. + syst."
    ### Expected breakdown
    combine -M MultiDimFit -t -1 --expectSignal 1 -m 125 -d ${wsm}.root --rMax 5 --rMin 0 \
    --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1  \
    -n SnapshotExpected --algo grid --points 100  2>&1 | tee $outsdir/ExpectedScanMultiDimFit.txt
    combine -M MultiDimFit -t -1 --expectSignal 1 -m 125 -d ${wsm}.root --rMax 5 --rMin 0 \
    --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1  \
    -n BestfitSnapshotExpected --saveWorkspace 2>&1 | tee $outsdir/ExpectedBestfitMultiDimFit.txt
    combine -M MultiDimFit -t -1 --expectSignal 1 higgsCombineBestfitSnapshotExpected.MultiDimFit.mH125.root -n ExpectedfreezeAll \
    -m 125 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
    --rMin 0 --rMax 5 --algo grid --points 100 \
    --freezeParameters allConstrainedNuisances --snapshotName MultiDimFit 2>&1 | tee $outsdir/ExpectedBreakdownMultiDimFit.txt
    
    plot1DScan.py higgsCombineSnapshotExpected.MultiDimFit.mH125.root --main-label "With systematics" --main-color 1 --others higgsCombineExpectedfreezeAll.MultiDimFit.mH125.root:"Stat-only":2 -o ExpectedBreakdown --breakdown Syst,Stat


    echo "Observed NLL scan, split to stat. + syst."
    ### Observed breakdown
    combine -M MultiDimFit  -m 125 -d ${wsm}.root --rMax 5 --rMin -5 \
    --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01  \
    -n SnapshotObserved --algo grid --points 100  2>&1 | tee $outsdir/ObservedMultiDimFit.txt
    combine -M MultiDimFit  -m 125 -d ${wsm}.root --rMax 5 --rMin -5 \
    --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
    -n BestfitSnapshotObserved --saveWorkspace 2>&1 | tee $outsdir/ObservedBestfitMultiDimFit.txt

    combine -M MultiDimFit higgsCombineBestfitSnapshotObserved.MultiDimFit.mH125.root -n ObservedfreezeAll \
    -m 125 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.03 \
    --rMin -5 --rMax 5 --algo grid --points 100 \
    --freezeParameters allConstrainedNuisances --snapshotName MultiDimFit 2>&1 | tee $outsdir/ObservedBreakdownMultiDimFit.txt
    
    plot1DScan.py higgsCombineSnapshotObserved.MultiDimFit.mH125.root --main-label "With systematics" --main-color 1 --others higgsCombineObservedfreezeAll.MultiDimFit.mH125.root:"Stat-only":2 -o ObservedBreakdown --breakdown Syst,Stat

fi

if [ $dfit_asimov = 1 ]; then

    echo "Fit Diagnostics"
    combine -M FitDiagnostics -m 125 -d $wsm.root \
    -t -1 --expectSignal=1 --saveWorkspace --saveToys -n Asimov --ignoreCovWarning \
    --saveShapes --saveNormalizations --saveWithUncertainties --saveOverallShapes 2>&1 | tee $outsdir/FitDiagnostics.txt

    # python diffNuisances.py fitDiagnostfitDiagnosticsAsimov.root --abs
    combineTool.py -M ModifyDataSet ${wsm}.root:w ${wsm}_asimov.root:w:toy_asimov -d higgsCombineAsimov.FitDiagnostics.mH125.123456.root:toys/toy_asimov

    echo "Fit Shapes"
    PostFitShapesFromWorkspace --dataset toy_asimov -w ${wsm}_asimov.root --output FitShapesAsimov.root \
    -m 125 -f fitDiagnosticsAsimov.root:fit_b --postfit --print 2>&1 | tee $outsdir/FitShapesAsimov.txt
fi


if [ $limits = 1 ]; then
    # echo "Expected limits"
    # combine -M AsymptoticLimits -m 125 -n "" -d ${wsm_snapshot}.root --snapshotName MultiDimFit -v 1 \
    # --saveWorkspace --saveToys --bypassFrequentistFit -s $seed \
    # --floatParameters r --toysFrequentist --run blind 2>&1 | tee $logsdir/AsymptoticLimits.txt

    combine -M AsymptoticLimits --run expected -d $wsm.root -t -1  -v 1 --expectSignal 1 2>&1 | tee $outsdir/AsymptoticLimits.txt
fi


if [ $impactsi = 1 ]; then

    echo "Initial fit for impacts"
    combineTool.py -M Impacts -d $wsm.root -t -1 --rMin -1 --rMax 2 -m 125 --robustFit 1 --doInitialFit --expectSignal 1 2>&1 | tee $outsdir/impact1.txt
    combineTool.py -M Impacts -d $wsm.root -t -1 --rMin -1 --rMax 2 -m 125 --robustFit 1 --doFits --expectSignal 1 --parallel 50 2>&1 | tee $outsdir/impact2.txt
    combineTool.py -M Impacts -d $wsm.root -t -1 --rMin -1 --rMax 2 -m 125 --robustFit 1 --output impacts.json --expectSignal 1 2>&1 | tee $outsdir/impact3.txt
    plotImpacts.py -i impacts.json -o impacts

fi

