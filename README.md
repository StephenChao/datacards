# datacards
Combine datacards for CMS boosted HWW analysis

## Setup Combine environment

run `setup_env.sh` to set up the CMSSSW & Combine environment, the Combine & Rhalphabet package are already integrated in the `root://cmseos.fnal.gov//store/user/rkansal/CMSSW_11_3_4.tgz`, thanks to @rkansal47 for this!

Since 0lepton datacards will need this Rhalphabet package to run, so this special environment is necessary. Note that an `slc7` cluster is needed.

## Data-cards for each channel

under `./cards`, there are datacards for each channel

## Scripts to run each channel's datacards

under `./scripts/bash`, there are bash scripts to run the Combine for each channel

E.g., run `source YOUR_PATH_TO/scripts/bash/run_0l.sh -wbsil` under `YOUR_PATH_TO/cards/0lepton` would run the workspace, NLL scan, significance, limits and initial fit for impacts.

## Run impact 
E.g.,
use `python3 -u YOUR_PATH_TO/datacards/scripts/python/run_impacts_locally.py --impacts` under `YOUR_PATH_TO/cards/0lepton` to run the impacts,
use `python3 -u YOUR_PATH_TO/datacards/scripts/python/run_impacts_locally.py --collect` under `YOUR_PATH_TO/cards/0lepton` to collect the impacts.

## Postfit plots
Use `./postfit/makeplots.ipynb`, currently only 0lepton & 1lepton postfit plots after combination is supported.