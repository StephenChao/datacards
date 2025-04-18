scripts_path="/home/pku/zhaoyz/Higgs/datacards_lxplus/datacards/scripts/bash/run_combined_1l_0l.sh" # to be replaced with your own path
# GoF of data
source ${scripts_path} -g

# GoF of toys, each run 60 toys, so total 5 * 6 = 300 toys, which will give sufficient statistics
source ${scripts_path} -t --seed 44
source ${scripts_path} -t --seed 45
source ${scripts_path} -t --seed 46
source ${scripts_path} -t --seed 47
source ${scripts_path} -t --seed 48
