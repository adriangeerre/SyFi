#!/bin/bash
#SBATCH --account CCRP_Data
#SBATCH --partition short,normal
#SBATCH --mem 8G
#SBATCH -t 01:00:00

# Running time
#-------------

# Direct variables
iso=${1}
cov=${2}
rep=${3}
threads=${4}

# Define values
if [ ${iso} == "Marinithermus" ]; then
	glen=2269167; nt=25943791; t=583; rlen=76
elif [ ${iso} == "Marinomonas" ]; then
	glen=3899940; nt=29132357; t=238; rlen=76
elif [ ${iso} == "Olsenella" ]; then
	glen=2051896; nt=19937742; t=1776; rlen=76
elif [ ${iso} == "Rhizobium" ]; then
	glen=4854518; nt=36782323; t=501; rlen=76
fi

# Exec
python running_time.py -i ${iso} -g ${glen} -a ${nt} -b ${t} -l ${rlen} -c ${cov} -r ${rep} -t ${threads}
