#!/bin/bash

# Isolates - Coverage - Replicate - Threads
for i in $(ls 00-RTData); do
	for j in 30 50 80 100 150 200; do
		for k in $(seq 10); do
			for c in 1 4 8 16 32; do
				sbatch -c ${c} sbatch_running_time.sh ${i} ${j} ${k} ${c}
			done
		done
	done
done