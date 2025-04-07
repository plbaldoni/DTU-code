#!/bin/bash
#SBATCH --mem=4g
#SBATCH --time=7-00:00:00
#SBATCH --partition=long

sbatch run-1-200.sh

sleep 4h

sbatch run-201-240.sh
