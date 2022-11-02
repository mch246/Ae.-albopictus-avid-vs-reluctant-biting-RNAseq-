#!/bin/bash
#SBATCH --job-name=trim --output=FastQC.out
#SBATCH --mail-type=END,FAIL --mail-user=mch284@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs FastQC #
#-----------------------------------------------------------------------------#

0

#- FastQC-------------------------------------------------------------$

/home/mch284/FastQC/fastqc -o  /home/mch284/Auto_biting/mRNA/trim_dir/fastQC  /home/mch284/Auto_biting/mRNA/trim_dir/*fastq.gz


