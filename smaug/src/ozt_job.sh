#!/bin/bash





##$ -j y
##$ -l arch=intel*
##$ -l gpu=1
##$ -l gpu_arch=nvidia-m2070
##$ -l gpu=1,gpu_arch=nvidia-k40m
##$ -P mhd
##$ -N tube1_128
##$ -l mem=12G
##$ -l rmem=12G
##$ -l h_rt=168:00:00







#submit using the command
# sbatch tftestjob.sh
# monitor queue using queue

#SBATCH --account=bdshe04  # Run job under project <project>
#SBATCH --time=1:0:0         # Run for a max of 1 hour

# Node resources:
# (choose between 1-4 gpus per node)

#SBATCH --partition=gpu    # Choose either "gpu" or "infer" node type
#SBATCH --nodes=1          # Resources from a single node
#SBATCH --gres=gpu:1       # One GPU per node (plus 25% of node CPU and RAM per GPU)

#SBATCH --mem=16G
#SBATCH --mail-user=m.griffiths@sheffield.ac.uk

#SBATCH --gpus=1


module load slurm
#module load libs/cuda/6.5.14
#module load gcc/9.3.1/cuda-10.1
#module load cuda/10.1
module load cuda/10.1.243

#cd include
#cp iosmaugparams_ot_1020.h iosmaugparams.h
#cp iosmaugparams_ot_256.h iosmaugparams.h
#cp iosmaugparams_tube1_128.h iosmaugparams.h
#cp iosmaugparams_jet_hydro.h iosmaugparams.h

#cd ..

#cd src
#make ot
#cp usersource_jet_hydro.cu usersource.cu
#cp boundary_jet_hydro.cu boundary.cu
#make clean
#make -f Makefile_k40 smaug
#make smaug
#make -f Makefile_3d smaug

#cd ..

bin/smaug


#Tensor flow test job batch file

#submit using the command
# sbatch tftestjob.sh
# monitor queue using queue

##SBATCH --account=bdshe01  # Run job under project <project>
##SBATCH --time=1:0:0         # Run for a max of 1 hour

# Node resources:
# (choose between 1-4 gpus per node)

##SBATCH --partition=gpu    # Choose either "gpu" or "infer" node type
##SBATCH --nodes=1          # Resources from a single node
##SBATCH --gres=gpu:1       # One GPU per node (plus 25% of node CPU and RAM per GPU)

##SBATCH --mem=16G
##SBATCH --mail-user=m.griffiths@sheffield.ac.uk

##SBATCH --gpus=1

#module load slurm
#module load Anaconda3/2020.02
#source activate tensorflow

nvidia-smi
# Replace my_program with the name of the program you want to run
#python tftest.py

echo test
