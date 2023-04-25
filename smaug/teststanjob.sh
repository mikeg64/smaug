#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1

source /opt/apps/testapps/el7/software/staging/Lmod/7.3/lmod/7.3/init/bash
module use /opt/apps/testapps/common/modules/all/
module use /opt/apps/testapps/el7/modules/staging/all/

module load CUDA/10.1.243

bin/smaug
