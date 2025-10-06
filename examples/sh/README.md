
## ⚙️ Shell Script Examples

### 💻 Running on a Local PC
```sh
#!/bin/bash

# Activate conda
export PATH="${HOME}/miniconda3/condabin/:${PATH}"
. ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate lcrx

# Run LCRX
input_toml="input.toml"
lcrx run -i ${input_toml}
```

### 🖥️ Running on PBS server
```sh
#!/bin/bash
#PBS -l nodes=1:ppn=16:gpus=1
#PBS -q default
#PBS -l walltime=24:00:00

test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR

# Load required modules (adjust as needed for your environment)
. /usr/share/Modules/init/sh
module load gcc/11.3.0 
module load cuda/12.6
module load gaussian16/c01

# Activate conda
export PATH="${HOME}/miniconda3/condabin/:${PATH}"
. ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate lcrx

# Run LCRX
input_toml="input.toml"
lcrx run -i ${input_toml}
```
