# To clone the base environment

conda create --name simon --clone base

# To activate the new environment called simon

source activate simon

# To see the available environments

conda info --envs

# To install in a given package

conda install -c omnia pymbar --name simon
