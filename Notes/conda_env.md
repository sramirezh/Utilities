# To clone the base environment

conda create --name simon --clone base

# To activate the new environment called simon

source activate simon

# To see the available environments

conda info --envs

# To install in a given package

conda install -c omnia pymbar --name simon

# To export an environment

conda env export > environment.yml

# To recreate it 

conda env create -f environment.yml

# To remove an environment

conda remove --name myenv --all

# when spyder does not start
conda update python spyder

# To install cloned packages 

source activate TARGET_ENVIRONMENT
(Below and example from https://stackoverflow.com/questions/22312665/install-python-packages-to-correct-anaconda-environment)
git clone git://github.com/pudo/dataset.git
pip install ./dataset.  (after conda install pip)

# To install the requirements file in the present environment 

conda install --file requirements.txt
