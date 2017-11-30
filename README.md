<img src="https://api.travis-ci.org/RUBi-ZA/MODE-TASK.svg?branch=master">
<img src="https://media.readthedocs.org/static/projects/badges/passing.svg">

# MODE-TASK

Collection of tools for analysing normal modes and performing principal component analysis

## Installation

*Download the project:*
```bash
git clone https://github.com/RUBi-ZA/MODE-TASK.git
cd MODE-TASK
```
*Install dependencies and set up Python 2.7 virtual environment:*
```bash
sudo apt install virtualenvwrapper python-dev 
virtualenv venv
source venv/bin/activate
pip install --upgrade pip
pip install numpy
pip install matplotlib
pip install cython
pip install mdtraj
```
*Compile C++ binaries:*
```
sudo apt install g++
g++ -I cpp/src/ ANM.cpp -o ANM
g++ -I cpp/input/ getEigenVectors.cpp -o getEigenVectors

```

## Usage

For more detailed documentation on installation and usage of the tool suite please see our [ReadTheDocs](http://mode-task.readthedocs.io/en/latest/index.html) site

