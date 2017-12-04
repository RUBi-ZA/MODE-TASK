<img src="https://api.travis-ci.org/RUBi-ZA/MODE-TASK.svg?branch=master" align="right"><img src="https://media.readthedocs.org/static/projects/badges/passing.svg" align="right"><img src="https://img.shields.io/badge/python-2.7%2C%203.6-blue.svg" align="right">

# MODE-TASK

Collection of tools for analysing normal modes and performing principal component analysis

## Installation

*Download the project:*
```bash
git clone https://github.com/RUBi-ZA/MODE-TASK.git
cd MODE-TASK
```
*Install dependencies and set up Python 2 virtual environment:*
```bash
sudo apt install virtualenvwrapper python-dev python-pip
virtualenv venv
source venv/bin/activate
pip install numpy
pip install cython
pip install mdtraj
pip install scipy
pip install pandas
pip install sklearn
pip install matplotlib
```
*Or optionally with a Python 3 virtual environment:*
```bash
sudo apt install virtualenvwrapper python3-dev python3-pip
virtualenv venv
source venv/bin/activate
pip3 install setuptools
pip3 install numpy
pip3 install cython
pip3 install mdtraj
pip3 install scipy
pip3 install pandas
pip3 install sklearn
pip3 install matplotlib
```
*Compile C++ binaries:*
```
sudo apt install g++
g++ -I cpp/src/ ANM.cpp -o ANM
g++ -I cpp/input/ getEigenVectors.cpp -o getEigenVectors

```

## Usage

For more detailed documentation on installation and usage of the tool suite please see our [ReadTheDocs](http://mode-task.readthedocs.io/en/latest/index.html) site

