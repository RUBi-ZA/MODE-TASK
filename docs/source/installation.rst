Installation
========================================

Platform compatibility
-------------------------------

A Linux-like operating system is recommended. However, MODE-TASK is compatible with most platforms which are able to run Python 2.7 or Python 3.6.

A compiler is required to compile C++ from source, in this instance we use g++.

Download the project
-------------------------------

MODE-TASK can be cloned from its GitHub repository ::

	git clone https://github.com/RUBi-ZA/MODE-TASK.git
	cd MODE-TASK

Installing dependencies 
--------------------------

**Ubuntu:**

Python 2.7 with pip and virtualenv ::

	sudo apt-get update
	sudo apt-get install python-pip virtualenv virtualenvwrapper
	virtualenv venv
	source venv/bin/activate
	pip install -r requirements.txt

Python 3.6 with pip and virtualenv ::

	sudo apt-get update
	sudo apt-get install python3-pip virtualenv virtualenvwrapper
	virtualenv venv
	source venv/bin/activate
	pip3 install -r requirements.txt

Conda ::

	conda create -n mode_task
	source activate mode_task
	conda install -c conda-forge numpy
	conda install -c conda-forge cython
	conda install -c omnia mdtraj
	conda install -c conda-forge scipy
	conda install -c conda-forge pandas
	conda install -c conda-forge sklearn-contrib-lightning
	conda install -c conda-forge matplotlib

To install conda follow their `documentation <https://conda.io/docs/user-guide/install/index.html#regular-installation>`_

**OSX:** ::

	brew update
	brew install python gcc
	pip install virtualenv virtualenvwrapper

**Windows:**

	Enable Windows Subsystem for Linux (WSL) by following `these instructions <https://msdn.microsoft.com/en-us/commandline/wsl/install_guide>`_.

	Install the system dependencies as with Ubuntu above.
