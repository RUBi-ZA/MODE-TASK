Installation
========================================

Platform compatibility
-------------------------------

MODE-TASK is compatible with most platforms which are able to run Python 2.7 and g++.


Install system dependencies
-----------------------------

**Ubuntu 14.04:** ::
	
	sudo apt-get update
	sudo apt-get install python-dev virtualenv virtualenvwrapper g++

**OSX:** ::

	brew update
	brew install python gcc
	pip install virtualenv virtualenvwrapper

**Windows:**

	Enable Windows Subsystem for Linux (WSL) by following `these instructions <https://msdn.microsoft.com/en-us/commandline/wsl/install_guide>`_.

	Install the system dependencies as with Ubuntu above.


Install Python dependencies
--------------------------------

It is recommended to create a Python virtual environment for installing and managing dependencies::

	virtualenv venv
	source venv/bin/activate
	pip install --upgrade pip
	pip install numpy
	pip install matplotlib
	pip install cython
	pip install mdtraj


Download the project
-------------------------------

MODE-TASK can be cloned from it's GitHub repository ::

	git clone https://github.com/RUBi-ZA/NMA-TASK.git
	cd NMA-TASK

Activate the virtual environment you created in the previous step when using MODE-TASK with::

	source venv/bin/activate
