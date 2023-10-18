# MDM-TASK-web

**MDM-TASK-web** is web server for MD-TASK and MODE-TASK and can be accessed here: [https://mdmtaskweb.rubi.ru.ac.za](https://mdmtaskweb.rubi.ru.ac.za)

If you use MDM-TASK-web in your work, please cite the following:
 - [MDM-TASK-web: MD-TASK and MODE-TASK web server for analyzing protein dynamics](https://doi.org/10.1016/j.csbj.2021.08.043)
 - [MODE-TASK: large-scale protein motion tools](https://academic.oup.com/bioinformatics/article/34/21/3759/5021681)

## Installation
Clone the MDM-TASK-web branch
```
git clone -b mdmtask https://github.com/RUBi-ZA/MODE-TASK
```
Go inside the cloned MODE-TASK folder and run the following to create the working environment named "mdmtaskweb"
```
conda env create -f environment.yml
```
Activate the conda environment
```
conda activate mdmtaskweb
```

In order to use ANM.py and getEigenVectors.py, ANM.cpp needs to be compiled. The following may be used to do so, from the MODE-TASK directory:
```
sudo apt install g++
g++ -I src/cpp/src/ src/ANM.cpp -o src/ANM

```
