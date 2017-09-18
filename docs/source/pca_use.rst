PCA Scripts
====================================

Principal component analysis (PCA) is a useful statistical technique that has found applications in detection of correlated motion in MD data. Protein dynamics is manifested as a change in molecular structure, or conformation over a time scale. PCA extracts most important motions from a protein's MD trajectory using a covariance/correlation matrix (C-matrix) constructed from atomic coordinates. Different types of coordinate systems (Cartesian or internal coordinates) can be employed to define atomic movement in each time frame of a trajectory. Modes describing the protein motions can be constructed by diagonalizing the C-matrix. It leads to a complete set of orthonormal (orthogonal and normalized) collective modes (eigenvectors) with eigenvalues (variance) that characterize the protein dynamics. The largest eigenvalues represents the most collective spatial motion. When the original mean centered data (MD trajectory) is projected on the eigenvectors, the results are called Principle Components (PC). Diagonalization of C-matrix could be done by Eigenvalue decomposition (EVD) or Singular value decomposition (SVD), with the latter being computationally efficient.  

As stated earlier, different representations of protein conformations can be used. One can choose Cartesian coordinates or internal coordinates such as the pairwise distance between atoms, 1-3 angle, torsional angles (psi or phi). Since decomposition of a C-matrix is memory intensive and very often the program will run out of memory, often a coarse graining is required such as selecting CA atoms. The user can select the subset of atoms from the trajectory for the analysis such as CA, backbone atoms or all protein's atoms. It is highly recommended that the user should strip the water from the trajectory before hand, as it would result in faster loading and alleviate the memory issues. 

PCA uses linear transformation which may not be sufficient in cases where variables are non-linearly related. Thus, the user has the option to perform Nonlinear generalization of PCA such as Kernel PCA (kPCA). Caution should be given while interpreting the kPCA results since it is mapped to a feature space which is inherently different than conformational space. Nevertheless, kPCA is useful in understanding the protein's functions in terms of its conformational dynamics.   


**General Usage:** 

To perform PCA on a protein's MD trajectory we need a sufficiently sampled MD trajectory and a corresponding topology file. This can be achieved by running the following command.  

**Command:** 
	``pca.py -t <MD trajectory> -p <topology file>``	

To see the all the available options run the following command: 
	``pca.py -h``

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Trajectory file *      | File       |``-t``              | MD trajectory input file    |
|                        |            |                    | (.xtc, .mdcrd etc.)         |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Topology file *        | File       |``-p``              | Topology file               |
|                        |            |                    | (.gro, .pdb etc)            |
+------------------------+------------+--------------------+-----------------------------+
| Output directory       | String     |``-out``            | Name of the output directory|
|                        |            |                    | .Default is out suffixed by |
|                        |            |                    | trajectory name             |
+------------------------+------------+--------------------+-----------------------------+
| Atom group             | String     |``-ag``             | Group of atom for PCA.      |
|                        |            |                    | Default is CA atoms.   |
| 			 |	      | 		   | Other options are:          |
|                        |            |                    | all= all atoms,             |
|                        |            |                    | backbone = backbone atoms,  |
|                        |            |                    | CA= C-alpha atoms,          |
|                        |            |                    | protein= protein's atoms    |
+------------------------+------------+--------------------+-----------------------------+
| Reference structure    | File       | ``-r``             | Reference structure for RMSD|
|                        |            |                    | Default: First frame of MD  |
|                        |            |                    | trajectory                  |
+------------------------+------------+--------------------+-----------------------------+
| PCA method             | String     | ``-pt``            | PCA method.                 |
|                        |            |                    | Default is svd (Single Value|
|                        |            |                    | Decomposition) PCA.         |
|                        |            |                    | Options are: evd, kpca, svd,|
|                        |            |                    | ipca. If svd is selected,   |
|                        |            |                    | additional arguments can be |
|                        |            |                    | passed by flag -st.         |
|                        |            |                    | If KernelPCA is selected    |
|                        |            |                    | kernel type can also be     |
|                        |            |                    | defined by flag -k          |
+------------------------+------------+--------------------+-----------------------------+
| Number of components   | Int        | ``-nc``		   | Number of components to keep|
|                        |            |                    | in a PCA object.            |
|                        |            |                    | Default all the components  |
|                        |            |                    | will be kept.               |
+------------------------+------------+--------------------+-----------------------------+
| Kernel Type            | String     | ``-kt``            | Type of kernel for          |
|                        |            |                    | KernalPCA.                  |
|                        |            |                    | default is linear.          |
|                        |            |                    | Options are :linear, poly,  |
|                        |            |                    | rbf, sigmoid, cosine,       |
|                        |            |                    | precomputed                 |
+------------------------+------------+--------------------+-----------------------------+
| SVD solver type        | String     | ``-st``            | Type of svd_solver for SVD  |
|                        |            |                    | (Single Value Decomposition)|
|                        |            |                    | PCA. Default is auto.       |
|                        |            |                    | Options are: auto, full,    |
|                        |            |                    | arpack, randomized          |
+------------------------+------------+--------------------+-----------------------------+
 
**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PC plots               | 2D Plot of first 3 PCs. It  |
|                        | is grace formatted text file|
+------------------------+-----------------------------+
| PC plots (.png)	 | 2D Plot of first 3 PCs. Same|
|                        | as above, but point are     |
|                        | color coded according to MD |
|                        | time                        |
+------------------------+-----------------------------+
| Scree plot 	         | Scree plot of contribution  |
|                        | of first 100 modes          |
|                        | (eigenvectors)              | 
+------------------------+-----------------------------+
| RMSD plot              | RMSD of selected atoms over |
|                        | the MD time                 |
+------------------------+-----------------------------+
| RMSD Modes             | Plot of contribution of     |
|                        | each resdiues toward first 3|
|                        | Mode (eigenvectors)         |
+------------------------+-----------------------------+

Beside the above-mentioned plots, it also prints useful information on the terminal such as, information about the trajectory, Kaiser-Meyer-Olkein (KMO) index of the trajectory, and cosine contents of the first few PCs. KMO value range from 0 to 1, 1 indicating that the MD has been sampled sufficiently. The cosine content of PCA projections can be used as an indicator if a simulation is converged. Squared cosine value should be more than 0.5.   


**Specific Examples:**

PCA on Cartesian coordinates
-------------------------------

Given a trajectory called ``trajectory.xtc`` and a topology file called ``complex.pdb``, the following command is used:

	``pca.py -t trajectory.xtc -p complex.pdb``

This will perform the singular value decomposition (SVD) based PCA on CA atoms by default. To use other method, see the following examples.


**SVD PCA**
^^^^^^^^^^^^^^^

To perform SVD PCA on CA atoms of a MD trajectory

**Command:** 
	``pca.py -t trajectory.xtc -p complex.pdb -ag CA -pt svd``

To perform the SVD PCA on backbone atoms

**Command:** 
	``pca.py -t trajectory.xtc -p complex.pdb -ag backbone -pt svd``



**Kernel PCA**
^^^^^^^^^^^^^^^^^^

To perform the Kernel PCA with linear kernel

**Command:** 
	``pca.py -t trajectory.xtc -p complex.pdb -ag CA -pt kpca -kt linear``

To perform the Kernel PCA with rbf kernel

**Command:** 
	``pca.py -t trajectory.xtc -p complex.pdb -ag CA -pt kpca -kt rbf``

**Incremental PCA** 
^^^^^^^^^^^^^^^^^^^^^^^

Incremental PCA (IPCA) is a variant of usual PCA, which uses low-rank approximation of the input MD trajectory. It uses the amount of memory to store the input trajectory which is independent of trajectory size. IPCA is very useful in case the size of trajectory is bigger than the available computer memory.

**Command:** 
	  ``pca.py -t trajectory.xtc -p complex.pdb -ag CA -pt ipca``

**Eigenvalue decomposition (EVD) PCA** 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To perform the PCA by eigenvalue decomposition

**Command:** 
	``pca.py -t trajectory.xtc -p complex.pdb -ag CA -pt evd``

**Detailed usage:** 

Run the following command to see the detailed usage and other options:
	``pca.py -h``


PCA on internal coordinates
-----------------------------

User can also perform the PCA on internal coordinates of a MD trajectory. Options are available for different types of internal coordinates such as:*pairwise distance between atoms*, *1-3 angle between backbone atoms*, *psi angle*, and *phi angle*.  

**General Usage:**

**Command:** 
	``internal_pca.py -t <MD trajectory> -p <topology file>``

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Trajectory file *      | File       |``-t``              | MD trajectory input file    |
|                        |            |                    | (.xtc, .mdcrd etc.)         |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Topology file *        | File       |``-p``              | Topology file               |
|                        |            |                    | (.gro, .pdb etc)            |
+------------------------+------------+--------------------+-----------------------------+
| Output directory       | String     |``-out``            | Name of the output directory|
|                        |            |                    | . Default is out suffixed by|
|                        |            |                    | trajectory name             |
+------------------------+------------+--------------------+-----------------------------+
| Atom group             | String     |``-ag``             | Group of atom for PCA.      |
|                        |            |                    | Default is CA atoms.   |
| 			 |	      | 		   | Other options are:          |
|                        |            |                    | all= all atoms,             |
|                        |            |                    | backbone = backbone atoms,  |
|                        |            |                    | CA= C-alpha atoms,          |
|                        |            |                    | protein= protein's atoms    |
+------------------------+------------+--------------------+-----------------------------+
| Coordinate Type        | String     | ``-ct``            | Internal cordinate type.    |
|                        |            |                    | Options are: distance,      |
|                        |            |                    | angles, phi and, psi        |
+------------------------+------------+--------------------+-----------------------------+

 
**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PC plots               | 2D Plot of first 3 PCs. It  |
|                        | is grace formatted text file|
+------------------------+-----------------------------+
| PC plots (.png)	 | 2D Plot of first 3 PCs. Same|
|                        | as above, but points are    |
|                        | color coded according to MD |
|                        | time                        |
+------------------------+-----------------------------+
| Scree plot 	         | Scree plot of contriution   |
|                        | of first 100 modes          |
|                        | (eigenvectors)              | 
+------------------------+-----------------------------+

**Specific Examples:** 

**PCA on pairwise distance between CA atoms:** 

To perform the PCA on pairwise distance between CA atoms of an MD trajectory ``trajectory.xtc`` and a topology file ``complex.pdb``

**Command:** 
	``internal_pca.py -t trajectory.xtc -p complex.pdb -ag CA -ct distance``	

**PCA on psi angles:** 

**Command:** 
	``internal_pca.py -t trajectory.xtc -p complex.pdb -ct psi``

**Detailed usage:** 

Run the following command to see the detailed usage and other options:
	``internal_pca.py -h``

MDS (Multi-dimentional scaling)  on MD trajectory
---------------------------------------------------

MDS is a tool to visualize the similarity or dissimilarity in a dataset. Two types of dissimilarity measures can be used in the case of a MD trajectory. The first is Euclidean distance between internal coordinates of a protein structure, the second is pairwise RMSD between a set of atoms over the frames of a MD trajectory. 

**General Usage:**

**command:**
	``mds.py -t <MD trajectory> -p <topology file>``

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Trajectory file *      | File       |``-t``              | MD trajectory input file    |
|                        |            |                    | (.xtc, .mdcrd etc.)         |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Topology file *        | File       |``-p``              | Topology file               |
|                        |            |                    | (.gro, .pdb etc)            |
+------------------------+------------+--------------------+-----------------------------+
| Output directory       | String     |``-out``            | Name of the output directory|
|                        |            |                    | . Default is out suffixed by|
|                        |            |                    | trajectory name             |
+------------------------+------------+--------------------+-----------------------------+
| Atom group             | String     |``-ag``             | Group of atom for MDS.      |
|                        |            |                    | Default is CA atoms.   |
|                        |            |                    | Other options are:          |
|                        |            |                    | all= all atoms,             |
|                        |            |                    | backbone = backbone atoms,  |
|                        |            |                    | CA= C-alpha atoms,          |
|                        |            |                    | protein= protein's atoms    |
+------------------------+------------+--------------------+-----------------------------+
| MDS type               | String     | ``-mt``            | Type of MDS. Options are    |
|                        |            |                    | nm=non-metric, metric=metric|
+------------------------+------------+--------------------+-----------------------------+
| Dissimilarity type     | String     | ``-dt``            | Type of dissimilarity matrix|
|                        |            |                    | to use. euc = Euclidean     |
|                        |            |                    | distance between internal   |
|                        |            |                    | coordinates, rmsd= pairwise |
|                        |            |                    | RMSD. Default is rmsd       |
+------------------------+------------+--------------------+-----------------------------+
| Coordinate type        | String     | ``-ct``            | Internal coordinates type.  |
|                        |            |                    | Default is pairwise distance|
|                        |            |                    | Only used if Dissimilarity  |
|                        |            |                    | type is euclidean           |
+------------------------+------------+--------------------+-----------------------------+
| Atom indices           | String     | ``-ai``            | Group of atom for pairwise  |
|                        |            |                    | distance. Default is CA|
|                        |            |                    | atoms. Other options are:   |
|                        |            |                    | all= all atoms,backbone =   |
|                        |            |                    | backbone atoms, alpha=      |
|                        |            |                    | C-alpha atoms,heavy= all non|
|                        |            |                    | hydrogen atoms, minimal=CA, |
|                        |            |                    | CB,C,N,O atoms              |
+------------------------+------------+--------------------+-----------------------------+

 
**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PC plots               | 2D Plot of the first 3 PCs. |
|                        | It                          |
|                        | is grace formatted text file|
+------------------------+-----------------------------+
| PC plots (.png)        | 2D Plot of the first 3 PCs. |
|                        | Same                        |
|                        | as above, but points are    |
|                        | color coded according to MD |
|                        | time                        |
+------------------------+-----------------------------+

**Specific Examples:**

**MDS on pairwise RMSD:**  

To perform MDS on the pairwise RMSD between CA atoms
	
**Command:** 
	``mds.py -t trajectory.xtc -p complex.pdb -dt rmsd -ag CA``

**MDS on internal coordinates:**  

To perform MDS on the pairwise distance between CA atoms 

**Command:** 
	``mds.py -t trajectory.xtc -p complex.pdb -dt euc -ag CA``

**Detailed usage:** 

Run the following command to see the detailed usage and other options:
	``mds.py -h``

t-SNE on MD trajectory
--------------------------------------------------------------------

t-distributed Stochastic Neighbor Embedding (t-SNE) is a tool for dimensionality reduction. It is a variant of stochastic  neighbor embedding technique. t-SNE uses a measure of dissimilarity, which in case of MD trajectory may be the Euclidean distance between internal coordinates or pairwise RMSD.    


**General Usage:**

**Command:**
	``tsne.py -t <MD trajectory> -p <topology file>``

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Trajectory file *      | File       |``-t``              | MD trajectory input file    |
|                        |            |                    | (.xtc, .mdcrd etc.)         |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Topology file *        | File       |``-p``              | Topology file               |
|                        |            |                    | (.gro, .pdb etc)            |
+------------------------+------------+--------------------+-----------------------------+
| Output directory       | String     |``-out``            | Name of the output directory|
|                        |            |                    | . Default is out suffixed by|
|                        |            |                    | trajectory name             |
+------------------------+------------+--------------------+-----------------------------+
| Atom group             | String     |``-ag``             | Group of atom for t-SNE.    |
|                        |            |                    | Default is CA atoms.   |
| 			 |	      | 		   | Other options are:          |
|                        |            |                    | all= all atoms,             |
|                        |            |                    | backbone = backbone atoms,  |
|                        |            |                    | CA= C-alpha atoms,          |
|                        |            |                    | protein= protein's atoms    |
+------------------------+------------+--------------------+-----------------------------+
| Coordinate type        | String     | ``-ct``            | Internal coordinates type.  |
|                        |            |                    | Default is pairwise distance|
|                        |            |                    | . Only used if Dissimilarity|
|                        |            |                    | type is euclidean           |
+------------------------+------------+--------------------+-----------------------------+
| Dissimilarity type     | String     | ``-dt``            | Type of dissimilarity matrix|
|                        |            |                    | to use. euc = Euclidean     |
|                        |            |                    | distance between internal   |
|                        |            |                    | coordinates, rmsd= pairwise |
|                        |            |                    | RMSD. Default is rmsd       |
+------------------------+------------+--------------------+-----------------------------+
| Atom indices           | String     | ``-ai``            | Group of atom for pairwise  |
|                        |            |                    | distance. Default is CA|
|                        |            |                    | atoms. Other options are:   |
|                        |            |                    | all= all atoms,backbone =   |
|                        |            |                    | backbone atoms, alpha=      |
|                        |            |                    | C-alpha atoms,heavy= all non|
|                        |            |                    | hydrogen atoms, minimal=CA, |
|                        |            |                    | CB,C,N,O atoms              |
+------------------------+------------+--------------------+-----------------------------+

 
**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PC plots               | 2D Plot of the first 3 PCs. |
|                        | It                          |
|                        | is grace formatted text file|
+------------------------+-----------------------------+
| PC plots (.png)	 | 2D Plot of the first 3 PCs. |
|                        | Same                        |
|                        | as above, but point are     |
|                        | color coded according to MD |
|                        | time                        |
+------------------------+-----------------------------+

**Specific Example:**

**t-SNE on CA atoms:**
To perform t-SNE using the pairwise RMSD of CA atoms as index of dissimilarity.

**command:**
	``tsne.py -t trajectory.xtc -p complex.pdb -ag CA -dt rmsd``

To perform t-SNE using the Euclidean space between pairwise distance of CA atoms as index of dissimilarity. 

**command:**
	``tsne.py -t trajectory.xtc -p complex.pdb -ag CA -dt euc -ai alpha``

**Detailed usage:**

Run the following command to see the detailed usage and other options:
	``tsne.py -h``
