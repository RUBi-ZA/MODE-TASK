NMA Scripts
====================================

Coarse grain
-------------------------------

Takes a protein structure or biological assembly and coarse grains to select a set amount of CB atoms
**Command:** ::
	
	coarseGrain.py <options> --pdbFile <pdb file> --atomType <string>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *        	 | File       |``--pdb``           | PDB structure to coarse     |
|                        |            |                    | grain. Can also accept      |
|                        |            |                    | biological assembly (.pdb1) |
+------------------------+------------+--------------------+-----------------------------+
| Atom type *            | String     |``--atomType``      | Specify the type of atom to |
|                        |            |                    | be selected in CG models.   |
|                        |            |                    | Only CA or CB accepted.     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Course grain level     | Comma      |``--cg``            | Level/Levels by which to    |
|                        | Separted   |                    | coarse grain protein.       |
|                        | String     |                    | Lower is less coarse grained|
|                        |            |                    | is less coarse grained.     |
|                        |            |                    | E.g --cg 4                  |
|                        |            |                    |     OR                      |
|                        |            |                    | E.g --cg 1,3,5              |
|                        |            |                    | Default: 4                  |
+------------------------+------------+--------------------+-----------------------------+
| Starting atom          | Integer    |``--startingAtom``  | Residue number of the    	 |
|                        |            |                    | starting atom.              |
|                        |            |                    | Default: 1                  |
+------------------------+------------+--------------------+-----------------------------+
| Output file            | File       |``--output``        | Specify a name for the PDB	 |
|                        |            |                    | output file.                |
|                        |            |                    | Default: ComplexCG.pdb      |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PDB file/s             | Coarse grained protein/s    |
|                        | or macromolecule/s          |
+------------------------+-----------------------------+

ANM
-------------------------------

Construct an elastic network model of a protein complex and solves for the eigenvalues and eigenvectors of the system. 

**Compile:** ::

    g++ -I cpp/src/ ANM.cpp -o ANM

**Command:** ::

	ANM <options> --pdb <pdb file> --atomType <atom type>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *             | File       |``--pdb``           | PDB input file              |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Atom type *            | String     |``--atomType``      | Specify the type of atom to |
|                        |            |                    | be selected in CG models.   |
|                        |            |                    | Only CA or CB accepted.     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Cutoff                 | Integer    |``--cutoff``        | Cuttoff radius in Å.        |
|                        |            |                    | Default: 15                 |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| W matrix               | Text file of :math:`3N`     |
|                        | eigenvalues                 |
+------------------------+-----------------------------+
| VT matrix              | Text file of :math:`3N`\ x\ |
|                        | :math:`3N` eigenvectors.    |
|                        | Printed in rows             |
+------------------------+-----------------------------+
| U matrix               | Text file of :math:`3N`\ x\ |
|                        | :math:`3N` eigenvectors.    |
|                        | Printed in columns          |
+------------------------+-----------------------------+

Get eigenvectors
-------------------------------

**Compile:** ::

	g++ -I cpp/input/ getEigenVectors.cpp -o getEigenVectors

**Command:** ::

	getEigenVectors <options> --vt <text file> --mode <mode index>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| VT matrix file *    	 | File       |``--vtMatrix``      | VT values from ANM script   |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Mode index *           | Integer    |``--mode``          | Specify the index of the    |
|                        |            |                    | mode you wish to target     |
+------------------------+------------+--------------------+-----------------------------+
| Direction              | Boolean    |``--direction``     | Direction of overlap        |
|                        | integer    |                    | correction                  |
|                        | (1 or -1)  |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| Eigenvector file       | Text file containing a      |
|                        | list of extracted           |
|                        | eigenvectors for a specific |
|                        | modes. Vectors are unit     |
|                        | vectors                     |
+------------------------+-----------------------------+

Mean square fluctuation
-------------------------------

Calculates and returns the diagonals of the correlation matrix for a given set of modes.

The user can also compare the msf between two protein complexes. Let's say that the user has performed NMA on two coarse grained models of the same protein, and now wants to compare
if the additional coarse graining decreased the accuracy. If we obtain the same mean square fluctuations for
each residue, then in each model we can say that the results are comparable regardless of the coarse graining
level. Obviously, we must compare only the residues that are common in each model. Hence we specify common residues

**Command:** ::

	meanSquareFluctuations.py <options> --pdb <PDB file> --wMatrix <text file> --vtMatrix <text file> --atomType <string>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *             | File       |``--pdb``           | PDB input file              |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| W matrix file *        | File       |``--wMatrix``	   | W values from ANM script    |
|                        |            |                    | for PDB                     |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``	   | VT values from ANM script   |
|                        |            |                    | for Comparison PDB          |
+------------------------+------------+--------------------+-----------------------------+
| Atom type *            | String     |``--atomType``      | Specify the type of atom to |
|                        |            |                    | be selected in CG models.   |
|                        |            |                    | Only CA or CB accepted.     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Comparison PDB         | File       |``--pdbC``          | When assigned, calculates   |
|                        |            |                    | mean square fluctautions    |
|                        |            |                    | based of common residues    |
|                        |            |                    | between the two proteins    |
+------------------------+------------+--------------------+-----------------------------+
| W matrix file          | File       |``--wMatrixC``	   | When assigned W values from |
| for pdbC               |            |                    | ANM for Comparison PDB      |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file         | File       |``--vtMatrixC``	   | When assigned VT values from|
| for pdbC               |            |                    | ANM for Comparison PDB      |
+------------------------+------------+--------------------+-----------------------------+
| Selected modes         | String     |``--modes``         | MSFs will be calculated     |
|                        |            |                    | over specified modes.       |
|                        |    OR      |                    | Options:                    | 
|	                 |            |                    | 1) Single mode E.g --modes 7|
|                        | Colon      |                    | 2) A range E.g --modes 7:20 |
|                        | Separated  |                    | 3) A list E.g --modes 8,9,11| 
|                        | String     |                    |                             |
|                        |            |                    | If unspecified MSFs will be |   
|                        |    OR      |                    | calculated for the first    |                           
|                        |            |                    | twenty slowest modes (7:27) |
|                        | Comma      |                    |                             | 
|                        | Separated  |                    |                             |
|                        | String     |                    |                             |
+------------------------+------------+--------------------+-----------------------------+ 

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| The following are generated for the PDB and          |
| Comparison PDB (if pdbC was assigned)                |
+------------------------+-----------------------------+
| MSF text file          | MSF for all residues,       |
|                        | calculated over all modes   |
+------------------------+-----------------------------+
| MSF modes text file    | MSF for all residues,       |
|                        | calculated for a specific   |
|                        | mode range                  |
+------------------------+-----------------------------+
| Common residue MSF     | MSF for all common          |
| text file              | residues, calculated over   |
|                        | all modes                   |
+------------------------+-----------------------------+
| Common residue MSF     | MSF for all common          |
| modes text file        | residues, calculated over a |
|                        | specific mode range         |
+------------------------+-----------------------------+

Assembly Covariance
-------------------------------

Calculates and plots Covariance matrices

The user can compare the Covariance between different regions in the biological assembly, or can calcaulate the Covariance across the full assembly complex.
The user also has the option to perform the calculation over a specified list of modes or a mode range. The function also has a zoom option that allows the
user create a Covraiance plot for a particular chain within a particular assymetric unit. 

**Command:** ::

	assemblyCovariance.py <options> --pdb <PDB file> --wMatrix <text file> --vtMatrix <text file> --atomType <string>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *             | File       |``--pdb``           | PDB input file              |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| W matrix file *        | File       |``--wMatrix``	   | W values from ANM script    |
|                        |            |                    | for PDB                     |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``	   | VT values from ANM script   |
|                        |            |                    | for Comparison PDB          |
+------------------------+------------+--------------------+-----------------------------+
| Atom type *            | String     |``--atomType``      | Specify the type of atom to |
|                        |            |                    | be selected in CG models.   |
|                        |            |                    | Only CA or CB accepted.     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Selected modes         | String     |``--modes``         | Covariance will be          |
|                        |            |                    | calculated over specified   |
|                        |    OR      |                    | modes                       |
|                        |            |                    |                             |
|                        | Colon      |                    | Options:                    |
|                        | Separated  |                    | 1) All modes E.g --modes all|
|                        | String     |                    | 2) Single mode E.g --modes 7|
|            	         |            |                    | 3) A range E.g --modes 7:20 |
|                        |    OR      |                    | 4) A list E.g --modes 8,9,11|
|                        |            |                    |                             |
|                        | Comma      |                    | If unspecified Covariance   |
|                        | Separated  |                    | will be  calculated for all |
|                        | String     |                    | modes.                      |
+------------------------+------------+--------------------+-----------------------------+
| Assymetric Units       | String     |``--aUnits``        | Covariance will be          | 
|                        |            |                    | calculated and plotted for  |
|                        |    OR      |                    | specified assyemtric units  |
|                        |            |                    | modes                       | 
|                        | Comma      |                    | Options:                    | 
|                        | Separated  |                    | 1) Single unit              |
|                        | String     |                    |    E.g --aUnits 5           |               
|                        |            |                    | 2) A list of units          |                  
|                        |            |                    |    E.g --aUnits 1,3         | 
|                        |            |                    |                             |
|                        |            |                    | If unspecified Covariance   | 
|                        |            |                    | will be calculated for the  |   
|                        |            |                    | first assymtereic unit in   |
|                        |            |                    | the assembly.               |
+------------------------+------------+--------------------+-----------------------------+
| Zoom                   | Comma      |``--zoom``          | If specified:Covariance will|
|                        | Separated  |                    | be calculated and plotted   |
|                        | String     |                    | for a specified chain in a  |
|                        |            |                    | specified unit.             |
|                        |            |                    | Only format accepts is:     |
|                        |            |                    | [Unit,Chain]                |
|                        |            |                    |    E.g --zoom 1,2           |
|                        |            |                    |        OR                   |
|                        |            |                    |    E.g --zoom 1,B           |
|                        |            |                    | (Chain specifier must match |
|                        |            |                    | chain label in PDB file)    |
|                        |            |                    | The above calculates the    |
|                        |            |                    | covairance for the second   |
|                        |            |                    | chain in the first          |
|                        |            |                    | assymetric unit.            |
+------------------------+------------+--------------------+-----------------------------+
| VMin                   | float      |``--vmin``          | Minimum axes value for plot |
|                        |            |                    | Default: -0.1               |
+------------------------+------------+--------------------+-----------------------------+
| VMax                   | float      |``--vmax``          | Maximum axes value for plot |
|                        |            |                    | Default:  0.1               |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| The following are generated for the PDB and          |
| Comparison PDB (if pdbC was assigned)                |
+------------------------+-----------------------------+
| Covariance Plots       | Covariance Matrices plotted |
|                        | as a Linear Segmented Color | 
|                        | map                         |
+------------------------+-----------------------------+
| Matrix text files      | Covariance Matrices printed |
|                        | in .txt format              | 
|                        |                             |
+------------------------+-----------------------------+

Conformation mode
-------------------------------

Identifies modes responsible for the conformational change of a molecule.

**Command:** ::

	conformationMode.py <options> --pdbConf <PDB file> --pdbANM <PDB file> --vtMatrix <text file> --atomType <string>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Unaligned PDB file *   | File       |``--pdbConf``       | PDB file of the             |
|                        |            |                    | conformational change       |
+------------------------+------------+--------------------+-----------------------------+
| PDB *                  | File       |``--pdbANM``        | PDB file that was useed to  |
|                        |            |                    | run ANM                     |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``      | Eigenavalues obtained from  |
|                        |            |                    | ANM script                  |
+------------------------+------------+--------------------+-----------------------------+
| Atom type *            | String     |``--atomType``      | Specify the type of atom to |
|                        |            |                    | be selected in CG models.   |
|                        |            |                    | Only CA or CB accepted.     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Output file            | File       |``--output``        | Specify a name for the PDB	 |
|                        |            |                    | output file. Default:       |
|                        |            |                    | ModesOfConfChange.pdb       |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| Conformation file      | Text file with the overlap  |
|                        | and correlation of each     |
|                        | mode                        |
+------------------------+-----------------------------+

Combination mode
-------------------------------

Calculates the combined overlap and correlation for specified set of modes to a known conformational change.

**Command:** ::

	combinationMode.py <options> --pdbConf <PDB file> --pdbANM <PDB file> --vtMatrix <text file> --modes <comma separated string> --atomType <string>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Unaligned PDB file *   | File       |``--pdbConf``       | PDB file of the             |
|                        |            |                    | conformational change       |
+------------------------+------------+--------------------+-----------------------------+
| PDB *                  | File       |``--pdbANM``        | PDB file that was useed to  |
|                        |            |                    | run ANM                     |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``      | Eigenavalues obtained from  |
|                        |            |                    | ANM script                  |
+------------------------+------------+--------------------+-----------------------------+
| Modes *                | Integer    |``--modes``         | Calculate the overlap for a |
|                        |            |                    | combination of specific     |
|                        |            |                    | modes. Numbers are          |
|                        |            |                    | separated by commas: 1,5,7  |
+------------------------+------------+--------------------+-----------------------------+
| Atom type *            | String     |``--atomType``      | Specify the type of atom to |
|                        |            |                    | be selected in CG models.   |
|                        |            |                    | Only CA or CB accepted.     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Output file            | File       |``--output``        | Specify a name for the PDB	 |
|                        |            |                    | output file. Default:       |
|                        |            |                    | ModesOfConfChange.pdb       |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| Combination file       | Text file with the overlap  |
|                        | and correlation of each     |
|                        | mode as well as the         |
|                        | combined overlap and        |
|                        | correlation for the modes   |
|                        | specified                   |
+------------------------+-----------------------------+

Mode visualisation
-------------------------------

Generates a trajectory with arrows that can be viewed in the tool VMD

**Command:** ::

	visualiseVector.py <options> --pdb <PDB file> --vectorFile <text file> --mode <int> --atomType <string>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Coarse grained PDB     | File       |``--pdb``           | Coarse grained PDB input    |
| file *                 |            |                    | file                        |
+------------------------+------------+--------------------+-----------------------------+
| Mode index value *     | Ingeter    |``--mode``          | Value specifying the index  |
|                        |            |                    | of the mode                 |
+------------------------+------------+--------------------+-----------------------------+
| Vector file *          | File       |``--vectorFile``    | File containing eigen       |
|                        |            |                    | vectors                     |
+------------------------+------------+--------------------+-----------------------------+
| Atom type *            | String     |``--atomType``      | Specify the type of atom to |
|                        |            |                    | be selected in CG models.   |
|                        |            |                    | Only CA or CB accepted.     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

Outputs are generated in output/VISUALISE directory by default.

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PDB file               | Output PDB to be opened in  |
|                        | VMD                         |
+------------------------+-----------------------------+
| Arrows file            | Tcl script that can be      |
|                        | copied into the VMD TK      |
|                        | console                     |
+------------------------+-----------------------------+

