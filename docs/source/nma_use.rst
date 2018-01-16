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
| Coarse grain level     | Comma      |``--cg``            | Level/Levels by which to    |
|                        | Separted   |                    | coarse grain protein.       |
|                        | String     |                    | Lower is less coarse grained|
|                        |            |                    |                             |
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

Constructs an elastic network model of a protein complex and solves for the eigenvalues and eigenvectors of the system. 

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

Mean square fluctuation
-------------------------------

Calculates and returns the diagonals of the correlation matrix for a given set of modes.

The user can also compare the msf between two protein complexes. Let's say that the user has performed NMA on two coarse grained models of the same protein, and now wants to compare
if the additional coarse graining decreased the accuracy. If we obtain the same mean square fluctuations for
each residue, then in each model we can say that the results are comparable regardless of the coarse graining
level. Obviously, we must compare only the residues that are common in each model. Hence, we specify common residues.

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
|                        |            |                    | mean square fluctuations    |
|                        |            |                    | based on common residues    |
|                        |            |                    | between the two systems.    |
+------------------------+------------+--------------------+-----------------------------+
| W matrix file          | File       |``--wMatrixC``	   | When assigned, calculates   |
| for pdbC               |            |                    | W values from ANM for       |
|                        |            |                    | Comparison PDB              |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file         | File       |``--vtMatrixC``	   | When assigned, calculates   |
| for pdbC               |            |                    | VT values from ANM for      |
|                        |            |                    | Comparison PDB              |
+------------------------+------------+--------------------+-----------------------------+
| Selected modes         | String     |``--modes``         | MSFs will be calculated     |
|                        |            |                    | over specified modes.       |
|                        |    OR      |                    | Options:                    | 
|	                 |            |                    | 1) Single mode E.g          |
|	                 | Colon      |                    | --modes 7;                  |
|                        | Separated  |                    | 2) A range E.g --modes 7:20;|
|                        | String     |                    | 3) A list E.g --modes 8,9,11| 
|                        |            |                    |                             |
|                        |    OR      |                    | If unspecified MSFs will be |   
|                        |            |                    | calculated for the first    |                           
|                        | Comma      |                    | twenty slowest modes (7:27) |
|                        | Separated  |                    |                             | 
|                        | String     |                    |                             |
|                        |            |                    |                             |
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

The user can compare the Covariance between different regions in the biological assembly, or can calculate the Covariance across the full assembly complex.
The user also has the option to perform the calculation over a specified list of modes or a mode range. The function also has a zoom option that allows the
user to create a Covariance plot for a particular chain within a particular asymmetric unit. 

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
|                        | Separated  |                    | 1) All modes E.g            |
|                        | String     |                    | --modes all;                |
|                        |            |                    | 2) Single mode E.g          |
|                        |    OR      |                    | --modes 7;                  |
|            	         |            |                    | 3) A range E.g --modes 7:20;|
|                        | Comma      |                    | 4) A list E.g --modes 8,9,11|
|                        | Separated  |                    |                             |
|                        | String     |                    | If unspecified, Covariance  |
|                        |            |                    | will be  calculated for all |
|                        |            |                    | modes.                      |
+------------------------+------------+--------------------+-----------------------------+
| Asymmetric Units       | String     |``--aUnits``        | Covariance will be          | 
|                        |            |                    | calculated and plotted for  |
|                        |    OR      |                    | specified asymmetric units  |
|                        |            |                    |                             | 
|                        | Comma      |                    | Options:                    | 
|                        | Separated  |                    | 1) Single unit              |
|                        | String     |          	   |    E.g --aUnits 5;          |               
|                        |            |                    | 2) A list of units          |                  
|                        |            |                    |    E.g --aUnits 1,3         | 
|                        |            |                    |                             |
|                        |            |                    | If unspecified, Covariance  | 
|                        |            |                    | will be calculated for the  |   
|                        |            |                    | first asymmetric unit in    |                         
|                        |            |                    | the assembly.               |                           
+------------------------+------------+--------------------+-----------------------------+ 
| Zoom                   | Comma      |``--zoom``          | If specified, Covariance    | 
|                        | Separated  |                    | will be calculated and      |
|                        | String     |                    | plotted for a specified     |
|                        |            |                    | chain in a specified unit.  | 
|                        |            |                    | Only format accepts is:     | 
|                        |            |                    | [Unit,Chain]                |
|                        |            |          	   |    E.g --zoom 1,2           |               
|                        |            |                    |        OR                   |                  
|                        |            |                    |    E.g --zoom 1,B           | 
|                        |            |                    |(Chain specifier must match  |
|                        |            |                    | chain label in PDB file)    |
|                        |            |                    | The above calculates the    | 
|                        |            |                    | covariance for the second   |   
|                        |            |                    | chain in the first          |                         
|                        |            |                    | asymmetric unit.            |                           
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
| PDB *                  | File       |``--pdbANM``        | PDB file that was used to   |
|                        |            |                    | run ANM                     |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``      | Eigenvalues obtained from   |
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

Calculates the combined overlap and correlation for specified set of modes to a known conformational change. This script also calculates the overlap and correlation per chain in each
asymmetric unit for the specified modes. This allows the user to determine which parts of the complex, in each mode, contribute the most to the overall conformational change.

**Command:** ::

	combinationMode.py <options> --pdbConf <PDB file> --pdbANM <PDB file> --vtMatrix <text file> --modes <comma separated string> --atomType <string>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Unaligned PDB file *   | File       |``--pdbConf``       | PDB file of the             |
|                        |            |                    | conformational change       |
+------------------------+------------+--------------------+-----------------------------+
| PDB *                  | File       |``--pdbANM``        | PDB file that was used to   |
|                        |            |                    | run ANM                     |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``      | Eigenvalues obtained from   |
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
| Break down per unit    | Text file with the overlap  |
| file                   | and correlation calculated  |
|                        | for each chain in each      |
|                        | asymmetric unit in the      |
|                        | complex. Calculations are   |
|                        | performed for each specified|
|                        | mode.                       |
+------------------------+-----------------------------+



Mode visualisation
-------------------------------

Generates a set of frames, where eigenvectors are plotted as a set of unit vectors multiplied by an increasing factor in each frame. Vectors are also plotted as arrows that can be viewed in the tool VMD

**Command:** ::

	visualiseVector.py <options> --pdb <PDB file> --vtMatrix <text file> --mode <int> --atomType <string> --direction <int>

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
| VT matrix file *    	 | File       |``--vtMatrix``      | VT values from ANM script   |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Atom type *            | String     |``--atomType``      | Specify the type of atom to |
|                        |            |                    | be selected in CG models.   |
|                        |            |                    | Only CA or CB accepted.     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Direction              | Boolean    |``--direction``     | Direction of overlap        |
|                        | integer    |                    | correction. Default = 1     |
|                        | (1 or -1)  |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Arrow head             | float      |``--head``          | Radius of cone that forms   |
|                        |            |                    | the head of each vector     |
|                        |            |                    | arrow                       |
+------------------------+------------+--------------------+-----------------------------+
| Arrow tail             | float      |``--tail``          | Radius of cylinder that     |
|                        |            |                    | forms the tail of each      |
|                        |            |                    | vector arrow                |
+------------------------+------------+--------------------+-----------------------------+
| Arrow length           | float      |``--arrowLength``   | Specify a factor by which   |
|                        |            |                    | to increase or decrease     |
|                        |            |                    | the length of each arrow    |
|                        |            |                    | E.g                         |
|                        |            |                    |  --arrowLength 2            |
|                        |            |                    |  doubles the default length |
|                        |            |                    |  and                        |
|                        |            |                    |  --arrowLength 0.5          |
|                        |            |                    |  halves the default length  |
+------------------------+------------+--------------------+-----------------------------+
| Colours                | Comma      |``--colourByChain`` | Colour the vectors arrows   |
|                        | Separted   |                    | of each chain.              |
|                        | String     |                    | E.g  for a two chain protein|
|                        |            |                    |   --colourByChain blue,red  | 
|                        |            |                    | will colour the arrows of   |
|                        |            |                    | Chain A as blue and         |
|                        |            |                    | Chain B as red              |
+------------------------+------------+--------------------+-----------------------------+
| Asymmetric Units       | String     |``--aUnits``        | Vector frames and arrows    | 
|                        |            |                    | will be plotted for         |
|                        |    OR      |                    | specified asymmetric units  |
|                        |            |                    |                             | 
|                        | Comma      |                    | Options:                    | 
|                        | Separated  |                    | 1) Single unit              |
|                        | String     |          	   |    E.g --aUnits 5           |               
|                        |            |                    | 2) A list of units          |                  
|                        |            |                    |    E.g --aUnits 1,3         | 
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Chain                  | String     |``--chain``         | Draws arrows only for the   |
|                        |            |                    | specified chain.            |
|                        |            |                    | This option only accepts    |
|                        |            |                    | a single chain              |                           
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

