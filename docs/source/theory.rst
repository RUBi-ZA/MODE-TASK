
Normal Mode Analysis (NMA)
--------------------------

NMA analyses the oscillations of a structure. For proteins, it is useful for studying the large amplitude motions for a selected conformation. The main assumption is that the motions are harmonic. Thus, each normal mode, which is a concerted motion of many atoms, acts as a simple harmonic oscillator, and it is independent of all the other normal modes.

For a harmonic oscillator with a mass :math:`m` supported on a spring with force constant :math:`k`, the potential energy of the system, :math:`V = kx^2`, for an extension :math:`x` leads to the restoring force,

.. math::
	\mathbf{F} = -\frac{dV}{dx} = -kx

By substituting this Hooke’s Law force into Newton’s Law, :math:`\mathbf{F} = m\mathbf{a}` leads to the differential equation,

.. math::
	-kx = m\frac{d^2x}{dt^2}

The solution is,

.. math::
	-kx = -4{\pi}^2v^2mx

with :math:`v` being the frequency of the vibration.

In three dimensions and for a set of :math:`N` atoms, one has the corresponding generalized equation,

.. math::
	-\mathbf{HX} = -4{\pi}^2v^2\mathbf{X}

where :math:`\mathbf{H}` is the :math:`3N`\ x\ :math:`3N` Hessian symmetric matrix of force constants and :math:`\mathbf{X}` is the :math:`3N`\ x\ :math:`1` vector of the positions of the atoms. The solutions of the equation may be cast in the form of an eigenvalue decomposition of :math:`\mathbf{H}` where the eigenvectors are the mass weighted displacements of the normal coordinates; i.e. the independent vibrational motions of the collection of atoms. The corresponding eigenvalues are the negative of the squared normal mode frequencies of each mode. :math:`\mathbf{H}` has exactly six zero eigenvalues for the translations and rotations of the molecule in three dimensional space.

In NMA of proteins, it is central to obtain a good representation of the protein for a proper analysis of the available modes of motion. One approach would be to obtain the Hessian matrix from the second derivatives of the energy for the conformation of interest following a very stringent minimization of the molecule. The latter is important as NMA is based on the harmonic assumption which is only a valid approximation at the bottom of the potential energy minima. For the complex potentials representing interactions in proteins in the water environment, these second derivatives may be obtained through numerical methods which are extremely costly.

Alternatively, :math:`\mathbf{H}` may be approximated by the variance-covariance matrix of the protein, obtained from a molecular dynamics (MD) simulation of suitable length. This is relevant as it may be shown through employing statistical mechanics that the average displacements of nodes from a mean structure for an atom :math:`i, \Delta \mathbf{R}_{i}` are related to those of all other atoms through the relationship, 

.. math::
	< \Delta \mathbf{R}_{i} \cdot \Delta \mathbf{R}_{j} > = 3k_{B}T \mathbf{H}^{-1}

It is important to select the length of the MD simulation for this procedure such that the molecule samples only the state of interest, [1]_ since when more than a single potential energy well along the conformational space of the molecule is sampled, the harmonic assumption would again fail.

As a third alternative for obtaining the Hessian, one may make use of the elastic network property of a folded protein. Here, the assumption is that, once a protein folds to a functionally relevant conformation, its total energy is represented by a simple harmonic potential such that residues are connected to their nearest neighbors via elastic springs. By further employing the assumption that the spring constants, :math:`\gamma`, are equivalent in all pairs of interactions, one arrives at the anisotropic network model (ANM). [2]_

:math:`\mathbf{H}` is thus composed of :math:`N`\ x\ :math:`N` super-elements, i.e., 

.. math::
	\mathbf{H} = 
	 \begin{bmatrix}
	  \mathbf{H}_{1,1} & \mathbf{H}_{1,2} & \cdots & \mathbf{H}_{1,\mathbf{N}} \\
	  \mathbf{H}_{2,1} &         &        & \mathbf{H}_{2,\mathbf{N}} \\
	  \vdots  &         &        & \vdots  \\
	  \mathbf{H}_{\mathbf{N},1} &         &        & \mathbf{H}_{\mathbf{N},\mathbf{N}} 
	 \end{bmatrix}

where each super-element :math:`\mathbf{H_{ij}}`   is a 3x3 matrix that holds the anisotropic information regarding the orientation of nodes :math:`i,j`:	

.. math::
	\mathbf{H_{ij}} = 
	 \begin{bmatrix}
      \partial ^2 V / \partial X_{i} \partial X_{j} & \partial ^2 V / \partial X_{i} \partial Y_{j} & \partial ^2 V / \partial X_{i} \partial Z_{j} \\
      \partial ^2 V / \partial Y_{i} \partial X_{j} & \partial ^2 V / \partial Y_{i} \partial Y_{j} & \partial ^2 V / \partial Y_{i} \partial Z_{j} \\
      \partial ^2 V / \partial Z_{i} \partial X_{j} & \partial ^2 V / \partial Z_{i} \partial Y_{j} & \partial ^2 V / \partial Z_{i} \partial Z_{j} \\
     \end{bmatrix}

Denoting the separation between nodes :math:`i` and :math:`j` by :math:`S_{ij}`, the elements of the off-diagonal super-elements :math:`\mathbf{H_{ij}}` are given by:	

.. math::
	\partial ^2 V / \partial X_{i} \partial Y_{j} = -\gamma (X_{j}-X_{i})(Y_{j}-Y_{i})/S^2_{ij}

and those of the diagonal super-elements :math:`\mathbf{H_{ij}}` are obtained via,

.. math::
	\partial ^2 V / \partial X^2_{i} = \gamma \sum_{j}(X_{j}-X_{i})^2/S^2_{ij} \text{ } \text{ } \text{ for the diagonal terms}

.. math::
	\partial ^2 V / \partial X_{i} \partial Y_{j} = \gamma \sum_{j} (X_{j}-X_{i})(Y_{j}-Y_{i})/S^2_{ij} \text{ } \text{ } \text{ for the off-diagonal terms.}

In this representation of the protein, the structure is coarse-grained at the level of a residue, usually through the coordinates of :math:`\mathbf{C_\alpha}` or :math:`\mathbf{C_\beta}` atoms obtained from the protein data bank. Moreover, pairs of nodes are assumed to interact if they are within a pre-selected cut-off distance, :math:`r_{c}`. For large structures such as viruses, further coarse graining may be shown to well describe the most collective modes of motion. 

The selection of :math:`r_{c}` has been the cause of much research. While it is clear that there is a lower bound where the condition of six zero eigenvalues of :math:`\mathbf{H}` should be satisfied, distances in the range 10-25 Å have been employed in the literature. It has been shown by a systematic study of increasing :math:`r_{c}` that the collective modes of motion do not change beyond a certain value for proteins; i.e. selection of too large :math:`r_{c}` is safer than a too small value. [3]_ The reason for this has been explained by expressing :math:`\mathbf{H}` as the sum of an essential and a trivial part. The essential part of :math:`\mathbf{H}` includes all the local information necessary for the correct representation of the modes. On the other hand, due to the symmetries in a protein originating from the orientational order of closely packed spheres, the effects from the long range neighbors cancel out.

.. rubric:: References

.. [1] C Atilgan, OB Okan, AR Atilgan, "Network-based Models as Tools Hinting at Non-evident Protein Functionality," Annual Review of Biophysics, 41, 205-225 (2012).

.. [2] AR Atilgan, SR Durell, RL Jernigan, MC Demirel, O Keskin, I Bahar, "Anisotropy of Fluctuation Dynamics of Proteins with an Elastic Network Model," Biophysical Journal, 80, 505-515.

.. [3] C Atilgan, OB Okan, AR Atilgan, "Orientational Order Governs Collectivity of Folded Proteins," Proteins: Structure, Function, Bioinformatics, 78, 3363-3375 (2010).


Principle Component Analysis (PCA)
----------------------------------

A molecular dynamics (MD) simulation of a protein provides the positional movements of each atom with
respect to a fixed reference frame at a given time. The mean squared positional fluctuations (variances) of
each atom are readily calculated once the total simulation and sampling times are set. Sufficiency of both
total observation period and the sampling rate are crucial in collecting the data so as to identify biologically
relevant motions.
Let us monitor the variance of each residue’s :math:`\mathbf{C_\alpha}` or :math:`\mathbf{C_\beta}` atom during a MD simulation of a protein. Suppose
that these variances do not change significantly in time, like a stationary process. This suggests that within
the period of observation we have recorded the motion about one (native) conformation. Though constant in
time for a given residue, the variances do change from one residue to another. It is important to distinguish
the differences between the variances of different parts of the protein and to explain the root cause of these
differences; e.g. while loop or disordered regions exhibit high values, relatively rigid parts, such as helices
or sheets display lower variances.

PCA [4]_ operates on the variance-covariance matrix, :math:`\mathbf{C}`, of the protein, obtained from a MD simulation of any
length; thus, the observed process need not be stationary. It is useful in distinguishing the different parts of
the energy landscape sampled during the MD simulation. To obtain :math:`\mathbf{C}`, first the protein coordinates are
superimposed on a reference structure, usually the initial coordinates, or the average coordinates. The
displacement vector for each residue (described by the :math:`\mathbf{C_\alpha}` or :math:`\mathbf{C_\beta}` coordinates of the residue :math:`i`) at a time point
:math:`t, \Delta \mathbf{R}_{i}(t)` is obtained. For a set of :math:`M` recorded coordinates, these are organized in the trajectory fluctuation
matrix of order :math:`3N`\ x\ :math:`M`:

.. math::
	\Delta \mathbf{R} = 
	 \begin{bmatrix}
       \Delta \mathbf{R}_{1}(t_{1}) & \Delta \mathbf{R}_{1}(t_{2}) & \cdot & \Delta \mathbf{R}_{1}(t_{M}) \\
       \Delta \mathbf{R}_{2}(t_{1}) & \Delta \mathbf{R}_{2}(t_{2}) & \cdot & \Delta \mathbf{R}_{2}(t_{M}) \\
       \Delta \mathbf{R}_{3}(t_{1}) & \Delta \mathbf{R}_{3}(t_{2}) & \cdot & \Delta \mathbf{R}_{3}(t_{M}) \\
       \cdot & \cdot & \cdot & \cdot \\
       \cdot & \cdot & \cdot & \cdot \\
       \Delta \mathbf{R}_{n}(t_{1}) & \Delta \mathbf{R}_{n}(t_{2}) &  & \Delta \mathbf{R}_{n}(t_{M}) \\
     \end{bmatrix}

The :math:`3N`\ x\ :math:`3N` :math:`\mathbf{C}` matrix is then obtained via the operation,

.. math::
	\mathbf{C} = \Delta \mathbf{R} \Delta \mathbf{R}^{\mathbf{T}}

If a single energy well along the potential energy surface of a protein is sampled, then :math:`\mathbf{C}` approximates the
inverse Hessian, :math:`\mathbf{H}^{-1}` , as the harmonic approximation applies in this case (see NMA for details). However, if
different parts of the landscape are sampled, the decomposition of :math:`\mathbf{C}` will carry information on all the
regions entered during the simulation. Thus, the diagonalization, 

.. math::
	\mathbf{C = U \Lambda U ^T}

yields the eigenvectors and the corresponding eigenvalues of the :math:`\mathbf{C}` matrix. :math:`\mathbf{\Lambda}` is the :math:`3N`\ x\ :math:`3N` diagonal 
matrix holding the eigenvalues :math:`\lambda_i` with six zero values corresponding to the translations and rotations of the
molecule. The :math:`i\mathrm{^{th}}`  row of the :math:`\mathbf{U}` matrix holding the eigenvector corresponding to the :math:`i\mathrm{^{th}}` eigenvalue. The
trajectory :math:`\Delta \mathbf{R}`  may be projected onto the eigenvectors to obtain the principal components, :math:`q_i`, which are the rows
of the :math:`3N`\ x\ :math:`M` :math:`\mathbf{Q}` matrix.

.. math::
	\mathbf{Q = U} \Delta \mathbf{R}

Since a few principal components usually carry the largest amount of information of the trajectory, the
different regions of the conformational space will manifest as more than one blob in a plot of :math:`q_i` versus :math:`q_j`
where :math:`i` and :math:`j` are small. Furthermore, the size of the blobs in the plots will provide information on the width
of the potential wells sampled. Finally, the time points when passage between different wells occur may be
pinpointed by this method.
The different implementations of the construction of the :math:`\mathbf{C}` matrix and the various ways of decomposing it
have been discussed in detail in the literature, [5]_ and implemented in MODE-TASK.

.. rubric:: References

.. [4] A Amadei, ABM Linssen, HJC Berendsen, “Essential Dynamics of Proteins,” Proteins: Structure, Function and Genetics, 17, 412-425 (1993).

.. [5] CC David, DJ Jacobs, “Principal component analysis: a method for determining the essential dynamics of proteins,” Methods in Molecular Biology, 1084, 193-226 (2014).


