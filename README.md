Semi-Monocoque Theory
===============

This project tries to reproduce the computations used in Semi-Monocoque theory.

The following elements are computed:

1. **Section Properties**:
	- Shear center
	- Torsional moment of inertia
	- Symmetry properties
2. ** Loads on elements**:
	- Axial loads on stringers
	- Shear fluxes on panels 
3. **Elements of corrective solutions**:
	- **L** and **H** matrices
	- $\beta^2$ and $\lambda$

### The project is composed of:

1. **Section.py** : a python class that performs all the computations
2. Some **Ipython notebooks** that show some examples of applications

### Note:

Most of computations are performed my means of the Python symbolic tool **Sympy** and it is **not** optimized for computation of *big* and/or *non symmetric* sections. 