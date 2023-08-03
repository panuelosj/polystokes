# PolyStokes

Source code for the SIGGRAPH 2023 paper [PolyStokes: A Polynomial Model Reduction Method for Viscous Fluid Simulation](https://jpanuelos.com/publications/PolyStokes.html). 

## Prerequisites

* Houdini 18.5+
* CMake 3.9+
* Eigen 3.3+

## Compiling

First make sure the Houdini dev environment is properly setup:
```
cd /opt/hfs.xx
source houdini_setup
```

Return to this project's root directory and build using CMake:
```
mkdir build
cd build
cmake ..
```

From here, either just run `make` or compile in visual studio.

## Usage

Compiling this project should make a new node in Houdini called `HDK Polynomial Stokes Solver`.
This node can be a drop-in replacement for the pressure solver in Houdini's flipsolver node. The default parameters for field names are setup to work as the flipsolver expects.

Parameter usage are as follows:
* *Do Solve:* Perform the Stokes solve. Running this box unchecked will output the same input velocity field, useful if you want just the cell classification data.
* *Export Matrices:* Saves out the main 'Ax=b' matrices and vectors to '.mtx' files.
* *Export Component Matrices:* Saves out component matrices such as the mass matrix and others to '.mtx' files.
* *Export Stats:* Saves out number of degrees of freedom used, runtime timing data.
* *Solver Tolerance:* Tolerance used for the iterative CG solver.
* *Max Solver Iterations:* Iteration cutoff for the CG solver.
* *Active Liquid Boundary:* Number of uniform cells used for the liquid-air surface. Minimum 2 required.
* *Solid Liquid Boundary:* Number of uniform cells used for the solid-liquid surface. Minimum 2 required.
* *Do Reduced Regions:* Run the reduced solve. Unchecking this box will perform a fully uniform Stokes solve.
* *Do Tile:* Split the reduced interior into tiles. Unchecking this box will use only a single large reduced region (not recommended).
* *Reduced Tile Size:* Width of the reduced tiles INCLUDING padding on one side.
* *Reduced Tile Padding:* Number of uniform cells used between reduced tiles.

## Other Content

* `results`: Contains a copy of the paper as well as renders, videos, and figures.
* `scenes`: Example houdini scene files.