## Software documentation file


### 1. Simulation software basics

Presented software is a general-purpose particle simulation toolkit optimized for execution on both GPUs and CPUs.
The code architecture was optimized to perform computations using GPUs, while CPU unit manages the host code.

The host code was written using C++ programming language.

This was achieved by including the OpenCL open-source library (Khronos Group, USA). The project was organized and compiled with the Microsoft Visual Studio Version 16.4.0 (Microsoft Corporation, USA).



>Warning
There are no guaranties that this software will run on your machine.   

### 2. Model file format

Model file contains information regarding basic data structures used by simulation software. Model file is written in plain text format. Each line of the text define one structure. Each line begins with *keyword*, which defines the type of data structure. *Keyword* is followed by set of *parameters* describing a given structure. A **tab** sign is used as a delimiter for *keywords* and *parameters*.

Basic data structures are:
* BEADS / PARTICLES
* SPRINGS
* ANGULAR BONDS


#### BEAD

The most basic model structure - a particle; can be a hard sphere or soft particle depending on force field used.

```text
BEAD  [type] [px] [py] [pz]

[type]         - string value indicating bead type; type must match one of
                 types declared in input file;
[px] [py] [pz] - position of bead along x,y and z axis of simulation box;

```

#### SPRING

Simple linear bond between two beads/particles. Can be permanent or break under tensile force.

```text
SPRING [bead_id_1] [bead_id_2] [rest_length] [stiffness]

[bead_id_1]      - id number of the first bead connected with linear bond;
                   bead id's are assigned automatically starting form 0 and
                   based on their order of appearance inside model file; first
                   bead has id = 0; second id = 1; etc.
[bead_id_2]      - id number of the second bead connected with linear bond;
[rest_length]    - rest length of the linear bond;
[stiffness]      - stiffness of the linear bond;
```

#### ANGLE

```text
  ANGLE [type] [bead_id_1] [bead_id_2] [bead_id_3] [spring_id_1] [spring_id_2] [rest_angle] [stiffness]

[type]         - integer indicating type of bond; not used at the moment
[bead_id_1]    - id of first bead in angular bond
[bead_id_2]    - id of second bead in angular bond
[bead_id_3]    - id of third bead in angular bond
[spring_id_1]  - id of spring bond between bead 1 and 2
[spring_id_2]  - id of spring bond between bead 2 and 3
[rest_angle]   - rest angle of the angular bond;
[stiffness]    - stiffness of the angular bond;
```



### 3. Simulation input file

Simulation input file contains all information required to start calculations. This includes input file name, output data file name, simulation restart file (optional), declarations of bead types and bead interactions, solver and simulation box settings. The order of appearance of data blocks is not important.

Input file commands:

**model**
```text
description:
  Name of the model file. Model file should be placed in the same directory with simulation software.

syntax:
  model    [model_file_path]

arguments:
  [model_file_path] - path to model file

example:
  model   WATER_BOX.txt
```
</br>

**output**			
```text
description:
  Defines the name of the data output file.

syntax:
  output   [output_file_name]

arguments:
  [output_file_name] - output file name

example:
  output   model_output
```
</br>

**output_freq**
```text
description:
  Defines the requency of output data dumps in number of steps.

syntax:
  output_freq   [number_of_steps]

arguments:
  [number_of_steps] - the number of simulation steps followed by a data dump

example:
  output_freq	5000
```
</br>


**export    type**
```text
description:
  Defines which type of molecules the data will be written from.

syntax:
  export    type   [bead_type]

arguments:
  [bead_type] - type of bead structure

example:
  export    type    W
```
</br>


**restart**
```text
description:
  Defines the name of simulation restart file. Restart file should be placed in the same directory with simulation software.

syntax:
  restart   [restart_file_name]

arguments:
  [restart_file_name] - name of the restart file

example:
  restart	model_large_restart.res
```
</br>


**run**
```text
description:
  Defines the number of simulation steps.

syntax:
  run   [simulation_steps]

arguments:
  [simulation_steps] - number of simulation steps to compute

example:
  run   10000
```
</br>

**minimize**
```text
description:
  Number of system energy minimalization steps to compute.

syntax:
  minimize   [number_of_steps]

arguments:
  [number_of_steps] - number of minimalization steps

example:
  minimize  80000
```
</br>


**step**
```text
description:
  Defines simulation time step.

syntax:
  step   [time_step_size]

arguments:
  [time_step_size] - size of simulation time step

example:
  step	0.0001
```
</br>

**sigma**
```text
description:
  Defines the DPD force field random force coefficient.

syntax:
  sigma   [value]

arguments:
  [value] - value of random force coefficient.

example:
  sigma 0.0
```
</br>

**gamma**
```text
description:
  Defines the DPD force field dissipative force coefficient.

syntax:
  gamma   [value]

arguments:
  [value] - dissipative force coefficient

example:
  gamma 0.0
```
</br>

**cutoff**
```text
description:
  Defines DPD force field cutoff radius. This value is used also as system base length.

syntax:
  cutoff   [range]

arguments:
  [range] - length of cutoff radius

example:
  cutoff	7.2
```
</br>

**lambda**
```text
description:
  XXXXX

syntax:
  lambda   [value]

arguments:
  [value] -

example:
  lambda	0.5
```
</br>

**box**
```text
description:
  XXXXX

syntax:
  box  [x] [y] [z]

arguments:
  [x] [y] [z] - size of simulation box along x,y,z axis, expressed in DPD reduced units

example:
  box	  7.2 7.2 7.2
```
</br>

**type**
```text
description:
  Defines the types and masses of beads used in model file.

syntax:
  type   [bead_type]    [bead_mass]

arguments:
  [bead_type]   - type of bead defined in model file
  [bead_mass]   - mass of bead

example:
  type   W  1.0
```
</br>

**pair**			
```text
description:
  Defines the coefficient of conservative force between two types of beads

syntax:
  pair   [bead_type_1] [bead_type_2] [c1]  [c2]

arguments:
  [bead_type_1] - type of bead
  [bead_type_2] - type of bead
  [c1] - conservative force coefficient
  [c2] - force cutoff radius

example:
  pair  W	W	30  7.2
```
</br>


modif_bead	C3	mass	0.0
modif_bead	C3	velocity	0.0	0.0	0.0

Input file example:
```text
model        WATER_BIG_BOX.txt
output       output_file

output_freq  200
minimize     2000
run          20000

export       type	W

sigma        3.000
gamma        4.500
lambda       0.5

cutoff      7.2

box         72.00	72.00	72.00

type        W       	1.0
type        GALUD	2.609
type        GALUN	2.623

pair        W   W       25.4	7.2
pair        W	GALUD	75.0	7.2
pair        W	GALUN	67.0	7.2

pair        GALUN	GALUN   58.50	7.2
pair        GALUD	GALUD	100.00	7.2
pair        GALUD	GALUN	65.00	7.2

step        0.001
```

### 4. Simulation output files

Currently simulation software exports all data into *.XYZ file format. The XYZ file format is a chemical file format. XYZ format specifies the molecule geometry by giving the number of atoms with Cartesian coordinates that will be read on the first line, a comment on the second, and the lines of atomic coordinates in the following lines.[1] The file format is used in computational chemistry programs for importing and exporting geometries.

The general standard of XYZ file used by this software is as follows:
```text
[number_of_beads_in_frame]
FRAME [simulation_step] [time_step] [box_size_x] [box_size_y] [box_size_z]
[bead_type] [X] [Y] [Z]
[bead_type] [X] [Y] [Z]
[bead_type] [X] [Y] [Z]
```

Example of output file contents:

Software allows to save intermediate or final simulation data, required to restart the simulation from the last completed step. The restart file is written as plain text file with *.res extension. The general standard of restart file used by this software is similar to XYZ file format, however it contains additional data fields:

```text
[number_of_beads_in_frame]
FRAME [simulation_step] [time_step] [box_size_x] [box_size_y] [box_size_z]
[bead_type] [X] [Y] [Z] [VX] [VY] [VZ] [FX] [FY] [FZ] [IX] [IY] [IZ]
[bead_type] [X] [Y] [Z] [VX] [VY] [VZ] [FX] [FY] [FZ] [IX] [IY] [IZ]
[bead_type] [X] [Y] [Z] [VX] [VY] [VZ] [FX] [FY] [FZ] [IX] [IY] [IZ]
```
Example of restart file contents:
```text
202988
FRAME	50000	0.0001	152	60	232
	C1	18.2936	34.432	184.938	18.2936	34.432	184.938	-0.0021	 0.005	 0.064	-0.006	-0.019	-0.000
	C1	18.6101	34.616	184.426	18.6101	34.616	184.426	-0.0014	-0.005	-0.060	-0.008	 0.041	 0.033
	C1	18.9245	34.770	183.902	18.9245	34.770	183.902	-0.0071	 0.005	-0.057	 0.026	-0.082	-0.035
	C1	19.2434	34.861	183.366	19.2434	34.861	183.366	-0.0060	-0.018	 0.050	-0.023	 0.088	 0.024
	C1	19.5726	34.914	182.832	19.5726	34.914	182.832	 0.0041	-0.002	-0.067	-0.020	-0.031	 0.047
	C1	19.9048	34.891	182.297	19.9048	34.891	182.297	 0.0002	-0.015	-0.067	 0.082	-0.032	-0.119
  ....
```
