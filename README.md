# Structural Bioinformatics Pipeline

This pipeline was developed for a university project aimed at simulating the pore domain of KCNQ1, a channel protein, embedded in a POPC membrane.

The system was prepared using CHARMM-GUI, after selecting with VMD the amino acids that constitute the pore domain and uploading the selected structure to the CHARMM-GUI web interface.

At the beginning of the pipeline, some switches are defined and initially set to `0`. When all switches are set to `0`, the pipeline does not run any block. By setting a
witch to `1`, the corresponding block of the pipeline is activated.
The pipeline is divided into three main blocks:

## 1. Simulations

This block handles the simulations, starting from equilibration and continuing with the NPT simulation, performed at fixed pressure, which is needed to determine the box volume to be used 

in the following step: the NVT simulation, performed at fixed volume.

The simulations were run using GROMACS. The SLURM scripts used for this step can be found in the `scripts_cluster` directory and are named:

- `equilibration.slurm`
- `npt_3ns.slurm`
- `nvt_3ns.slurm`

This block can be interrupted and restarted from three different checkpoints:

- start of equilibration;
- end of equilibration / start of NPT;
- end of NPT / start of NVT.


## 2. Energy estimation

This block extracts, from the NPT and NVT `.edr` files, the values of total energy, potential energy, temperature, and pressure over time for both simulations.

The `.edr` files store the parameters and energetic information from the previous simulations. This step helps to understand whether GROMACS handled the given 

parameters correctly before starting the full analysis.

In this part, the analysis is performed by the Python script `script1.py`, located in the `1a` directory.

## 3. Structural analysis

This block performs the analysis of four main topics:

- COM, center of mass
- Volume, calculated using the convex hull
- RMSD
- Lipid atoms located within 5 Å of the channel pore

The script used by this part of the pipeline can be found in the `1b` directory and is named `full_1b.py`.

This section adds more information to help understand whether the observed results are coherent or whether problems occurred during the simulation process.

Possible problems can be divided into two main types:

- technical problems, for example something went wrong in the commands or scripts;
- biological or interpretation-related problems, for example the analysis may not be correctly capturing the feature of interest.

## Note on the convex hull

The convex hull is a function available in Python through packages such as MDAnalysis. Its definition is:

> The convex hull is the tightest convex boundary enclosing all points in a given set.

In this case, the given set of points is represented by the coordinates of the outer atoms of the pore.

Since we are analyzing a channel protein, part of its internal volume is empty. However, the convex hull still represents a reasonable approximation of the protein volume.

I chose this method instead of the radius of gyration because the radius of gyration approximates the structure as a sphere. This sphere can have a diameter related to the 
width or height of the protein in a selected region, leading to a much rougher approximation of the protein volume.
