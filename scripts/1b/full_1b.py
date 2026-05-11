#!/usr/bin/env python3

import argparse
import os

import MDAnalysis as mda
from MDAnalysis.analysis import rms

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial import ConvexHull


def analyze_com(u, protein, out_dir, sim_name):
    """
    Analize COM position of the protein over time.
    Save:
    - CSV with COM coordinates and displacement relative to the initial frame
    - COM x/y/z plot
    - COM displacement plot
    """

    time = []
    com_x = []
    com_y = []
    com_z = []
    displacement = []

    com_initial = None

    for ts in u.trajectory:
        """
        The for loop iterates over the trajectory frames,
        calculating the COM of the protein and its displacement
        from the initial position, then appending the results to
        the respective lists.
        """
        time_ps = ts.time
        com = protein.center_of_mass()
        """the comand protein.center_of_mass() cames from MDAnalysis
        and it's part of the set of commands taht can be applied to 
        AtomGroup objects, in this case the 'protein' selection."""

        if com_initial is None:
            com_initial = com.copy()

        com_displacement = np.linalg.norm(com - com_initial)
        """The displacement is calculated as the Euclidean distance
          between the current COM and the initial COM position,
          the np.linalg.norm function computes this distance by taking 
          the square root of the sum of the squared differences 
          between the current and initial COM coordinates."""

        time.append(time_ps)
        com_x.append(com[0])
        com_y.append(com[1])
        com_z.append(com[2])
        displacement.append(com_displacement)

    mean_displacement = np.mean(displacement)
    std_displacement = np.std(displacement, ddof=1)

    df = pd.DataFrame({
        "Time_ps": time,
        "COM_x_A": com_x,
        "COM_y_A": com_y,
        "COM_z_A": com_z,
        "COM_displacement_A": displacement,
        "Mean_COM_displacement_A": [mean_displacement] * len(time),
        "Std_COM_displacement_A": [std_displacement] * len(time)
    })

    csv_path = os.path.join(out_dir, f"{sim_name}_com.csv")
    df.round(4).to_csv(csv_path, index=False)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    ax[0].plot(time, com_x, label="COM X")
    ax[0].plot(time, com_y, label="COM Y")
    ax[0].plot(time, com_z, label="COM Z")
    ax[0].set_xlabel("Time (ps)")
    ax[0].set_ylabel("COM coordinate (Å)")
    ax[0].set_title(f"{sim_name.upper()} - Protein COM coordinates")
    ax[0].legend()

    ax[1].plot(time, displacement, label="COM displacement", alpha=0.7)
    ax[1].plot(time, [mean_displacement] * len(time), label="Mean displacement", linestyle="--")
    ax[1].set_xlabel("Time (ps)")
    ax[1].set_ylabel("COM displacement (Å)")
    ax[1].set_title(f"{sim_name.upper()} - COM displacement")
    ax[1].legend()

    plt.tight_layout()

    fig_path = os.path.join(out_dir, f"{sim_name}_com.png")
    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

    return df


def analyze_volume(u, protein, out_dir, sim_name):
    """
    Computes the approximate volume of the protein using Convex Hull.
    Saves:
    - CSV with volume and delta volume
    - plot of volume and delta volume
    """

    time = []
    volume = []

    for ts in u.trajectory:
        """
        The for cycle iterates over each frame of the trajectory.
        At every frame, the protein atomic coordinates are updated,
        so the volume is calculated on the current protein conformation.
        """
        time_ps = ts.time
        coords = protein.positions
        """
        protein.positions returns the Cartesian coordinates of all atoms
        belonging to the protein AtomGroup in the current trajectory frame.
        These coordinates are used as the spatial points from which the
        protein volume is estimated.
        """

        hull = ConvexHull(coords)
        """
        ConvexHull is a scipy function that builds the smallest convex
        polyhedron able to contain all the selected atomic coordinates.
        In this script, it is used as an approximation of the protein volume,
        because the hull encloses the outermost atoms of the protein structure.
        """

        time.append(time_ps)
        volume.append(hull.volume)
        """
        hull.volume returns the volume enclosed by the convex hull.
        This value is appended frame by frame, producing a time series
        that describes how the approximate protein volume changes during
        the simulation.
        """

    delta_volume = [0]

    for i in range(1, len(volume)):
        delta = volume[i] - volume[i - 1]
        """
        Delta volume is calculated as the difference between the current
        volume and the volume of the previous frame. This gives the
        frame-to-frame variation of the protein approximate volume.
        """
        delta_volume.append(delta)

    mean_volume = np.mean(volume)
    std_volume = np.std(volume, ddof=1)

    mean_delta_volume = np.mean(delta_volume)
    std_delta_volume = np.std(delta_volume, ddof=1)

    df = pd.DataFrame({
        "Time_ps": time,
        "Volume_A3": volume,
        "Mean_volume_A3": [mean_volume] * len(time),
        "Std_volume_A3": [std_volume] * len(time),
        "Delta_volume_A3": delta_volume,
        "Mean_delta_volume_A3": [mean_delta_volume] * len(time),
        "Std_delta_volume_A3": [std_delta_volume] * len(time)
    })

    csv_path = os.path.join(out_dir, f"{sim_name}_volume.csv")
    df.round(4).to_csv(csv_path, index=False)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    ax[0].plot(time, volume, label="Convex Hull volume", alpha=0.7)
    ax[0].plot(time, [mean_volume] * len(time), label="Mean volume", linestyle="--")
    ax[0].set_xlabel("Time (ps)")
    ax[0].set_ylabel("Convex Hull volume (Å³)")
    ax[0].set_title(f"{sim_name.upper()} - Protein volume")
    ax[0].legend()

    ax[1].plot(time, delta_volume, label="Delta volume", alpha=0.7)
    ax[1].plot(time, [mean_delta_volume] * len(time), label="Mean delta volume", linestyle="--")
    ax[1].set_xlabel("Time (ps)")
    ax[1].set_ylabel("Delta volume (Å³)")
    ax[1].set_title(f"{sim_name.upper()} - Delta volume")
    ax[1].legend()

    plt.tight_layout()

    fig_path = os.path.join(out_dir, f"{sim_name}_volume.png")
    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

    return df


def analyze_rmsd(tpr_file, xtc_file, out_dir, sim_name):
    """
    Calculates RMSD of the protein's backbone atoms with respect to the initial frame.
    Saves:
    - CSV with RMSD
    - plot of RMSD
    """

    u = mda.Universe(tpr_file, xtc_file)
    """
    mda.Universe loads the topology and trajectory files into an MDAnalysis
    Universe object. The Universe contains both the structural information
    and the time-dependent atomic coordinates of the simulation.
    """
    ref = mda.Universe(tpr_file, xtc_file)
    """
    A second Universe is created to be used as reference structure.
    In this case, the RMSD will be calculated by comparing each trajectory
    frame against the initial frame of this reference Universe.
    """

    backbone = u.select_atoms("protein and backbone")
    """
    select_atoms is an MDAnalysis command used to extract a specific group
    of atoms from the Universe. Here, the selection keeps only the backbone
    atoms of the protein, which are commonly used for RMSD analysis because
    they represent the global fold of the protein.
    """

    if backbone.n_atoms == 0:
        raise ValueError("La selezione 'protein and backbone' non contiene atomi.")

    R = rms.RMSD(
        u,
        ref,
        select="protein and backbone",
        ref_frame=0
    )
    """
    rms.RMSD is an MDAnalysis analysis class that calculates the Root Mean
    Square Deviation of a selected group of atoms over time. RMSD measures
    how much the structure deviates from a reference conformation, after
    optimal structural alignment. Here, each frame is compared to frame 0.
    """

    R.run()
    """
    The run() method executes the RMSD analysis over the trajectory.
    It iterates through the frames, aligns the selected atoms to the
    reference structure, calculates the RMSD value, and stores the results.
    """

    rmsd_data = R.results.rmsd
    """
    R.results.rmsd contains the numerical output of the RMSD analysis.
    The resulting array stores, for each frame, the frame index, the
    simulation time, and the RMSD value of the selected atoms.
    """
    frames = rmsd_data[:, 0]
    time = rmsd_data[:, 1]
    rmsd_values = rmsd_data[:, 2]
    """
    The RMSD result array is divided into separate columns.
    The first column contains frame numbers, the second contains simulation
    time, and the third contains the RMSD values calculated for the
    protein backbone.
    """
    mean_rmsd = np.mean(rmsd_values)
    std_rmsd = np.std(rmsd_values, ddof=1)

    df = pd.DataFrame({
        "Frame": frames,
        "Time_ps": time,
        "Backbone_RMSD_A": rmsd_values,
        "Mean_Backbone_RMSD_A": [mean_rmsd] * len(time),
        "Std_Backbone_RMSD_A": [std_rmsd] * len(time)
    })

    csv_path = os.path.join(out_dir, f"{sim_name}_rmsd.csv")
    df.round(4).to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(time, rmsd_values, label="Backbone RMSD")
    ax.plot(time, [mean_rmsd] * len(time), linestyle="--", label="Mean RMSD")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Backbone RMSD (Å)")
    ax.set_title(f"{sim_name.upper()} - Backbone RMSD")
    ax.legend()

    plt.tight_layout()

    fig_path = os.path.join(out_dir, f"{sim_name}_rmsd.png")
    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

    return df


def analyze_near_lipids(u, out_dir, sim_name, lipid_resname="POPC", cutoff=5.0):
    """
    counts the number of lipids near the protein over time.
    By default, it selects POPC residues within 5 Å of the protein.
    Saves:
    - CSV with the number of nearby lipids
    - plot of nearby lipids
    """

    lipids = u.select_atoms(f"resname {lipid_resname}")
    """
    This command selects all atoms belonging to residues with the specified
    lipid residue name. In MDAnalysis, resname is used to identify residue
    types, so this selection extracts the lipid population that will be
    monitored around the protein.
    """
    if lipids.n_atoms == 0:
        raise ValueError(f"Nessun lipide trovato con resname {lipid_resname}.")

    time = []
    near_lipids = []

    selection = f"resname {lipid_resname} and around {cutoff} protein"
    """
    The selection string combines the lipid residue name with the MDAnalysis
    'around' selection keyword. The 'around' command selects atoms located
    within a given distance from another selection, in this case the protein.
    Therefore, this expression identifies lipid atoms close to the protein.
    """
    for ts in u.trajectory:
        """
        The for cycle iterates over each frame of the trajectory.
        Since atomic positions change at every frame, the nearby lipid selection
        is recalculated through time to follow changes in protein-lipid contacts.
        """
        time_ps = ts.time

        near_lipid_atoms = u.select_atoms(selection)
        """
        The selection is evaluated again at each trajectory frame.
        This is important because the atoms within the cutoff distance
        can change over time as both the protein and lipids move during
        the molecular dynamics simulation.
        """
        near_lipid_residues = near_lipid_atoms.residues
        """
        The residues attribute converts the selected lipid atoms into their
        corresponding residues. This avoids counting individual atoms and
        instead allows the script to count whole lipid molecules near the
        protein.
        """ 
        n_near_lipids = near_lipid_residues.n_residues
        """
        n_residues returns the number of residues in the selected group.
        Since each lipid molecule is represented as a residue, this value
        corresponds to the number of lipid molecules found near the protein
        in the current frame.
        """
        time.append(time_ps)
        near_lipids.append(n_near_lipids)
        
    mean_near_lipids = np.mean(near_lipids)
    std_near_lipids = np.std(near_lipids, ddof=1)

    df = pd.DataFrame({
        "Time_ps": time,
        "Near_lipids": near_lipids,
        "Mean_near_lipids": [mean_near_lipids] * len(time),
        "Std_near_lipids": [std_near_lipids] * len(time)
    })

    csv_path = os.path.join(out_dir, f"{sim_name}_near_lipids.csv")
    df.round(4).to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(time, near_lipids, label=f"{lipid_resname} near protein", alpha=0.7)
    ax.plot(time, [mean_near_lipids] * len(time), linestyle="--", label="Mean near lipids")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(f"Number of {lipid_resname} residues near protein")
    ax.set_title(f"{sim_name.upper()} - Lipids within {cutoff} Å from protein")
    ax.legend()

    plt.tight_layout()

    fig_path = os.path.join(out_dir, f"{sim_name}_near_lipids.png")
    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Analysis script for MD trajectory: COM, volume, RMSD and near lipids."
    )

    parser.add_argument("tpr_file", help="Input .tpr file")
    parser.add_argument("xtc_file", help="Input .xtc file")
    parser.add_argument("out_dir", help="Output directory")
    parser.add_argument("sim_name", help="Simulation name, e.g. npt or nvt")

    parser.add_argument(
        "--lipid-resname",
        default="POPC",
        help="Lipid residue name used for near-lipid analysis. Default: POPC"
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="Cutoff distance in Å for near-lipid analysis. Default: 5.0"
    )

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("======================================")
    print(f"Running analysis for: {args.sim_name}")
    print(f"TPR file: {args.tpr_file}")
    print(f"XTC file: {args.xtc_file}")
    print(f"Output directory: {args.out_dir}")
    print("======================================")

    if not os.path.isfile(args.tpr_file):
        raise FileNotFoundError(f"TPR file not found: {args.tpr_file}")

    if not os.path.isfile(args.xtc_file):
        raise FileNotFoundError(f"XTC file not found: {args.xtc_file}")

    u = mda.Universe(args.tpr_file, args.xtc_file)

    protein = u.select_atoms("protein")

    if protein.n_atoms == 0:
        raise ValueError("La selezione 'protein' non contiene atomi.")

    print(f"Total atoms: {u.atoms.n_atoms}")
    print(f"Protein atoms: {protein.n_atoms}")
    print(f"Protein residues: {protein.residues.n_residues}")
    print("")

    print("Analyzing protein COM...")
    analyze_com(u, protein, args.out_dir, args.sim_name)

    print("Analyzing protein volume...")
    analyze_volume(u, protein, args.out_dir, args.sim_name)

    print("Analyzing backbone RMSD...")
    analyze_rmsd(args.tpr_file, args.xtc_file, args.out_dir, args.sim_name)

    print("Analyzing near lipids...")
    analyze_near_lipids(
        u,
        args.out_dir,
        args.sim_name,
        lipid_resname=args.lipid_resname,
        cutoff=args.cutoff
    )

    print("")
    print("Analysis completed.")
    print(f"Results saved in: {args.out_dir}")


if __name__ == "__main__":
    main()