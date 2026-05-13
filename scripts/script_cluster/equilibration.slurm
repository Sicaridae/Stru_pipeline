#!/bin/bash
#SBATCH --job-name=kcnq1_eq
###esample path for the output and error files, to be changed with the real path of the user
#SBATCH --output=/home/studmol/yy2025/jfedeli/data/slurm/kcnq1_eq.out
#SBATCH --error=/home/studmol/yy2025/jfedeli/data/slurm/kcnq1_eq.err  
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

module load molecular-dynamics
module load gromacs/2023.5
set -e

BASE_DIR="/home/pesce/bio_strutturale/pipe_line/chermm_steps"
GROMACS_DIR="$BASE_DIR/gromacs"

cd "$BASE_DIR"
### gromppp prepare the file .tpr (binary input file for mdrun)
###parameters:
# -f: mdp file with the parameters of the simulation
# -c: structure file to start the simulation
# -r: reference structure for position restraints
# -p: topology file
# -n: index file
# -o: output file name for the .tpr file
gmx grompp -f step6.0_minimization.mdp -c step5_input.gro -r step5_input.gro -p "$GROMACS_DIR/topol.top" -n "$GROMACS_DIR/index.ndx" -o step6.0_minimization.tpr
### mdrun run the simulation
### parameters:
# -deffnm: base name for output files
# -ntomp: number of threads to use for the simulation
gmx mdrun -deffnm step6.0_minimization -ntomp 8

echo "Equilibration step 0 completed. Output files: step6.0_minimization.gro, step6.0_minimization.cpt"

gmx grompp -f step6.1_equilibration.mdp -c step6.0_minimization.gro -r step5_input.gro -p "$GROMACS_DIR/topol.top" -n "$GROMACS_DIR/index.ndx" -o step6.1_equilibration.tpr
gmx mdrun -deffnm step6.1_equilibration -ntomp 8

echo "Equilibration step 1 completed. Output files: step6.1_equilibration.gro, step6.1_equilibration.cpt"

gmx grompp -f step6.2_equilibration.mdp -c step6.1_equilibration.gro -r step5_input.gro -p "$GROMACS_DIR/topol.top" -n "$GROMACS_DIR/index.ndx" -o step6.2_equilibration.tpr
gmx mdrun -deffnm step6.2_equilibration -ntomp 8

echo "Equilibration step 2 completed. Output files: step6.2_equilibration.gro, step6.2_equilibration.cpt"

gmx grompp -f step6.3_equilibration.mdp -c step6.2_equilibration.gro -r step5_input.gro -p "$GROMACS_DIR/topol.top" -n "$GROMACS_DIR/index.ndx" -o step6.3_equilibration.tpr
gmx mdrun -deffnm step6.3_equilibration -ntomp 8

echo "Equilibration step 3 completed. Output files: step6.3_equilibration.gro, step6.3_equilibration.cpt"

gmx grompp -f step6.4_equilibration.mdp -c step6.3_equilibration.gro -r step5_input.gro -p "$GROMACS_DIR/topol.top" -n "$GROMACS_DIR/index.ndx" -o step6.4_equilibration.tpr
gmx mdrun -deffnm step6.4_equilibration -ntomp 8

echo "Equilibration step 4 completed. Output files: step6.4_equilibration.gro, step6.4_equilibration.cpt"

gmx grompp -f step6.5_equilibration.mdp -c step6.4_equilibration.gro -r step5_input.gro -p "$GROMACS_DIR/topol.top" -n "$GROMACS_DIR/index.ndx" -o step6.5_equilibration.tpr
gmx mdrun -deffnm step6.5_equilibration -ntomp 8

echo "Equilibration step 5 completed. Output files: step6.5_equilibration.gro, step6.5_equilibration.cpt"

gmx grompp -f step6.6_equilibration.mdp -c step6.5_equilibration.gro -r step5_input.gro -p "$GROMACS_DIR/topol.top" -n "$GROMACS_DIR/index.ndx" -o step6.6_equilibration.tpr
gmx mdrun -deffnm step6.6_equilibration -ntomp 8

echo "Equilibration step 6 completed. Output files: step6.6_equilibration.gro, step6.6_equilibration.cpt"

### often, the last step of the equilibration, keep some leves of restrains

NPT_INPUT_DIR="/home/pesce/bio_strutturale/pipe_line/output/buffer/npt_input"
NVT_INPUT_DIR="/home/pesce/bio_strutturale/pipe_line/output/buffer/nvt_input"

mkdir -p "$NPT_INPUT_DIR" "$NVT_INPUT_DIR"
### check of pressure and temperature on npt_3ns.mdp and nvt_3ns.mdp, if they are not set to
### 310 K and 1 atm, set them to these values, and make a backup of the original files before 
### modifying them, in case the user want to check the original files or restore them"
cp "$GROMACS_DIR/npt_3ns.mdp" "$GROMACS_DIR/npt_3ns.mdp.bak"
cp "$GROMACS_DIR/nvt_3ns.mdp" "$GROMACS_DIR/nvt_3ns.mdp.bak"
echo "checking and setting ref_t and ref_p in npt_3ns.mdp"
grep "^ref_t" "$GROMACS_DIR/npt_3ns.mdp" || true ### the || true is to avoid the error message in case the parameter is not found, but in this case i 
grep "^ref_p" "$GROMACS_DIR/npt_3ns.mdp" || true ### t should be there because it's a copy of the original file, but just in case, if the parameter is not found, it will be added with the correct value

echo "Setting ref_t and ref_p in npt_3ns.mdp to 310 K and 1 atm"

sed -i -E 's/^ref_t\s*=.*/ref_t                   = 310 310 310/' "$GROMACS_DIR/npt_3ns.mdp"
sed -i -E 's/^ref_p\s*=.*/ref_p                   = 1.0     1.0/' "$GROMACS_DIR/npt_3ns.mdp"

echo "Backup file: $GROMACS_DIR/npt_3ns.mdp.bak"

grep "^ref_t" "$GROMACS_DIR/nvt_3ns.mdp" || true
grep "^pcoupl" "$GROMACS_DIR/nvt_3ns.mdp" || true

sed -i -E 's/^ref_t\s*=.*/ref_t                   = 310 310 310/' "$GROMACS_DIR/nvt_3ns.mdp"
sed -i -E 's/^pcoupl\s*=.*/pcoupl                  = no/' "$GROMACS_DIR/nvt_3ns.mdp"

echo "Backup file: $GROMACS_DIR/nvt_3ns.mdp.bak"
### move all the files needed for the NPT and NVT simulations to the respective input directories
cp step6.6_equilibration.gro "$NPT_INPUT_DIR/"
cp step6.6_equilibration.cpt "$NPT_INPUT_DIR/"
cp "$GROMACS_DIR/topol.top" "$NPT_INPUT_DIR/"
cp "$GROMACS_DIR/index.ndx" "$NPT_INPUT_DIR/"
cp "$GROMACS_DIR/npt_3ns.mdp" "$NPT_INPUT_DIR/"

cp "$GROMACS_DIR/nvt_3ns.mdp" "$NVT_INPUT_DIR/"
### check of the existence of the essential files for the NPT simulation, if they are missing, exit with an error message
if [ ! -s "$NPT_INPUT_DIR/step6.6_equilibration.gro" ] || [ ! -s "$NPT_INPUT_DIR/step6.6_equilibration.cpt" ]; then
    echo "ERROR: missing essential NPT input files in $NPT_INPUT_DIR"
    exit 1
fi

echo "SUCCESS: essential NPT input files found."
echo "Equilibration completed. Input files for NPT and NVT simulations have been prepared."
echo "NPT input files are located in: $NPT_INPUT_DIR"
echo "NVT input files are located in: $NVT_INPUT_DIR"