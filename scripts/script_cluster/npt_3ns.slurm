#!/bin/bash
#SBATCH --job-name=kcnq1_npt
#SBATCH --output=/home/studmol/yy2025/jfedeli/data/slurm/kcnq1_npt.out
#SBATCH --error=/home/studmol/yy2025/jfedeli/data/slurm/kcnq1_npt.err
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

set -e

module load molecular-dynamics
module load gromacs/2023.5


buffer="/home/pesce/bio_strutturale/pipe_line/output/buffer/npt_input"
output="/home/pesce/bio_strutturale/pipe_line/output/buffer/nvt_input"
GROMACS_DIR="/home/pesce/bio_strutturale/pipe_line/chermm_steps/gromacs"
DIR_1A="/home/pesce/bio_strutturale/pipe_line/input/1a"
DIR_1B="/home/pesce/bio_strutturale/pipe_line/input/1b"

mkdir -p "$DIR_1A" "$DIR_1B"
mkdir -p "$output"
cd "$buffer"

gmx grompp \
  -f "$buffer/npt_3ns.mdp" \
  -o "$buffer/npt_3ns.tpr" \
  -c "$buffer/step6.6_equilibration.gro" \
  -r "$buffer/step6.6_equilibration.gro" \
  -t "$buffer/step6.6_equilibration.cpt" \
  -p "$buffer/topol.top" \
  -n "$buffer/index.ndx" -maxwarn 1

gmx mdrun -deffnm npt_3ns -ntomp 8

mkdir -p "$DIR_1A/{npt,nvt}" "$DIR_1B/{npt,nvt}"

cp -n "$buffer/npt_3ns.edr" "$DIR_1A/npt/"
cp -n "$buffer/npt_3ns.tpr" "$DIR_1B/npt/"
cp -n "$buffer/npt_3ns.xtc" "$DIR_1B/npt/"

cp -n "$buffer/npt_3ns.gro" "$output/"
cp -n "$buffer/npt_3ns.cpt" "$output/"
cp -n "$GROMACS_DIR/topol.top" "$output/"
cp -n "$GROMACS_DIR/index.ndx" "$output/"
cp -n "$GROMACS_DIR/nvt_3ns.mdp" "$output/"

### check of the existence of the essential files for the NPT simulation, if they are missing, exit with an error message
if [ ! -s "$output/npt_3ns.gro" ] || [ ! -s "$output/npt_3ns.cpt" ]; then
    echo "ERROR: missing essential NPT input files in $output"
    exit 1
fi

echo "SUCCESS: essential NVT input files found."
echo "NPT completed. Input files for NVT simulations have been prepared."
echo "NVT input files are located in: $output"