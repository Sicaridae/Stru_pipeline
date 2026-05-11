#!/bin/bash
#SBATCH --job-name=kcnq1_nvt
#SBATCH --output=/home/studmol/yy2025/jfedeli/data/slurm/kcnq1_nvt.out
#SBATCH --error=/home/studmol/yy2025/jfedeli/data/slurm/kcnq1_nvt.err
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

set -e

module load molecular-dynamics
module load gromacs/2023.5

buffer="/home/pesce/bio_strutturale/pipe_line/output/buffer/nvt_input"

DIR_1A="/home/pesce/bio_strutturale/pipe_line/input/1a"
DIR_1B="/home/pesce/bio_strutturale/pipe_line/input/1b"
#### would be awesome to add a cp -n comand that copi from the gromacs directory 
#### index and nvt.mdp file, but now i'm tired

mkdir -p "$DIR_1A" "$DIR_1B"
cd "$buffer"

gmx grompp \
  -f ${buffer}/nvt_3ns.mdp \
  -o ${buffer}/nvt_3ns.tpr \
  -c ${buffer}/npt_3ns.gro \
  -r ${buffer}/npt_3ns.gro \
  -t ${buffer}/npt_3ns.cpt \
  -p ${buffer}/topol.top \
  -n ${buffer}/index.ndx \
  -maxwarn 1

gmx mdrun \
  -deffnm nvt_3ns \
  -ntomp 8



####need to send into 1a and 1b inputs

mkdir -p "$DIR_1A/{npt,nvt}" "$DIR_1B/{npt,nvt}"

cp -n "$buffer/nvt_3ns.edr" "$DIR_1A/nvt/"
cp -n "$buffer/nvt_3ns.tpr" "$DIR_1B/nvt/"
cp -n "$buffer/nvt_3ns.xtc" "$DIR_1B/nvt/"

if [ ! -s "$DIR_1A/nvt/nvt_3ns.edr" ] || [ ! -s "$DIR_1B/nvt/nvt_3ns.tpr" ] || [ ! -s "$DIR_1B/nvt/nvt_3ns.xtc" ]; then
    echo "ERROR: missing essential input files in $DIR_1A or $DIR_1B"
    exit 1
fi

echo "SUCCESS: essential input files found."
echo "NVT completed. Input files for next stage have been prepared."
echo "NVT input files are located in: $DIR_1A and $DIR_1B"