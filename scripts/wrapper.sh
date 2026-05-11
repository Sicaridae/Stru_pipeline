#!/bin/bash



EDR_DIR="/home/pesce/bio_strutturale/pipe_line/input/1a"
script_dirA="/home/pesce/bio_strutturale/pipe_line/scripts/1a"
out_dirA="/home/pesce/bio_strutturale/pipe_line/results/1a"

INPUT_1B="/home/pesce/bio_strutturale/pipe_line/input/1b"
script_dirB="/home/pesce/bio_strutturale/pipe_line/scripts/1b"
data_dirB="/home/pesce/bio_strutturale/pipe_line/input/1b"
out_dirB="/home/pesce/bio_strutturale/pipe_line/results/1b"
script_slurm="/home/pesce/bio_strutturale/pipe_line/scripts/script_cluster"

DIR_1A="/home/pesce/bio_strutturale/pipe_line/input/1a"
DIR_1B="/home/pesce/bio_strutturale/pipe_line/input/1b"


slurm=0
slurm_e=0
slurm_npt=0  

slurm_nvt=0
A1=0
B1=0

### this blok is for running the simulation. 
### as input takes the .mdp, index and topology files, and the initial structure.
### the user can choose to run the whole pipeline (from equilibration to nvt) or only some steps.
### the pipeline alow the user to stops and check the outputs

if [ "$slurm" -eq 1 ]; then

    if [ "$slurm_e" -eq 1 ] && [ "$slurm_npt" -eq 1 ] && [ "$slurm_nvt" -eq 1 ]; then
        echo "Launching equilibration on SLURM"

        jid_eq=$(sbatch "$script_slurm/equilibration.slurm" | awk '{print $4}')

        if [ -z "$jid_eq" ]; then
            echo "ERROR: equilibration job submission failed"
            exit 1
        fi

        echo "Equilibration submitted with job ID: $jid_eq"

        echo "Launching NPT on SLURM after equilibration"
        jid_npt=$(sbatch --dependency=afterok:$jid_eq "$script_slurm/npt_3ns.slurm" | awk '{print $4}')

        if [ -z "$jid_npt" ]; then
            echo "ERROR: NPT job submission failed"
            exit 1
        fi

        echo "NPT submitted with job ID: $jid_npt"
        echo "NPT will start only if equilibration finishes correctly."

        echo "Launching NVT on SLURM after NPT"
        jid_nvt=$(sbatch --dependency=afterok:$jid_npt "$script_slurm/nvt_3ns.slurm" | awk '{print $4}')

        if [ -z "$jid_nvt" ]; then
            echo "ERROR: NVT job submission failed"
            exit 1
        fi

        echo "NVT submitted with job ID: $jid_nvt"
        echo "NVT will start only if NPT finishes correctly."

    elif [ "$slurm_e" -eq 0 ] && [ "$slurm_npt" -eq 1 ] && [ "$slurm_nvt" -eq 1 ]; then
        echo "Launching NPT on SLURM"

        jid_npt=$(sbatch "$script_slurm/npt_3ns.slurm" | awk '{print $4}')

        if [ -z "$jid_npt" ]; then
            echo "ERROR: NPT job submission failed"
            exit 1
        fi

        echo "NPT submitted with job ID: $jid_npt"

        echo "Launching NVT on SLURM after NPT"
        jid_nvt=$(sbatch --dependency=afterok:$jid_npt "$script_slurm/nvt_3ns.slurm" | awk '{print $4}')

        if [ -z "$jid_nvt" ]; then
            echo "ERROR: NVT job submission failed"
            exit 1
        fi

        echo "NVT submitted with job ID: $jid_nvt"

    elif [ "$slurm_e" -eq 0 ] && [ "$slurm_npt" -eq 0 ] && [ "$slurm_nvt" -eq 1 ]; then
        echo "Launching NVT on SLURM"

        jid_nvt=$(sbatch "$script_slurm/nvt_3ns.slurm" | awk '{print $4}')

        if [ -z "$jid_nvt" ]; then
            echo "ERROR: NVT job submission failed"
            exit 1
        fi

        echo "NVT submitted with job ID: $jid_nvt"
    fi
fi
    

#### this block is for the analysis of total energy, potential energy, temperature and pressure
#### It takes also the extraction from the nvt .edr of this parameters and the parsing of the .xvg filesautomatically
#### a small control is within the python script to check if ther're some empty lines in the parsed file, if so the script will stop and print an error message.

if [ "$A1" -eq 1 ]; then
    echo "Running block 1A"
    
    for file in "$EDR_DIR"/*.edr; do
        base="${file%.edr}"
        echo -e "11\n13\n15\n16\n0" | gmx energy -f "$file" -o "${base}_energy.xvg"
    done
    
    cd "$EDR_DIR"

    for file in "$EDR_DIR"/*.xvg; do
        base="${file%.xvg}"
        parsed_file="${base}_parsed.txt"
        echo "$parsed_file"
        grep -v '^[#@]' "$file" > "$parsed_file"
        python3 "$script_dirA/script1.py" "$parsed_file" "$out_dirA"    #### 
done
fi

#### this block takes care of th structural analysis of the simulations (NPT and NVT)
#### the python script is made of a chain of functions that take care of the analysis of the trajectory


if [ "$B1" -eq 1 ]; then
    echo "Running block 1B"

    for sim_dir in "$data_dirB"/*; do

        if [ -d "$sim_dir" ]; then
            sim_name="$(basename "$sim_dir")"
            sim_out_dir="$out_dirB/$sim_name"

            mkdir -p "$sim_out_dir"

            tpr_file="$sim_dir/${sim_name}_3ns.tpr"
            xtc_file="$sim_dir/${sim_name}_3ns.xtc"

            echo "Analyzing simulation: $sim_name"
            echo "TPR file: $tpr_file"
            echo "XTC file: $xtc_file"
            echo "Output folder: $sim_out_dir"

            python3 "$script_dirB/full_1b.py" "$tpr_file" "$xtc_file" "$sim_out_dir" "$sim_name"
        fi
    done
fi
