#!/bin/bash
shopt -s expand_aliases
INPUT_DIR=$1
OUTPUT_DIR=$2

alias gmx="/software/gromacs-2024.3/build/bin/gmx_mpi"

# Files
INPUT_PDB=$INPUT_DIR/datasets/2wwx_clean.pdb
ION_ADD=$INPUT_DIR/required_files/ions.mdp
EM_ADD=$INPUT_DIR/required_files/em.mdp
NVT_BASE=$INPUT_DIR/required_files/NVT.mdp
NPT_BASE=$INPUT_DIR/required_files/NPT.mdp
NVT_DIR=$INPUT_DIR/required_files/NVTs
NPT_DIR=$INPUT_DIR/required_files/NPTs
MD_1NS=$INPUT_DIR/required_files/md_1ns.mdp
MD_10NS=$INPUT_DIR/required_files/md_10ns.mdp
MD_50NS=$INPUT_DIR/required_files/md_50ns.mdp

# Temperature and replicate settings
TEMPERATURES=(300 320 280)
REPLICATES=1

# Create output directory structure
mkdir -p $OUTPUT_DIR

# Run replicates as independent simulations
for rep in $(seq 1 $REPLICATES); do
    REP_DIR=${OUTPUT_DIR}/replicate_${rep}
    mkdir -p ${REP_DIR}/{1_env_set,2_define_box,3_solvent_box,4_add_ions,5_neutralization,6_em,7_nvt,8_npt,9_md}
    
    echo "Starting replicate $rep"
    
    # Create a working directory for this replicate
    WORK_DIR=${REP_DIR}/working
    mkdir -p $WORK_DIR
    cd $WORK_DIR
    
    # Initial setup (same for all temperatures)
    echo "8" | gmx pdb2gmx -f $INPUT_PDB -o ${REP_DIR}/1_env_set/complex_processed.gro -water spce -ignh -missing
    
    # Copy topology files to the env_set directory
    cp *.top *.itp ${REP_DIR}/1_env_set/
    
    gmx editconf -f ${REP_DIR}/1_env_set/complex_processed.gro -o ${REP_DIR}/2_define_box/complex_box.gro -c -d 1.0 -bt cubic
    
    gmx solvate -cp ${REP_DIR}/2_define_box/complex_box.gro -cs spc216.gro -o ${REP_DIR}/3_solvent_box/complex_solv.gro -p topol.top
    
    # Copy updated topology to solvent_box directory
    cp topol.top ${REP_DIR}/3_solvent_box/
    
    gmx grompp -f $ION_ADD -c ${REP_DIR}/3_solvent_box/complex_solv.gro -p topol.top -o ${REP_DIR}/4_add_ions/ions.tpr -maxwarn 2
    
    echo "13" | gmx genion -s ${REP_DIR}/4_add_ions/ions.tpr -o ${REP_DIR}/5_neutralization/complex_ions.gro -p topol.top -neutral 
    
    # Copy updated topology to neutralization directory
    cp topol.top ${REP_DIR}/5_neutralization/
    
    # Energy minimization
    gmx grompp -f $EM_ADD -c ${REP_DIR}/5_neutralization/complex_ions.gro -p topol.top -o ${REP_DIR}/6_em/em.tpr -maxwarn 2
    
    gmx mdrun -v -deffnm ${REP_DIR}/6_em/em -nb gpu -ntomp 8 -pin on -pinoffset 0
    
    # Copy updated topology to em directory
    cp topol.top ${REP_DIR}/6_em/
    
    
    # Run each temperature
    for temp in ${TEMPERATURES[@]}; do
        echo "Running temperature $temp K for replicate $rep"
        
        # Create temperature-specific directories
        mkdir -p ${REP_DIR}/7_nvt/${temp}K
        mkdir -p ${REP_DIR}/8_npt/${temp}K
        mkdir -p ${REP_DIR}/9_md/${temp}K/{1ns,10ns,50ns}
        mkdir -p ${REP_DIR}/10_analysis/${temp}K/{1ns,10ns,50ns}
        
        # NVT equilibration - 200ps
        gmx grompp -f ${NVT_DIR}/NVT_${temp}K.mdp -c ${REP_DIR}/6_em/em.gro -r ${REP_DIR}/6_em/em.gro -p topol.top -o ${REP_DIR}/7_nvt/${temp}K/nvt.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/7_nvt/${temp}K/nvt -nb gpu -ntomp 8 -pin on -pinoffset 0
        
        # Copy topology to nvt directory
        cp topol.top ${REP_DIR}/7_nvt/${temp}K/
        
        echo "NVT equilibration completed for $temp K"
        
        # NPT equilibration (500ps)
        gmx grompp -f ${NPT_DIR}/NPT_${temp}K.mdp -c ${REP_DIR}/7_nvt/${temp}K/nvt.gro -r ${REP_DIR}/7_nvt/${temp}K/nvt.gro -t ${REP_DIR}/7_nvt/${temp}K/nvt.cpt -p topol.top -o ${REP_DIR}/8_npt/${temp}K/npt.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/8_npt/${temp}K/npt -nb gpu -ntomp 8 -pin on -pinoffset 0
        
        # Copy topology to npt directory
        cp topol.top ${REP_DIR}/8_npt/${temp}K/
        
        # Production MD runs at different time scales
        # 1ns run
        gmx grompp -f $MD_1NS -c ${REP_DIR}/8_npt/${temp}K/npt.gro -t ${REP_DIR}/8_npt/${temp}K/npt.cpt -p topol.top -o ${REP_DIR}/9_md/${temp}K/1ns/md_1ns.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/9_md/${temp}K/1ns/md_1ns -nb gpu -ntomp 8 -pin on -pinoffset 0
        
        # Copy topology to md directory
        cp topol.top ${REP_DIR}/9_md/${temp}K/1ns/
        
        # Process trajectory for analysis (remove PBC)
        echo "1 0" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/1ns/md_1ns.tpr -f ${REP_DIR}/9_md/${temp}K/1ns/md_1ns.xtc -o ${REP_DIR}/9_md/${temp}K/1ns/md_1ns_noPBC.xtc -pbc mol -center

        # 10ns run (only if needed - can comment out to save time initially)
        gmx grompp -f $MD_10NS -c ${REP_DIR}/9_md/${temp}K/1ns/md_1ns.gro -t ${REP_DIR}/9_md/${temp}K/1ns/md_1ns.cpt -p topol.top -o ${REP_DIR}/9_md/${temp}K/10ns/md_10ns.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/9_md/${temp}K/10ns/md_10ns -nb gpu -ntomp 10 -pin on -pinoffset 0
        
        # Copy topology to 10ns directory
        cp topol.top ${REP_DIR}/9_md/${temp}K/10ns/
        
        # Process trajectory for 10ns analysis
        echo "1 0" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/10ns/md_10ns.tpr -f ${REP_DIR}/9_md/${temp}K/10ns/md_10ns.xtc -o ${REP_DIR}/9_md/${temp}K/10ns/md_10ns_noPBC.xtc -pbc mol -center
        
        # 50ns run (only if needed - can comment out to save time initially)
        gmx grompp -f $MD_50NS -c ${REP_DIR}/9_md/${temp}K/10ns/md_10ns.gro -t ${REP_DIR}/9_md/${temp}K/10ns/md_10ns.cpt -p topol.top -o ${REP_DIR}/9_md/${temp}K/50ns/md_50ns.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/9_md/${temp}K/50ns/md_50ns -nb gpu -ntomp 10 -pin on -pinoffset 0
        
        # Copy topology to 50ns directory
        cp topol.top ${REP_DIR}/9_md/${temp}K/50ns/
        
        # Process trajectory for 50ns analysis
        echo "1 0" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/50ns/md_50ns.tpr -f ${REP_DIR}/9_md/${temp}K/50ns/md_50ns.xtc -o ${REP_DIR}/9_md/${temp}K/50ns/md_50ns_noPBC.xtc -pbc mol -center

        echo "Completed all simulations for $temp K, replicate $rep"
    done
    
    # Clean up working directory
    cd $OUTPUT_DIR
    
    echo "Completed replicate $rep"
done

echo "All simulations completed successfully!"