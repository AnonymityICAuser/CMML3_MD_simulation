#!/bin/bash
shopt -s expand_aliases
INPUT_DIR=$1
OUTPUT_DIR=$2


alias gmx="/software/gromacs-2024.3/build/bin/gmx_mpi"

# Files
INPUT_PDB=$INPUT_DIR/datasets/2wwx_clean.pdb
ION_ADD=$INPUT_DIR/required_files/ions.mdp
EM_ADD=$INPUT_DIR/required_files/em.mdp
NVT_DIR=$INPUT_DIR/required_files/NVTs
NPT_DIR=$INPUT_DIR/required_files/NPTs

# Temperature and replicate settings
TEMPERATURES=(320 300 280)
REPLICATES=1

# Create output directory structure
mkdir -p $OUTPUT_DIR



# Run replicates as independent simulations
for rep in $(seq 1 $REPLICATES); do
    REP_DIR=${OUTPUT_DIR}/replicate_${rep}
    mkdir -p ${REP_DIR}/{1_env_set,2_define_box,3_solvent_box,4_add_ions,5_neutralization,6_em,7_nvt,8_npt,9_analysis}
    
    echo "Starting replicate $rep"
    
    # Initial setup (same for all temperatures)
    gmx pdb2gmx -f $INPUT_PDB -o ${REP_DIR}/1_env_set/complex_processed.gro -water spce -ff charmm27 -ignh -missing
    
    gmx editconf -f ${REP_DIR}/1_env_set/complex_processed.gro -o ${REP_DIR}/2_define_box/complex_box.gro -c -d 1.0 -bt cubic
    
    gmx solvate -cp ${REP_DIR}/2_define_box/complex_box.gro -cs spc216.gro -o ${REP_DIR}/3_solvent_box/complex_solv.gro -p topol.top
    
    gmx grompp -f $ION_ADD -c ${REP_DIR}/3_solvent_box/complex_solv.gro -p topol.top -o ${REP_DIR}/4_add_ions/ions.tpr -maxwarn 2
    
    gmx genion -s ${REP_DIR}/4_add_ions/ions.tpr -o ${REP_DIR}/5_neutralization/complex_ions.gro -p topol.top -pname NA -nname CL -neutral 
    
    # Energy minimization
    gmx grompp -f $EM_ADD -c ${REP_DIR}/5_neutralization/complex_ions.gro -p topol.top -o ${REP_DIR}/6_em/em.tpr -maxwarn 2
    
    gmx mdrun -v -deffnm ${REP_DIR}/6_em/em -nb gpu -ntomp 8 -pin on -pinoffset 0
    
    # Run each temperature
    for temp in ${TEMPERATURES[@]}; do
        echo "Running temperature $temp K for replicate $rep"
        
        # Create temperature-specific directories
        mkdir -p ${REP_DIR}/7_nvt/${temp}K
        mkdir -p ${REP_DIR}/8_npt/${temp}K/{1ns,10ns,50ns}
        mkdir -p ${REP_DIR}/9_analysis/${temp}K/{1ns,10ns,50ns}
        
        # NVT equilibration - 200ps
        # Create a modified NVT mdp file with 200ps duration
        cp ${NVT_DIR}/NVT_${temp}K.mdp ${REP_DIR}/7_nvt/${temp}K/NVT_${temp}K_200ps.mdp
        sed -i 's/nsteps.*=.*/nsteps = 100000 ; 200 ps (0.002 * 100000)/g' ${REP_DIR}/7_nvt/${temp}K/NVT_${temp}K_200ps.mdp
        
        gmx grompp -f ${REP_DIR}/7_nvt/${temp}K/NVT_${temp}K_200ps.mdp -c ${REP_DIR}/6_em/em.gro -r ${REP_DIR}/6_em/em.gro -p topol.top -o ${REP_DIR}/7_nvt/${temp}K/nvt.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/7_nvt/${temp}K/nvt -nb gpu -ntomp 8 -pin on -pinoffset 0
        
        echo "NVT equilibration completed for $temp K"
        
        # STAGE 1: NPT short run (500ps)
        gmx grompp -f ${NPT_DIR}/NPT_${temp}K.mdp -c ${REP_DIR}/7_nvt/${temp}K/nvt.gro -r ${REP_DIR}/7_nvt/${temp}K/nvt.gro -t ${REP_DIR}/7_nvt/${temp}K/nvt.cpt -p topol.top -o ${REP_DIR}/8_npt/${temp}K/1ns/npt_500ps.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/8_npt/${temp}K/1ns/npt_500ps -nb gpu -ntomp 8 -pin on -pinoffset 0

        
        # Analysis for 1ns simulation
        ANALYSIS_DIR=${REP_DIR}/9_analysis/${temp}K/1ns
        
        # Process trajectory (remove PBC, center protein)
        echo "Protein System" | gmx trjconv -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.trr \
            -o ${ANALYSIS_DIR}/npt_noPBC.xtc -pbc mol -center -ur compact
        
        # RMSD analysis (structural stability)
        echo "Backbone Backbone" | gmx rms -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/rmsd_backbone.xvg -tu ns
        
        # RMSF analysis (residue flexibility)
        echo "Backbone" | gmx rmsf -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/rmsf.xvg -res
        
        # Radius of gyration (protein compactness)
        echo "Protein" | gmx gyrate -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/gyrate.xvg
        
        # Hydrogen bond analysis
        echo "Protein Protein" | gmx hbond -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -num ${ANALYSIS_DIR}/hbnum.xvg
        
        # Secondary structure analysis
        echo "Protein" | gmx do_dssp -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/ss.xpm -sc ${ANALYSIS_DIR}/scount.xvg
        
        # Solvent accessible surface area
        echo "Protein" | gmx sasa -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/sasa.xvg
        
        # Energy analysis
        echo "Potential" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.edr -o ${ANALYSIS_DIR}/potential.xvg
        echo "Temperature" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.edr -o ${ANALYSIS_DIR}/temperature.xvg
        echo "Pressure" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.edr -o ${ANALYSIS_DIR}/pressure.xvg
        echo "Density" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.edr -o ${ANALYSIS_DIR}/density.xvg
        
        # Create custom index groups if needed
        echo -e "r 1-10\nname 1 ActiveSite\nq" | gmx make_ndx -f ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -o ${ANALYSIS_DIR}/index.ndx
        
        # Principal Component Analysis
        mkdir -p ${ANALYSIS_DIR}/pca
        echo "Backbone" | gmx covar -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/pca/eigenvalues.xvg -v ${ANALYSIS_DIR}/pca/eigenvectors.trr \
            -xpma ${ANALYSIS_DIR}/pca/covariance.xpm
        
        echo "Backbone Backbone" | gmx anaeig -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -v ${ANALYSIS_DIR}/pca/eigenvectors.trr -first 1 -last 2 \
            -proj ${ANALYSIS_DIR}/pca/projection.xvg
        
        # Cluster analysis
        mkdir -p ${ANALYSIS_DIR}/clusters
        echo "Backbone" | gmx cluster -s ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -method gromos -cutoff 0.2 -o ${ANALYSIS_DIR}/clusters/clusters.xpm \
            -cl ${ANALYSIS_DIR}/clusters/cluster_structures.pdb \
            -clid ${ANALYSIS_DIR}/clusters/cluster_ids.xvg
        
        echo "1ns analysis completed. Check results before continuing..."
        echo "Press Enter to continue with 10ns simulation or Ctrl+C to stop..."
        read -t 10 || true  # Wait for user input with 10 second timeout
        
        # STAGE 2: NPT medium run (10ns)
        cp ${NPT_DIR}/NPT_${temp}K.mdp ${REP_DIR}/8_npt/${temp}K/10ns/NPT_${temp}K_10ns.mdp
        sed -i 's/nsteps.*=.*/nsteps = 5000000 ; 10 ns (0.002 * 5000000)/g' ${REP_DIR}/8_npt/${temp}K/10ns/NPT_${temp}K_10ns.mdp
        
        gmx grompp -f ${REP_DIR}/8_npt/${temp}K/10ns/NPT_${temp}K_10ns.mdp -c ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.gro -r ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.gro -t ${REP_DIR}/8_npt/${temp}K/1ns/npt_1ns.cpt -p topol.top -o ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns -nb gpu -ntomp 8 -pin on -pinoffset 0
        
        echo "10ns NPT simulation completed for $temp K"
        
        # Analysis for 10ns simulation
        ANALYSIS_DIR=${REP_DIR}/9_analysis/${temp}K/10ns
        
        # Process trajectory (remove PBC, center protein)
        echo "Protein System" | gmx trjconv -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.trr \
            -o ${ANALYSIS_DIR}/npt_noPBC.xtc -pbc mol -center -ur compact
        
        # RMSD analysis (structural stability)
        echo "Backbone Backbone" | gmx rms -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/rmsd_backbone.xvg -tu ns
        
        # RMSF analysis (residue flexibility)
        echo "Backbone" | gmx rmsf -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/rmsf.xvg -res
        
        # Radius of gyration (protein compactness)
        echo "Protein" | gmx gyrate -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/gyrate.xvg
        
        # Hydrogen bond analysis
        echo "Protein Protein" | gmx hbond -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -num ${ANALYSIS_DIR}/hbnum.xvg
        
        # Secondary structure analysis
        echo "Protein" | gmx do_dssp -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/ss.xpm -sc ${ANALYSIS_DIR}/scount.xvg
        
        # Solvent accessible surface area
        echo "Protein" | gmx sasa -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/sasa.xvg
        
        # Energy analysis
        echo "Potential" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.edr -o ${ANALYSIS_DIR}/potential.xvg
        echo "Temperature" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.edr -o ${ANALYSIS_DIR}/temperature.xvg
        echo "Pressure" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.edr -o ${ANALYSIS_DIR}/pressure.xvg
        echo "Density" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.edr -o ${ANALYSIS_DIR}/density.xvg
        
        # Create custom index groups if needed
        echo -e "r 1-10\nname 1 ActiveSite\nq" | gmx make_ndx -f ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -o ${ANALYSIS_DIR}/index.ndx
        
        # Principal Component Analysis
        mkdir -p ${ANALYSIS_DIR}/pca
        echo "Backbone" | gmx covar -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/pca/eigenvalues.xvg -v ${ANALYSIS_DIR}/pca/eigenvectors.trr \
            -xpma ${ANALYSIS_DIR}/pca/covariance.xpm
        
        echo "Backbone Backbone" | gmx anaeig -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -v ${ANALYSIS_DIR}/pca/eigenvectors.trr -first 1 -last 2 \
            -proj ${ANALYSIS_DIR}/pca/projection.xvg
        
        # Cluster analysis
        mkdir -p ${ANALYSIS_DIR}/clusters
        echo "Backbone" | gmx cluster -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -method gromos -cutoff 0.2 -o ${ANALYSIS_DIR}/clusters/clusters.xpm \
            -cl ${ANALYSIS_DIR}/clusters/cluster_structures.pdb \
            -clid ${ANALYSIS_DIR}/clusters/cluster_ids.xvg
        
        # Contact map analysis
        echo "Backbone" | gmx mdmat -s ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -mean ${ANALYSIS_DIR}/contact_map.xpm
        
        echo "10ns analysis completed. Check results before continuing..."
        echo "Press Enter to continue with 50ns simulation or Ctrl+C to stop..."
        read -t 10 || true  # Wait for user input with 10 second timeout
        
        # STAGE 3: NPT full run (50ns)
        cp ${NPT_DIR}/NPT_${temp}K.mdp ${REP_DIR}/8_npt/${temp}K/50ns/NPT_${temp}K_50ns.mdp
        sed -i 's/nsteps.*=.*/nsteps = 25000000 ; 50 ns (0.002 * 25000000)/g' ${REP_DIR}/8_npt/${temp}K/50ns/NPT_${temp}K_50ns.mdp
        
        gmx grompp -f ${REP_DIR}/8_npt/${temp}K/50ns/NPT_${temp}K_50ns.mdp -c ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.gro -r ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.gro -t ${REP_DIR}/8_npt/${temp}K/10ns/npt_10ns.cpt -p topol.top -o ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -maxwarn 2
        
        gmx mdrun -v -deffnm ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns -nb gpu -ntomp 8 -pin on -pinoffset 0
        
        echo "50ns NPT simulation completed for $temp K"
        
        # Analysis for 50ns simulation
        ANALYSIS_DIR=${REP_DIR}/9_analysis/${temp}K/50ns
        
        # Process trajectory (remove PBC, center protein)
        echo "Protein System" | gmx trjconv -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.trr \
            -o ${ANALYSIS_DIR}/npt_noPBC.xtc -pbc mol -center -ur compact
        
        # RMSD analysis (structural stability)
        echo "Backbone Backbone" | gmx rms -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/rmsd_backbone.xvg -tu ns
        
        # RMSF analysis (residue flexibility)
        echo "Backbone" | gmx rmsf -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/rmsf.xvg -res
        
        # Radius of gyration (protein compactness)
        echo "Protein" | gmx gyrate -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/gyrate.xvg
        
        # Hydrogen bond analysis
        echo "Protein Protein" | gmx hbond -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -num ${ANALYSIS_DIR}/hbnum.xvg
        
        # Secondary structure analysis
        echo "Protein" | gmx do_dssp -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/ss.xpm -sc ${ANALYSIS_DIR}/scount.xvg
        
        # Solvent accessible surface area
        echo "Protein" | gmx sasa -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/sasa.xvg
        
        # Energy analysis
        echo "Potential" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.edr -o ${ANALYSIS_DIR}/potential.xvg
        echo "Temperature" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.edr -o ${ANALYSIS_DIR}/temperature.xvg
        echo "Pressure" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.edr -o ${ANALYSIS_DIR}/pressure.xvg
        echo "Density" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.edr -o ${ANALYSIS_DIR}/density.xvg
        
        # Create custom index groups if needed
        echo -e "r 1-10\nname 1 ActiveSite\nq" | gmx make_ndx -f ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -o ${ANALYSIS_DIR}/index.ndx
        
        # Principal Component Analysis
        mkdir -p ${ANALYSIS_DIR}/pca
        echo "Backbone" | gmx covar -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -o ${ANALYSIS_DIR}/pca/eigenvalues.xvg -v ${ANALYSIS_DIR}/pca/eigenvectors.trr \
            -xpma ${ANALYSIS_DIR}/pca/covariance.xpm
        
        echo "Backbone Backbone" | gmx anaeig -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -v ${ANALYSIS_DIR}/pca/eigenvectors.trr -first 1 -last 2 \
            -proj ${ANALYSIS_DIR}/pca/projection.xvg
        
        echo "Backbone" | gmx anaeig -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -v ${ANALYSIS_DIR}/pca/eigenvectors.trr -first 1 -last 1 -extr ${ANALYSIS_DIR}/pca/extreme_pc1.pdb
        
        # Cluster analysis
        mkdir -p ${ANALYSIS_DIR}/clusters
        echo "Backbone" | gmx cluster -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -method gromos -cutoff 0.2 -o ${ANALYSIS_DIR}/clusters/clusters.xpm \
            -cl ${ANALYSIS_DIR}/clusters/cluster_structures.pdb \
            -clid ${ANALYSIS_DIR}/clusters/cluster_ids.xvg
        
        # Contact map analysis
        echo "Backbone" | gmx mdmat -s ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.tpr -f ${ANALYSIS_DIR}/npt_noPBC.xtc \
            -mean ${ANALYSIS_DIR}/contact_map.xpm
        
        # Create symbolic links for the final result in the main NPT directory for convenience
        ln -sf ${REP_DIR}/8_npt/${temp}K/50ns/npt_50ns.* ${REP_DIR}/8_npt/${temp}K/
        
        echo "Temperature $temp K completed for replicate $rep"
    done
    
    echo "Replicate $rep completed"
done

# Optional: Generate combined analysis across replicates and temperatures
mkdir -p ${OUTPUT_DIR}/combined_analysis/{1ns,10ns,50ns}/{rmsd,rmsf,gyrate,hbond,energy}

# Process each time point separately
for duration in 1ns 10ns 50ns; do
    # Extract and combine RMSD data from all replicates
    for temp in ${TEMPERATURES[@]}; do
        touch ${OUTPUT_DIR}/combined_analysis/${duration}/rmsd/combined_rmsd_${temp}K.dat
        echo "# Time(ns) RMSD(nm) Replicate" > ${OUTPUT_DIR}/combined_analysis/${duration}/rmsd/combined_rmsd_${temp}K.dat
        
        for rep in $(seq 1 $REPLICATES); do
            # Extract time and RMSD columns, add replicate number
            grep -v "@" ${OUTPUT_DIR}/replicate_${rep}/9_analysis/${temp}K/${duration}/rmsd_backbone.xvg | \
            grep -v "#" | awk -v rep=$rep '{print $1, $2, rep}' >> ${OUTPUT_DIR}/combined_analysis/${duration}/rmsd/combined_rmsd_${temp}K.dat
        done
    done
    
    # Extract and combine RMSF data
    for temp in ${TEMPERATURES[@]}; do
        touch ${OUTPUT_DIR}/combined_analysis/${duration}/rmsf/combined_rmsf_${temp}K.dat
        echo "# ResID RMSF_Rep1 RMSF_Rep2 RMSF_Rep3" > ${OUTPUT_DIR}/combined_analysis/${duration}/rmsf/combined_rmsf_${temp}K.dat
        
        # Get residue IDs from first replicate
        grep -v "@" ${OUTPUT_DIR}/replicate_1/9_analysis/${temp}K/${duration}/rmsf.xvg | \
        grep -v "#" | awk '{print $1}' > ${OUTPUT_DIR}/combined_analysis/${duration}/rmsf/resids.tmp
        
        # Get RMSF values for each replicate
        for rep in $(seq 1 $REPLICATES); do
            grep -v "@" ${OUTPUT_DIR}/replicate_${rep}/9_analysis/${temp}K/${duration}/rmsf.xvg | \
            grep -v "#" | awk '{print $2}' > ${OUTPUT_DIR}/combined_analysis/${duration}/rmsf/rmsf_rep${rep}.tmp
        done
        
        # Combine data
        paste ${OUTPUT_DIR}/combined_analysis/${duration}/rmsf/resids.tmp \
              ${OUTPUT_DIR}/combined_analysis/${duration}/rmsf/rmsf_rep*.tmp \
              >> ${OUTPUT_DIR}/combined_analysis/${duration}/rmsf/combined_rmsf_${temp}K.dat
        
        # Clean up temporary files
        rm ${OUTPUT_DIR}/combined_analysis/${duration}/rmsf/*.tmp
    done
    
    # Generate average properties across replicates
    for temp in ${TEMPERATURES[@]}; do
        # Average radius of gyration
        touch ${OUTPUT_DIR}/combined_analysis/${duration}/gyrate/avg_gyrate_${temp}K.dat
        echo "# Time(ns) Avg_Rg StdDev_Rg" > ${OUTPUT_DIR}/combined_analysis/${duration}/gyrate/avg_gyrate_${temp}K.dat
        
        # Process gyration data for each time point
        for rep in $(seq 1 $REPLICATES); do
            grep -v "@" ${OUTPUT_DIR}/replicate_${rep}/9_analysis/${temp}K/${duration}/gyrate.xvg | \
            grep -v "#" | awk '{print $1, $2}' > ${OUTPUT_DIR}/combined_analysis/${duration}/gyrate/gyrate_rep${rep}.tmp
        done
        
        # Use awk to calculate average and standard deviation for each time point
        awk '
        BEGIN {rep_count='"$REPLICATES"'}
        NR==FNR {time[$1]=$1; rg[1][$1]=$2; next}
        {
            file_num=ARGIND;
            if (file_num <= rep_count) {
                rg[file_num][$1]=$2;
            }
        }
        END {
            for (t in time) {
                sum=0; sumsq=0;
                for (r=1; r<=rep_count; r++) {
                    sum += rg[r][t];
                    sumsq += rg[r][t]^2;
                }
                avg = sum/rep_count;
                stddev = sqrt(sumsq/rep_count - avg^2);
                printf "%f %f %f\n", t, avg, stddev;
            }
        }
        ' ${OUTPUT_DIR}/combined_analysis/${duration}/gyrate/gyrate_rep*.tmp | sort -n > ${OUTPUT_DIR}/combined_analysis/${duration}/gyrate/avg_gyrate_${temp}K.dat
        
        # Clean up temporary files
        rm ${OUTPUT_DIR}/combined_analysis/${duration}/gyrate/*.tmp
    done
    
