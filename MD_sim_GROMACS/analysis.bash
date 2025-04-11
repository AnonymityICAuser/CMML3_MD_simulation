#!/bin/bash
# Run this after simulation - Combined whole complex and chain-specific analysis with robust PBC handling and XVG formatting
shopt -s expand_aliases
alias gmx="/software/gromacs-2024.3/build/bin/gmx_mpi"

BASE_DIR="work_space"  
REP_DIR="replicate_1"         

# Create analysis directory structure
mkdir -p ${REP_DIR}/result_analysis
mkdir -p ${REP_DIR}/result_analysis/complex/{280K,300K,320K}/{1ns,10ns,50ns}
mkdir -p ${REP_DIR}/result_analysis/chains/ChainH/{280K,300K,320K}/{1ns,10ns,50ns}
mkdir -p ${REP_DIR}/result_analysis/chains/ChainD/{280K,300K,320K}/{1ns,10ns,50ns}
mkdir -p ${REP_DIR}/result_analysis/chains/interface/{280K,300K,320K}/{1ns,10ns,50ns}

# Create index files with chain identifiers
cat > make_ndx_input.txt << EOF
chain A
name 10 ChainH
chain B
name 11 ChainD
r 1-999 & a CA
name 12 C_alpha
10 & a CA
name 13 ChainH_CA
11 & a CA
name 14 ChainD_CA
q
EOF

# Create the index file with chain-based groups
gmx make_ndx -f ${BASE_DIR}/datasets/2wwx_clean.pdb -o ${REP_DIR}/chains.ndx < make_ndx_input.txt
rm make_ndx_input.txt

echo "11 0" | gmx energy -f ${REP_DIR}/6_em/em.edr -o ${REP_DIR}/result_analysis/potential.xvg

# Function to modify XVG files
modify_xvg() {
    local xvg_file=$1
    local title=$2
    local xlabel=$3
    local ylabel=$4
    
    # Use sed to modify the XVG file
    sed -i "s/@    title.*/@    title \"${title}\"/" "${xvg_file}"
    sed -i "s/@    xaxis  label.*/@    xaxis  label \"${xlabel}\"/" "${xvg_file}"
    sed -i "s/@    yaxis  label.*/@    yaxis  label \"${ylabel}\"/" "${xvg_file}"
    
    echo "Modified: ${xvg_file}"
}

# Function to annotate high RMSF residues
annotate_rmsf() {
    local rmsf_file=$1
    local output_file=$2
    local threshold=0.3 # nm
    
    awk -v thresh=$threshold '{
        if (NF >= 2 && $2 > thresh && $1 !~ /^[@#]/) 
            print "Residue "$1" RMSF: "$2" nm"
    }' "${rmsf_file}" > "${output_file}"
    
    echo "High RMSF residues annotated in: ${output_file}"
}

# Modify potential energy XVG
modify_xvg "${REP_DIR}/result_analysis/potential.xvg" "Potential Energy During Energy Minimization" "Time (ps)" "Potential Energy (kJ/mol)"

for temp in  300 320; do  
    echo "Analyzing temperature ${temp}K data..."
    
    # NVT temperature
    echo "17 0" | gmx energy -f ${REP_DIR}/7_nvt/${temp}K/nvt.edr -o ${REP_DIR}/result_analysis/complex/${temp}K/temperature.xvg
    modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/temperature.xvg" "System Temperature at ${temp}K" "Time (ps)" "Temperature (K)"
    
    # NPT pressure and density
    echo "19 0" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/npt.edr -o ${REP_DIR}/result_analysis/complex/${temp}K/pressure.xvg
    modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/pressure.xvg" "System Pressure at ${temp}K" "Time (ps)" "Pressure (bar)"
    
    echo "25 0" | gmx energy -f ${REP_DIR}/8_npt/${temp}K/npt.edr -o ${REP_DIR}/result_analysis/complex/${temp}K/density.xvg
    modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/density.xvg" "System Density at ${temp}K" "Time (ps)" "Density (kg/m³)"
    
    for sim_time in 1 10 50; do
        echo "Processing ${sim_time}ns simulation data for ${temp}K..."

        # Improved trajectory processing with robust PBC handling
        # 1. Make molecules whole using "System" group (Group 0)
        echo "1
0" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.xtc \
            -o ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_whole.xtc \
            -pbc atom -center
        
        # 2. Remove jumps across periodic boundaries using "Protein" group (Group 1)
        echo "1
1" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_whole.xtc \
            -o ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_nojump.xtc \
            -pbc nojump
        
        # 3. Cluster the trajectory to ensure specific parts remain intact
        echo "1
1" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_nojump.xtc \
            -o ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_cluster.xtc \
            -pbc cluster
        
        # 4. Center the trajectory with the reference group
        echo "1
1
1" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_cluster.xtc \
            -o ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -center -fit rot+trans
        
        # Debugging step: Check trajectory integrity
        gmx check -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc
        
        # 5. Generate protein-only GRO file from the first frame
        echo "1" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -o ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein.gro \
            -dump 0
        
        # 6. Generate protein-only GRO file from the last frame
        echo "1" | gmx trjconv -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -o ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein_final.gro \
            -dump ${sim_time}
        
        gmx editconf -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein_final.gro -o ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein_final.pdb 

        # 7. Create protein-only TPR file
        echo "1" | gmx convert-tpr -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -o ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein.tpr
        
        ###################
        # COMPLEX ANALYSIS
        ###################
        echo "Analyzing whole complex..."
        
        # RMSD of whole complex
        echo "4 4" | gmx rms -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -o ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/rmsd.xvg -tu ns
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/rmsd.xvg" "RMSD of Protein Complex at ${temp}K (${sim_time}ns)" "Time (ns)" "RMSD (nm)"
        
        echo "4 4" | gmx rms -s ${REP_DIR}/6_em/em.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -o ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/rmsd_xtal.xvg -tu ns
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/rmsd_xtal.xvg" "RMSD vs Initial Structure at ${temp}K (${sim_time}ns)" "Time (ns)" "RMSD (nm)"
        
        # Gyration radius of whole complex
        echo "1" | gmx gyrate -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -o ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/gyrate.xvg

        # RMSF of whole complex 
        echo "1" | gmx rmsf -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/rmsf.xvg -oq ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/rmsf_bfac.pdb -res
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/rmsf.xvg" "RMSF Per Residue at ${temp}K (${sim_time}ns)" "Residue Number" "RMSF (nm)"
        annotate_rmsf "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/rmsf.xvg" "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/high_rmsf.txt"
        
        # Hydrogen bonds analysis (all and backbone)
        printf "1\n1\n" | gmx hbond -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -tu ns \
            -num ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/hb_all.xvg
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/hb_all.xvg" "Total Hydrogen Bonds at ${temp}K (${sim_time}ns)" "Time (ns)" "Number of H-bonds"
        
        printf "7\n7\n" | gmx hbond -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -tu ns \
            -num ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/hb_bb.xvg
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/hb_bb.xvg" "Backbone Hydrogen Bonds at ${temp}K (${sim_time}ns)" "Time (ns)" "Number of H-bonds"
        
        # Secondary structure analysis (DSSP)
        echo "7" | gmx dssp -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -tu ns \
            -hmode dssp \
            -o ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/dssp.dat \
            -num ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/num.xvg
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/num.xvg" "Secondary Structure Elements at ${temp}K (${sim_time}ns)" "Time (ns)" "Number of Residues"
        
        #################
        # CHAIN ANALYSIS
        #################
        
        # Chain H analysis (Group 10)
        echo "Analyzing Chain H..."
        
        # RMSD for Chain H
        echo "10 
        10" | gmx rms -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/rmsd.xvg -tu ns
        modify_xvg "${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/rmsd.xvg" "RMSD of Chain H at ${temp}K (${sim_time}ns)" "Time (ns)" "RMSD (nm)"
        
        echo "10 
        10" | gmx rms -s ${REP_DIR}/6_em/em.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/rmsd_xtal.xvg -tu ns
        modify_xvg "${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/rmsd_xtal.xvg" "RMSD of Chain H vs Initial Structure at ${temp}K (${sim_time}ns)" "Time (ns)" "RMSD (nm)"
        
        # Gyration radius for Chain H
        echo "10" | gmx gyrate -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/gyrate.xvg

        # RMSF for Chain H (CA atoms only, per residue)
        echo "13" | gmx rmsf -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/rmsf.xvg -oq ${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/rmsf_bfac.pdb -res
        modify_xvg "${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/rmsf.xvg" "RMSF Per Residue of Chain H at ${temp}K (${sim_time}ns)" "Residue Number" "RMSF (nm)"
        annotate_rmsf "${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/rmsf.xvg" "${REP_DIR}/result_analysis/chains/ChainH/${temp}K/${sim_time}ns/high_rmsf.txt"
        
        # Chain D analysis (Group 11)
        echo "Analyzing Chain D..."
        
        # RMSD for Chain D
        echo "11 
        11" | gmx rms -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/rmsd.xvg -tu ns
        modify_xvg "${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/rmsd.xvg" "RMSD of Chain D at ${temp}K (${sim_time}ns)" "Time (ns)" "RMSD (nm)"
        
        echo "11 
        11" | gmx rms -s ${REP_DIR}/6_em/em.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/rmsd_xtal.xvg -tu ns
        modify_xvg "${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/rmsd_xtal.xvg" "RMSD of Chain D vs Initial Structure at ${temp}K (${sim_time}ns)" "Time (ns)" "RMSD (nm)"
        
        # Gyration radius for Chain D
        echo "11" | gmx gyrate -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/gyrate.xvg

        # RMSF for Chain D (CA atoms only, per residue)
        echo "14" | gmx rmsf -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -o ${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/rmsf.xvg -oq ${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/rmsf_bfac.pdb -res
        modify_xvg "${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/rmsf.xvg" "RMSF Per Residue of Chain D at ${temp}K (${sim_time}ns)" "Residue Number" "RMSF (nm)"
        annotate_rmsf "${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/rmsf.xvg" "${REP_DIR}/result_analysis/chains/ChainD/${temp}K/${sim_time}ns/high_rmsf.txt"
        
        # Interface analysis
        echo "Analyzing chain interface..."
        
        # Calculate minimum distances between chains
        echo "10 
        11" | gmx mindist -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc -n ${REP_DIR}/chains.ndx -od ${REP_DIR}/result_analysis/chains/interface/${temp}K/${sim_time}ns/mindist.xvg -pi
        modify_xvg "${REP_DIR}/result_analysis/chains/interface/${temp}K/${sim_time}ns/mindist.xvg" "Minimum Distance Between Chains at ${temp}K (${sim_time}ns)" "Time (ns)" "Distance (nm)"
        
        # Interface hydrogen bonds
        printf "10\n11\n" | gmx hbond -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_protein.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -n ${REP_DIR}/chains.ndx \
            -tu ns \
            -num ${REP_DIR}/result_analysis/chains/interface/${temp}K/${sim_time}ns/hbnum.xvg \
            -hbn ${REP_DIR}/result_analysis/chains/interface/${temp}K/${sim_time}ns/hbmap.ndx \
            -hbm ${REP_DIR}/result_analysis/chains/interface/${temp}K/${sim_time}ns/hbmat.xpm
        modify_xvg "${REP_DIR}/result_analysis/chains/interface/${temp}K/${sim_time}ns/hbnum.xvg" "Hydrogen Bonds at Interface at ${temp}K (${sim_time}ns)" "Time (ns)" "Number of H-bonds"
        
        # Contact map between chains
        echo "10 
        11" | gmx mdmat -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -n ${REP_DIR}/chains.ndx \
            -frames ${REP_DIR}/result_analysis/chains/interface/${temp}K/${sim_time}ns/contact_frames.xpm \
            -mean ${REP_DIR}/result_analysis/chains/interface/${temp}K/${sim_time}ns/contact_mean.xpm \
            -t 0.5 \
            -nlevels 20
        
        # Principal component analysis (PCA) for the complex
        echo "Performing principal component analysis..."
        
        # Create a covariance matrix
        echo "12 
        12" | gmx covar -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -n ${REP_DIR}/chains.ndx \
            -o ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/eigenval.xvg \
            -v ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/eigenvec.trr
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/eigenval.xvg" "PCA Eigenvalues at ${temp}K (${sim_time}ns)" "Eigenvector index" "Eigenvalue (nm²)"
        
        # Project trajectory on eigenvectors
        echo "12" | gmx anaeig -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -n ${REP_DIR}/chains.ndx \
            -v ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/eigenvec.trr \
            -proj ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/proj_ev.xvg \
            -extr ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/extreme_structures.pdb \
            -first 1 -last 3
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/proj_ev.xvg" "Projection on Principal Components at ${temp}K (${sim_time}ns)" "Time (ps)" "Projection (nm)"
        
        # Calculate free energy landscape
        echo "12" | gmx anaeig -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -n ${REP_DIR}/chains.ndx \
            -v ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/eigenvec.trr \
            -2d ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/2d_proj.xvg \
            -first 1 -last 2
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/2d_proj.xvg" "Free Energy Landscape at ${temp}K (${sim_time}ns)" "PC1 (nm)" "PC2 (nm)"
        
        # Create extreme structures for visualization
        echo "12" | gmx anaeig -s ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns.tpr \
            -f ${REP_DIR}/9_md/${temp}K/${sim_time}ns/md_${sim_time}ns_centered.xtc \
            -n ${REP_DIR}/chains.ndx \
            -v ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/eigenvec.trr \
            -eig ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/eigenval.xvg \
            -comp ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/ev_components.xvg \
            -rmsf ${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/ev_rmsf.xvg \
            -first 1 -last 3
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/ev_components.xvg" "Eigenvector Components at ${temp}K (${sim_time}ns)" "Atom Index" "Component (nm)"
        modify_xvg "${REP_DIR}/result_analysis/complex/${temp}K/${sim_time}ns/ev_rmsf.xvg" "RMSF from Eigenvectors at ${temp}K (${sim_time}ns)" "Atom Index" "RMSF (nm)"
        
        echo "Analysis for ${temp}K ${sim_time}ns completed."
    done
done