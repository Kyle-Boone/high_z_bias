export OUTFILE_GRID="/mnt/nvme1/kboone/Sims/BigFiles/Sim_Matter/Matter_grid.dat"

export NGRID=2048
export NGRIDCIC=2048

export LBOX=2000.
export STARTI=0
export ENDI=4000

export DK=0.006283185307179587 # 2*kF # 0.025, divide by factor of 20/3 due to larger box
export KMIN=0.0015707963267948967 # 0.5*kF # 0.025
export KMAX=0.20263272615654168 # 64.5*kF # 2.025

for i in {0..24}; do

    export POS_FILE=$(printf "../../Data/Gala_z8_Positions/Gala_ph0%02d.h5" "$i")
    export FILE_ROOT=$(printf "/mnt/nvme1/kboone/Data/BigFiles/Gala_8/Gala_ph0%02d/Gala" "$i")
    export SNAP=${FILE_ROOT}.bigfile

    python Pos_To_Mesh.py $POS_FILE $SNAP $LBOX $NGRID

    # if [ "$i" -eq 0 ]; then
    #     mpirun -np 32 python ../../../../bskit/scripts/measure/measure_bs_slow.py bigfile_grid $SNAP $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE_GRID --k_bin_info $KMIN $KMAX $DK --triangle_type all --meas_type grid_info
        
    # fi
    
    export OUTFILE_UNNORM=${FILE_ROOT}_fine_unnormbs.dat
    export OUTFILE_FINAL=${FILE_ROOT}_fine_comb.dat
    
    mpirun -np 8 python ../../../../bskit/scripts/measure/measure_bs_fast.py bigfile_grid $SNAP $NGRID $NGRIDCIC $LBOX $STARTI $ENDI $OUTFILE_UNNORM --k_bin_info $KMIN $KMAX $DK  --triangle_type all
    
    python ../../../../bskit/scripts/process/process_fast_bs_measurement.py $OUTFILE_GRID $OUTFILE_UNNORM $KMAX $OUTFILE_FINAL

    rm -r "$SNAP"

done
