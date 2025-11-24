for i in 2 4; do
    export DENS_FILE=$(printf "/mnt/nvme1/kboone/Sims/kN%01d_N.npy" "$i")
    export FILE_ROOT=$(printf "/mnt/nvme1/kboone/Sims/BigFiles/Sim_Matter%01d_N/Matter" "$i")
    
    export IN_FILE=${FILE_ROOT}.bigfile
    export NGRID=512
    export NGRIDCIC=512
    export BOX_SIZE=2000.
    
    export DK_IN_KF=1.
    export KMIN_IN_KF=0.5
    
    export OUT_JSON_FILE=${FILE_ROOT}_ps.json
    export OUT_DAT_FILE=${FILE_ROOT}_ps.dat

    python Dens_Array_To_Mesh.py $DENS_FILE $IN_FILE $BOX_SIZE
    
    mpirun -np 8 python ../../../../bskit/scripts/measure/measure_ps_from_bigfile.py $IN_FILE $OUT_JSON_FILE $OUT_DAT_FILE $NGRID $NGRIDCIC $BOX_SIZE $DK_IN_KF $KMIN_IN_KF

    rm -r "$IN_FILE"

done