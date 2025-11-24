for i in {0..24}; do
    export POS_FILE=$(printf "../../Data/7E10_Gala_z8_Positions/Gala_ph0%02d.h5" "$i")
    export FILE_ROOT=$(printf "../../Data/BigFiles/7E10_Gala_8/N_1024/Gala_ph0%02d/Gala" "$i")
    # export FILE_ROOT=../../Data/BigFiles/Gala_8/Gala_ph024/Gala
    
    export IN_FILE=${FILE_ROOT}.bigfile
    export NGRID=1024
    export NGRIDCIC=1024
    export BOX_SIZE=2000.
    
    export DK_IN_KF=1.
    export KMIN_IN_KF=0.5
    
    export OUT_JSON_FILE=${FILE_ROOT}_ps.json
    export OUT_DAT_FILE=${FILE_ROOT}_ps.dat

    python Pos_To_Mesh.py $POS_FILE $IN_FILE $BOX_SIZE $NGRID
    
    mpirun -np 8 python ../../../../bskit/scripts/measure/measure_ps_from_bigfile.py $IN_FILE $OUT_JSON_FILE $OUT_DAT_FILE $NGRID $NGRIDCIC $BOX_SIZE $DK_IN_KF $KMIN_IN_KF

    rm -r "$IN_FILE"

done