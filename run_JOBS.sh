#!/bin/bash

for dir in JOBS/*
do 
    filename=${dir##*/}
    echo "WORKING ON ${dir##*/}"

    if [ -d ".simu_JOBS/$filename" ]
    then
        echo "Directory simu_JOBS/$filename exists."
        exist
    fi

    mkdir -p simu_JOBS/$filename/input
    cp -r JOBS/$filename/* ./simu_JOBS/$filename/input/
    cp ./axisem3d ./simu_JOBS/$filename/
    cd simu_JOBS/$filename
    mpirun -np 16 ./axisem3d > LOGS_$filename.txt
    cd ../..
done