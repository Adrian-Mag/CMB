#!/bin/bash

echo "Name the simulation directory:"
read filename

if [ -d "./simu3D/$filename" ] 
then
    echo "Directory simu3D/$filename exists." 
    exit
fi

# create simualtion dir
mkdir -p simu3D/$filename/input

# copy input files
cp -r input3D/* ./simu3D/$filename/input/

# copy binary
cp axisem3d ./simu3D/$filename/

# run
cd simu3D/$filename
mpirun -np 8 axisem3d
cd ../..
