#!/bin/bash

echo "Name the simulation directory:"
read filename

if [ -d "./simu3D_CMB/$filename" ] 
then
    echo "Directory simu3D_CMB/$filename exists." 
    exit
fi

# create simualtion dir
mkdir -p simu3D_CMB/$filename/input

# copy input files
cp -r input3D_CMB/* ./simu3D_CMB/$filename/input/

# copy binary
cp axisem3d ./simu3D_CMB/$filename/

# run
cd simu3D_CMB/$filename
mpirun -np 16 axisem3d
cd ../..
