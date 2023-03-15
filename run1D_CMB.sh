#!/bin/bash

echo "Name the simulation directory:"
read filename

if [ -d "./simu1D_CMB/$filename" ] 
then
    echo "Directory simu1D_CMB/$filename exists." 
    exit
fi

# create simualtion dir
mkdir -p simu1D_CMB/$filename/input

# copy input files
cp -r input1D_CMB/* ./simu1D_CMB/$filename/input/

# copy binary
cp axisem3d ./simu1D_CMB/$filename/

# run
cd simu1D_CMB/$filename
mpirun -np 8 axisem3d
cd ../..
