#!/bin/bash

echo "Name the simulation directory:"
read filename

if [ -d "./simu1D/$filename" ] 
then
    echo "Directory simu1D/$filename exists." 
    exit
fi

# create simualtion dir
mkdir -p simu1D/$filename/input

# copy input files
cp -r input1D/* ./simu1D/$filename/input/

# copy binary
cp axisem3d ./simu1D/$filename/

# run
cd simu1D/$filename
mpirun -np 8 axisem3d
cd ../..
