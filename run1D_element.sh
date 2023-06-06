#!/bin/bash

echo "Name the simulation directory:"
read filename

if [ -d "./simu1D_element/$filename" ] 
then
    echo "Directory simu1D_element/$filename exists." 
    exit
fi

# create simualtion dir
mkdir -p simu1D_element/$filename/input

# copy input files
cp -r input1D_element/* ./simu1D_element/$filename/input/

# copy binary
cp axisem3d ./simu1D_element/$filename/

# run
cd simu1D_element/$filename
mpirun -np 16 axisem3d
cd ../..
