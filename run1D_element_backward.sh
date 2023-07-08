#!/bin/bash

echo "Name the simulation directory:"
read filename

if [ -d "./simu1D_element_backward/$filename" ] 
then
    echo "Directory simu1D_element_backward/$filename exists." 
    exit
fi

# create simualtion dir
mkdir -p simu1D_element_backward/$filename/input

# copy input files
cp -r input1D_element_backward/* ./simu1D_element_backward/$filename/input/

# copy binary
cp axisem3d ./simu1D_element_backward/$filename/

# run
cd simu1D_element_backward/$filename
mpirun -np 16 axisem3d
cd ../..
