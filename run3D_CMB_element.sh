#!/bin/bash

echo "Name the simulation directory:"
read filename

if [ -d "./simu3D_CMB_element/$filename" ] 
then
    echo "Directory simu3D_CMB_element/$filename exists." 
    exit
fi

# create simualtion dir
mkdir -p simu3D_CMB_element/$filename/input

# copy input files
cp -r input3D_CMB_element/* ./simu3D_CMB_element/$filename/input/

# copy binary
cp axisem3d ./simu3D_CMB_element/$filename/

# run
cd simu3D_CMB_element/$filename
mpirun -np 16 axisem3d > LOG.txt
cd ../..