# create simualtion dir
mkdir -p simu3D_CMB/input

# copy input files
cp -r input3D_CMB/* ./simu3D_CMB/input/

# copy binary
cp axisem3d ./simu3D_CMB/

# run
cd simu3D_CMB
mpirun -np 16 axisem3d
cd ..
