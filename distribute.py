import os 
from distutils.dir_util import copy_tree

for dir in os.walk('simu_JOBS'):
    for subdir in dir[1]:
        if subdir.find('_1D_CMB') != -1:
            copy_tree(os.path.abspath('simu_JOBS/' + subdir), os.path.abspath('simu1D_CMB/' + subdir))
        elif subdir.find('_1D') != -1:
            copy_tree(os.path.abspath('simu_JOBS/' + subdir), os.path.abspath('simu1D/' + subdir))
        elif subdir.find('_3D_CMB') != -1:
            copy_tree(os.path.abspath('simu_JOBS/' + subdir), os.path.abspath('simu3D_CMB/' + subdir))
        elif subdir.find('_3D') != -1:
            copy_tree(os.path.abspath('simu_JOBS/' + subdir), os.path.abspath('simu3D/' + subdir))
    break
