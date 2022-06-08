#!/bin/bash
wget https://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz
tar xvzf VTK-7.1.1.tar.gz 
cd VTK-7.1.1/ 
mkdir build
cd build
cmake -DVTK_Group_Rendering=OFF -DVTK_BUILD_ALL_MODULES_FOR_TESTS:BOOL=OFF -DVTK_Group_StandAlone=OFF -DModule_vtkCommonCore:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkIOXML:BOOL=ON .. 
make all -j 4
sudo make install

cd ../..
rm -rf VTK-7.1.1
rm VTK-7.1.1.tar.gz


