#!/bin/bash
git clone -n https://github.com/assimp/assimp
cd assimp 
git checkout tags/v5.2.3 -b v5.2.3
mkdir build 
cd build 
cmake .. -DASSIMP_BUILD_ZLIB=ON 
make all -j 4 
sudo make install 
# sudo ln -s /assimp/build/bin/libassimp.so.5 /usr/lib/x86_64-linux-gnu/
cd ../..
rm -rf assimp

