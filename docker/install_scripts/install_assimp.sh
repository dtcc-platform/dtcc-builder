#!/bin/bash
git clone https://github.com/assimp/assimp && cd assimp && git checkout 1427e67b54906419e9f83cc8625e2207fbb0fcd5 
mkdir build && cd build && cmake -D CMAKE_CXX_FLAGS="-Wno-error=tautological-compare -Wno-error=class-memaccess" .. && make all -j 4 && sudo make install 
#sudo ln -s /assimp/build/bin/libassimp.so.5 /usr/lib/x86_64-linux-gnu/

