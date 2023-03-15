git clone https://gitlab.com/dtcc-platform/dtcc-io

cd dtcc-io
git checkout develop
git submodule update --init --recursive
cd dtcc-model
git checkout develop
git pull
cd ..

mkdir build && cd build
cmake ..
make -j4 && make install
cd ../..
rm -rf dtcc-io


