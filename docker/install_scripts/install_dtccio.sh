git clone https://gitlab.com/dtcc-platform/dtcc-io

cd dtcc-io
git submodule update --init --recursive
cd dtcc-model
git checkout develop
cd ..

mkdir build && cd build
cmake ..
make -j4 && make install
cd ..
python3 -m pip install .
