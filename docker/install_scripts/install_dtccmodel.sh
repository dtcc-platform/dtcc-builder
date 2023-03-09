git clone https://github.com/dtcc-platform/dtcc-model.git
cd dtcc-model
git checkout develop

# make sure we rebuild the protobuf files using the
# same verion of protoc as installed in the docker image
rm -rf src/cpp/protobuf/dtcc.pb.*
mkdir build
cd build
cmake ..
make && make install
