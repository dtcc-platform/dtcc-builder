git clone https://github.com/PointCloudLibrary/pcl.git
cd pcl
git checkout tags/pcl-1.12.1 -b pcl-1.12.1
mkdir build && cd build
cmake -DBUILD_visualization=OFF -DBUILD_stereo=OFF -DBUILD_features=OFF -DBUILD_registration=OFF -DBUILD_io=OFF -DBUILD_tracking=off -DBUILD_ml=OFF -DBUILD_recognition=OFF -DWITH_CUDA=OFF ..
make -j4 && make install