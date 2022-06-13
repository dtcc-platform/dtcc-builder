git clone https://github.com/kthohr/optim
cd optim
git submodule update --init
export EIGEN_INCLUDE_PATH=/usr/include/eigen3
./configure --header-only-version
mkdir /usr/local/include/optim
cp -r header_only_version/ /usr/local/include/optim/
