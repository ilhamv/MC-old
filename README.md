# MC
A Monte Carlo Neutron Transport Simulator

# Build and Run
mkdir build;
cd build;
cmake -Doptimize=ON ..;
make;
./MC examples/infinite_GCR_TRMM;
cd examples/infinite_GCR_TRMM/;
python plot;
