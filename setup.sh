source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
export LD_LIBRARY_PATH=/home/edesanti/FCCk4hep/k4geo/build/lib:$LD_LIBRARY_PATH
pwd=$(pwd)
cd k4geo/InstallArea
source bin/thisk4geo.sh
cd ..
k4_local_repo InstallArea/
cd ${pwd}