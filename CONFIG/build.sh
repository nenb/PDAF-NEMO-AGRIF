#! /bin/sh -

module ()
{
    eval `/usr/local/Modules/bin/modulecmd bash $*`
}

module unload cray-mpich
module unload craype-network-ofi
module load craype-network-ucx
module load cray-mpich-ucx
module load libfabric
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
module load gcc

touch prep_build

{
    module list
} &> prep_build


{
    echo y | ./makenemo -n NEMO_Build -r ORCA2_LIM -m X86_ARCHER2-Cray -j 0 clean_config
} >> prep_build 2>&1

{
    ./makenemo -n NEMO_Build -r ORCA2_LIM -m X86_ARCHER2-Cray -j 0
} >> prep_build 2>&1

{
    cp cpp_NEMO_Build.fcm NEMO_Build/
} >> prep_build 2>&1

{
    cp -r MY_SRC/ NEMO_Build/
} >> prep_build 2>&1

{
    ./makenemo -n NEMO_Build -m X86_ARCHER2-Cray -j 8
}
