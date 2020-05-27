#! /bin/sh -

module ()
{
    eval `/opt/modules/3.2.10.6/bin/modulecmd bash $*`
}

module unload PrgEnv-cray
module load PrgEnv-intel
module unload ncview
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel

touch prep_build

{
    module list
} &> prep_build


{
    echo y | ./makenemo -n NEMO_Build -r ORCA2_LIM -m XC_ARCHER_INTEL -j0 clean_config
} >> prep_build 2>&1

{
    ./makenemo -n NEMO_Build -r ORCA2_LIM -m XC_ARCHER_INTEL -j0
} >> prep_build 2>&1

{
    cp cpp_NEMO_Build.fcm NEMO_Build/
} >> prep_build 2>&1

{
    cp -r MY_SRC/ NEMO_Build/
} >> prep_build 2>&1

{
    ./makenemo -n NEMO_Build -m XC_ARCHER_INTEL -j8
}

module unload cray-hdf5-parallel
module unload cray-netcdf-hdf5parallel
module load ncview
