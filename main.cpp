#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <mpi.h>

#include "fdtd3d.h"

constexpr int Num_Individual { 24 };

int main( int argc, char** argv ){

    MPI::Init(argc, argv);
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();
    const int assigned_num = Num_Individual / size;
    int name_length = 256;
    char* name = new char[name_length];
    MPI::Get_processor_name( name, name_length );

    std::ofstream ofs;

    if(rank == 0){
      ofs.open("./time_file/elapsed_time_multi_node.dat");
    }

    int* start_idx = new int[size];
    int* end_idx = new int[size];
    for(int myrank = 0; myrank < size; myrank++){
        start_idx[myrank] = myrank * assigned_num;
        end_idx[myrank] = (myrank + 1) * assigned_num;
    }

    double time_0(0.0), time_1(0.0);

    /* Set Perturbation Infomation */
    perturbation P_info;
    P_info.set_alpha( 10.0 );
    P_info.set_center( 74, 25, Nphi/2 );
    P_info.set_sigma( 2.0e3, 30.0e3 );

    /* Set Y(Year)M(Month)D(Date) */
    date ymd;
    ymd.set_ymd(2016, 3, 1);
    ymd.set_h( 9.0 );

    /* Set geocoordinate points */
    geocoordinate lla_info;
    lla_info.set_point( 32.0, 135.0, (Alt_lower_ionosphere/1.0e3) );

    /* Observation Points on propagation path */
    int Num_obs = ( Nphi - 2*L ) - k_s;
    geocoordinate* obs_p = new geocoordinate[Num_obs + 1];
    for( int k = 0; k <= Num_obs; k++ ){
        obs_p[k].set_obs( 0, 50, k + k_s );
    }

    /* Magnitude */
    double *Magnitude = new double[Num_obs + 1];

    std::cout << name << " rank : " << rank << std::endl;

    MPI::COMM_WORLD.Barrier();
    time_0 = MPI::Wtime();

    for( int i = start_idx[rank]; i < end_idx[rank]; i++ ){

        fdtd_calc(P_info, ymd, lla_info, Num_obs, obs_p, Magnitude, rank, 0);
        MPI::COMM_WORLD.Barrier();
    }

    MPI::COMM_WORLD.Barrier();

    time_1 = MPI::Wtime();

    if( rank == 0 ) ofs << size << " " << time_1 - time_0 << std::endl;

    MPI::Finalize();

    std::cout << " Erapsed times : " << time_1 - time_0 << std::endl;

    return 0;
}











