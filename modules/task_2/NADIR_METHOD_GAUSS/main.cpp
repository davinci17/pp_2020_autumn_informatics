// Copyright 2020 Nadir Mohammed
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<mpi.h>
#include "./Methodgauss.h"



 TEST(GAUSS_PARALLEL_MPI, 2x3)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int row = 2, col = 3;
    double* array= new double[row * col];
    array[0] = 1; array[1] = 2;
    array[2] = 3; array[3] = 4;
    array[4] = 5; array[5] = 6;
    double* sub_solution1 = new double[row];

    methodGaussParallel(array, sub_solution1, row,col);
    if (rank == 0) {double* sub_solution2 = new double[row];
        methodGauss(array, sub_solution2, row);
        for (int i = 0; i < row; i++) { 
         
           ASSERT_NEAR(sub_solution1[i], sub_solution2[i], 0.00001); 
        }
       
        delete[] sub_solution2;
    }
    delete[] array;
    delete[] sub_solution1;
}
TEST(GAUSS_PARALLEL_MPI, 3x4)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int row = 3, col = 4;
    double* array = new double[row * col];
    array[0] = 1; array[1] = 2; array[2] = 3; 
    array[3] = 4; array[4] = 5; array[5] = 6;
    array[6] = 7; array[7] = 8; array[8] = 9;
    array[9] = 10; array[10] = 11; array[11] = 12;
    double* sub_solution1 = new double[row];

    methodGaussParallel(array, sub_solution1,row ,col);
    if (rank == 0) {double* sub_solution2 = new double[row];
        methodGauss(array, sub_solution2, row);
        for (int i = 0; i < row; i++) {
          

                ASSERT_NEAR(sub_solution1[i], sub_solution2[i], 0.00001);
            }
        delete[] sub_solution2;
    }
    delete[] array;
    delete[] sub_solution1;
}
