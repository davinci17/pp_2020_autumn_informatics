// Copyright 2020 Nadir Mohammed
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<mpi.h>
#include "./methodGauss.h"


TEST(GAUSS_PARALLEL_MPI, 2x3)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int row = 2, col = 3;

    double* array1{ new double[row * col]{ 1,2,3,4,5,6 } };
    double* sub_solution1=methodGaussParallel(array1, row,col);



    //**************************************
    std::vector<std::vector<double>>array2{ {1,2,3} ,{4,5,6} };
    std::vector<std::vector<double>>array_sol = methodGauss(array2, 2);
    double* sub_solution2 = new double[row * col];
    convertVector(array_sol, sub_solution2);
    for (int i = 0; i < row * col; i++) 
    {
        if (i > 0) { sub_solution2[i] /= 4; }
    }

    if (rank == 0) {
        
        for (int i = 0; i < row*col; i++)
        { 
            ASSERT_NEAR(sub_solution1[i], sub_solution2[i], 0.00001); 
        }
       
       
    }
    delete[] array1;
    delete[] sub_solution1;
    delete[] sub_solution2;
}
TEST(GAUSS_PARALLEL_MPI, 3x4)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int row = 3, col = 4;

    double* array1{ new double[row * col]{ 1,2,3,4,5,6,7,8,9,10,11,12 } };
    double* sub_solution1 = methodGaussParallel(array1, row, col);



    //**************************************
    std::vector<std::vector<double>>array2{ {1,2,3,4} ,{5,6,7,8},{9,10,11,12} };
    std::vector<std::vector<double>>array_sol = methodGauss(array2, 3);
    double* sub_solution2 = new double[row * col];
    convertVector(array_sol, sub_solution2);
    for (int i = 0; i < row * col; i++)
    {
        if (i > 0) { sub_solution2[i] /= 5; }
    }

    if (rank == 0) {

        for (int i = 0; i < row * col; i++)
        {
            ASSERT_NEAR(sub_solution1[i], sub_solution2[i], 0.00001);
        }
    }
    delete[] array1;
    delete[] sub_solution1;
    delete[] sub_solution2;
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
