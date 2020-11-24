#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<mpi.h>
#include "../../../modules/task_2/nadir_method_gauss/methodGauss.h"

TEST(GAUSS_PARALLEL_MPI, 2x3)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int row = 2, col = 3;
    double* array{ new double[row * col]{ 1,2,3,4,5,6} };
    double* sub_solution1 = new double[row];

    methodGaussParallel(array, sub_solution1);
    if (rank == 0) {double* sub_solution2 = new double[row];
        methodGauss(array, sub_solution2);
        for (int i = 0; i < row; i++) {ASSERT_EQ(sub_solution1[i], sub_solution2[i]);}
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
    double* array{ new double[row * col]{ 1,2,3,4,5,6,7,8,9,10,11,12 } };
    double* sub_solution1 = new double[row];

    methodGaussParallel(array, sub_solution1);
    if (rank == 0) {double* sub_solution2 = new double[row];
        methodGauss(array, sub_solution2);
        for (int i = 0; i < row; i++) { ASSERT_EQ(sub_solution1[i], sub_solution2[i]); }
        delete[] sub_solution2;
    }
    delete[] array;
    delete[] sub_solution1;
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
