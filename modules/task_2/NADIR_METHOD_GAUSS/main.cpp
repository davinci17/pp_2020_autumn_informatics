#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<mpi.h>
#include "../../../modules/task_2/nadir_method_gauss/methodGauss.h"


TEST(GAUSS_PARALLEL_MPI, 2x3)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int row = 2;
    int col = 3;
    double* array{ new double[row * col]{ 1,2,3,4,5,6} };
    double* solution = new double[row];
    double* sub_vec2 = new double[row * col];
    sub_vec2 = methodGaussParallel(array, solution);



    std::vector<std::vector<double>>array1{ {1,2,3},{4,5,6} };
    std::vector<std::vector<double>> sub_vec1 = methodGauss(array1);

    ASSERT_EQ(sub_vec1, sub_vec2);

}
TEST(GAUSS_PARALLEL_MPI, 3x4)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int row = 3;
    int col = 4;
    double* array{ new double[row * col]{ 1,2,3,4,5,6,7,8,9,10,11,12 } };
    double* solution = new double[row];
    double* sub_vec2 = new double[row * col];

    sub_vec2 = methodGaussParallel(array, solution);



    std::vector<std::vector<double>>array1{ {1,2,3,4} ,{5,6,7,8},{9,10,11,12} };
    std::vector<std::vector<double>> sub_vec1 = methodGauss(array1);

    ASSERT_EQ(sub_vec1, sub_vec2);

}

TEST(GAUSS_PARALLEL_MPI, 4x5)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int row = 4;
    int col = 5;
    double* array{ new double[row * col]{ 5,8,9,4,7,5,6,7,8,10,9,10,11,12,15,13,14,15,16,20} };
    double* solution = new double[row];
    double* sub_vec2 = new double[row * col];
    sub_vec2 = methodGaussParallel(array, solution);


    std::vector<std::vector<double>>array1{ {5,8,9,4,7} ,{5,6,7,8,10},{9,10,11,12,15},{13,14,15,16,20} };
    std::vector<std::vector<double>> sub_vec1 = methodGauss(array1);

    ASSERT_EQ(sub_vec1, sub_vec2);

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
