#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<mpi.h>
#include "../../../modules/task_2/nadir_method_gauss/methodGauss.h"


TEST(GAUSS_PARALLEL_MPI, Random)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::vector<double>> vec = GenerateRandom();

    std::vector<std::vector<double>> sub_vec1 = methodGauss(vec);
    if (rank == 0) {
        std::vector<std::vector<double>> sub_vec2 = methodGaussParallel(vec);
        ASSERT_EQ(sub_vec1, sub_vec2);
    }
}
TEST(GAUSS_PARALLEL_MPI, VEC_3x3)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::vector<double>> vec = { {1,2,3} ,{4,5,6}, {7,8,9} };

    std::vector<std::vector<double>> sub_vec1 = methodGauss(vec);
    if (rank == 0) {
        std::vector<std::vector<double>> sub_vec2 = methodGaussParallel(vec);
        ASSERT_EQ(sub_vec1, sub_vec2);
    }
}
TEST(GAUSS_PARALLEL_MPI, VEC_4x4)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::vector<double>> vec = { {1,2,3,4} ,{5,6,7,8}, {9,10,11,12} };

    std::vector<std::vector<double>> sub_vec1 = methodGauss(vec);
    if (rank == 0) {
        std::vector<std::vector<double>> sub_vec2 = methodGaussParallel(vec);
        ASSERT_EQ(sub_vec1, sub_vec2);
    }
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
