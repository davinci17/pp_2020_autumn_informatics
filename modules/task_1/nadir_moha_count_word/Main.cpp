#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <stdio.h>
#include<string>
#include<mpi.h>
#include<ctime>
#include <random>
#include "CountWord.h"



TEST(Count_Word, space) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int result = CountWordMultiP("");
    if (rank == 0) {
        ASSERT_EQ(result, 0);
    }
}
TEST(Count_Word, First_Generate_False) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int result = CountWordMultiP("Take It Easy");
    if (rank == 0) {
        ASSERT_EQ(result, 0);
    }
   
}
TEST(Count_Word, First_Generate_True) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int result = CountWordMultiP("I Try To Do It");
    if (rank == 0) {
        ASSERT_EQ(result, 5);
    }
}
TEST(Count_Word, Second_Generate_SeveralSpace) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int result = CountWordMultiP("The End Of  The   Game");
    if (rank == 0) {
        ASSERT_EQ(result, 5);
    }
}

TEST(Count_Words, Random_String) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int res = CountWordMultiP(GenerateRandom());
    if (rank == 0) {
        ASSERT_EQ(res, 1);
    }
}

int main
(int argc, char** argv) {
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