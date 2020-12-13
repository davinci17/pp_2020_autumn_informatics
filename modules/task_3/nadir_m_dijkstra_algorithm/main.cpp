// Copyright 2020 Nadir Mohammed
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<iostream>
#include "./dijkstra_algorithm.h"

TEST(FIRST_TEST, 8X8) {
    int size, rank;
    const int r_col = 9;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int graph[r_col][r_col] = { { 0, 4, 0, 0, 0, 0, 0, 8, 0 },
                        { 4, 0, 8, 0, 0, 0, 0, 11, 0 },
                        { 0, 8, 0, 7, 0, 4, 0, 0, 2 },
                        { 0, 0, 7, 0, 9, 14, 0, 0, 0 },
                        { 0, 0, 0, 9, 0, 10, 0, 0, 0 },
                        { 0, 0, 4, 14, 10, 0, 2, 0, 0 },
                        { 0, 0, 0, 0, 0, 2, 0, 1, 6 },
                        { 8, 11, 0, 0, 0, 0, 1, 0, 7 },
                        { 0, 0, 2, 0, 0, 0, 6, 7, 0 } };
    int* arr_resu = new int[r_col];
    dijkstra(graph, 4, arr_resu);
    //**************************************************************
    int graph1[r_col * r_col] = { 0, 4, 0, 0, 0, 0, 0, 8, 0 ,
                  4, 0, 8, 0, 0, 0, 0, 11, 0 ,
                  0, 8, 0, 7, 0, 4, 0, 0, 2 ,
                  0, 0, 7, 0, 9, 14, 0, 0, 0 ,
                  0, 0, 0, 9, 0, 10, 0, 0, 0 ,
                  0, 0, 4, 14, 10, 0, 2, 0, 0 ,
                  0, 0, 0, 0, 0, 2, 0, 1, 6 ,
                  8, 11, 0, 0, 0, 0, 1, 0, 7 ,
                  0, 0, 2, 0, 0, 0, 6, 7, 0 };
    int* arr1_resu = new int[r_col];
    getParallelDijkstras(graph1, 4, arr1_resu);
        if (rank == 0) {
            for (int i = 0; i < r_col; i++) {
                ASSERT_EQ(arr1_resu[i], arr_resu[i]);
            }
        }
}

TEST(FIRST_TEST, 9X9) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int r_col = 9;
    int graph[r_col][r_col] = { { 0, 4, 0, 0, 0, 0, 0, 8, 0 },
                        { 4, 0, 8, 0, 0, 0, 0, 11, 0 },
                        { 0, 8, 0, 7, 0, 4, 0, 0, 2 },
                        { 0, 0, 7, 0, 9, 14, 0, 0, 0 },
                        { 0, 0, 0, 9, 0, 10, 0, 0, 0 },
                        { 0, 0, 4, 14, 10, 0, 2, 0, 0 },
                        { 0, 0, 0, 0, 0, 2, 0, 1, 6 },
                        { 8, 11, 0, 0, 0, 0, 1, 0, 7 },
                        { 0, 0, 2, 0, 0, 0, 6, 7, 0 } };
    int* arr_resu = new int[r_col];

    dijkstra(graph, 4, arr_resu);
    //**************************************************************
    int graph1[r_col * r_col] = { 0, 4, 0, 0, 0, 0, 0, 8, 0,
                  4, 0, 8, 0, 0, 0, 0, 11, 0,
                  0, 8, 0, 7, 0, 4, 0, 0, 2,
                  0, 0, 7, 0, 9, 14, 0, 0, 0,
                  0, 0, 0, 9, 0, 10, 0, 0, 0,
                  0, 0, 4, 14, 10, 0, 2, 0, 0,
                  0, 0, 0, 0, 0, 2, 0, 1, 6,
                  8, 11, 0, 0, 0, 0, 1, 0, 7,
                  0, 0, 2, 0, 0, 0, 6, 7, 0};
    int* arr1_resu = new int[r_col];
    getParallelDijkstras(graph1, 4, arr1_resu);
    if (rank == 0) {
        for (int i = 0; i < r_col; i++) {
            ASSERT_EQ(arr1_resu[i], arr_resu[i]);
        }
    }
}
int main(int argc, char* argv[]) {
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
