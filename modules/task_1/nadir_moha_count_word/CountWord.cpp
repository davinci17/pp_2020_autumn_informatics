#include "CountWord.h"
#include <iostream>
#include<random>
#include<cctype>
#include<string>
#include<mpi.h>

using namespace std;
string GenerateRandom() 
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> distr1(1, 50);
    string word;
    for (int i = 0; i < distr1(gen); i++)
    {
         word = "";
        std::uniform_int_distribution<unsigned int> distr(1, 25);

        for (int j = 0; j < distr(gen); j++) { word += 'a' + distr(gen); }

        int space = (distr(gen) % 2) + 1;
        for (int k = 0; k < space; k++) { word += " "; }
    }
    return word ;
}
int CountWordOneP(std::string Phrases)
{
    int count = 0;
    for (unsigned int i = 0; i < Phrases.size(); i++)
    {
        if (isspace(Phrases[i]) && Phrases[i] != '\0')
        {
            if (isspace(Phrases[i + 1])) { continue; }
            else count++;
        }
        else continue;
    }
    count++;
    return count;
}
int CountWordMultiP(string Phrases)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
        int MiniTask, RestMiniTask;
        MiniTask = Phrases.size() / size;
        RestMiniTask = Phrases.size() % size;

        if (rank == 0) 
        {
            for (int NumProc = 1; NumProc < size; NumProc++) 
            {
                int Buf = RestMiniTask + MiniTask * NumProc;
                MPI_Send(&Buf, 1, MPI_INT, NumProc, 0, MPI_COMM_WORLD);
            }
        }
        else { cout << "from another process"; }
        MPI_Bcast(&MiniTask, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int Count;
        int TotalCount = 0;
        std::string Unity = "";
        int Taken;
       

        if (rank == 0) 
        {
            Unity = Phrases.substr(0, RestMiniTask + MiniTask);
        }
        else {
            MPI_Status status;
            MPI_Recv(&Taken, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            for (int Pro = Taken + 1; Pro <= Taken + MiniTask; Pro++) {
                Unity += Phrases[Pro];
            }
        } 
        TotalCount = CountWordOneP(Unity);
        MPI_Reduce(&Count, &TotalCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        return TotalCount++;
    
}



