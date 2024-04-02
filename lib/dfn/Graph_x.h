
#ifndef GRAPH_X_H
#define GRAPH_X_H

#include <iostream>
#include <ctime>
#include <cmath>
#include <queue>
#include <mechsys/linalg/matvec.h>

class Graph_i
{
public:
    size_t NumOfFractures; ///< variable with the number of fractures
    //Array<Array<size_t>> G;
    Array<Array<size_t>> Adj; ///< adjacent list
    Array<bool> Visited;      ///< history: if this fracture has been visited or not

public:
    Graph_i(const size_t NumOfFractures_1, const Array<size_t> Connections);
    void DFS(size_t V, Array<size_t> &B, size_t &Tag_r);
    void CreateGraph_i(Array<Array<size_t>> &S);
};

inline Graph_i::Graph_i(const size_t NumOfFractures_1, const Array<size_t> Connections)
{
    //std::cout << "debug_0.01\n";
    NumOfFractures = NumOfFractures_1;

    /*G.resize(NumOfFractures);

    for (size_t i = 0; i < NumOfFractures; ++i)
        G[i].resize(NumOfFractures);

    for (size_t i = 0; i < NumOfFractures; ++i)
    {
        for (size_t j = 0; j < NumOfFractures; ++j)
            G[i][j] = 0;
    }*/
    Adj.resize(NumOfFractures_1);
    //std::cout << "debug_0.02\n";
    for (size_t i = 0; i < NumOfFractures; ++i)
    {
        Adj[i].resize(1);
        Adj[i][0] = i;
    }

    //std::cout << "debug_0.03\n";

    for (size_t i = 0; i < Connections.Size() / 2; ++i)
    {
        //G[Connections[2 * i]][Connections[2 * i + 1]] = 1;
        //G[Connections[2 * i + 1]][Connections[2 * i]] = 1;
        Adj[Connections[2 * i]].Push(Connections[2 * i + 1]);
        Adj[Connections[2 * i + 1]].Push(Connections[2 * i]);
    }

    Visited.resize(NumOfFractures);
    for (size_t i = 0; i < NumOfFractures; ++i)
    {
        Visited[i] = false;
    }
}

inline void Graph_i::DFS(size_t V, Array<size_t> &B /*one array to record one cluster*/, size_t &Tag_r)
{
    Visited[V] = true;
    B[Tag_r] = V;
    //B.Push(V);

    /*for (size_t i = 0; i < NumOfFractures; ++i)
    {
        if (Visited[i] == false && G[V][i] == 1)
        {
            Tag_r++;
            DFS(i, B, Tag_r);
        }
    }*/

    for (size_t i = 0; i < Adj[V].Size(); i++)
    {
        if (!Visited[Adj[V][i]])
        {
            Tag_r++;
            DFS(Adj[V][i], B, Tag_r);
        }
    }
}

inline void Graph_i::CreateGraph_i(Array<Array<size_t>> &S)
{
    //std::cout << "graph_x\n";
    if (S.Size() != 0)
    {
        std::cout << "Error! ListOfClusters should be empty array!\n";
        exit(0);
    }
    for (size_t i = 0; i < NumOfFractures; ++i)
    {
        if (!Visited[i])
        {
            size_t Tag_r = 0;
            //Array<size_t> A;
            Array<size_t> A(NumOfFractures);
            //std::cout << "debug_0.5\n";
            DFS(i, A, Tag_r); // so, at least, the first element of A will form a cluster
            //std::cout << "debug_1\n";
            //need a function to determine true size of one cluster, i.e., provided A = [0,1,2,3,0,0,0], the size of this cluster is four

            Array<size_t> new_A(Tag_r + 1);

            S.Push(new_A);
            //std::cout << "debug_2\n";
            for (size_t i = 0; i < Tag_r + 1; ++i)
                S[S.Size() - 1][i] = A[i];
            //std::cout << "debug_3\n";
            //S.Push(A);
            //S.Push(new_A); //
        }
    }
}
#endif
