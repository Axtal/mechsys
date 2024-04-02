/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_TREE_H
#define MECHSYS_TREE_H

// STL
#include <stdlib.h> // for calloc
#include <cstdio>

// iGraph
extern "C"
{
    #include <igraph.h>
}

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>


namespace Util
{

class Tree
{
public:
    // Constructors
    Tree (Array<int> const & NodesPairs, bool Directed=false) : _directed(Directed) { _init_graph(NodesPairs,Directed); _init_path(); } ///< Size(NodesPairs)/2 == number of edges

    // Destructor
    ~Tree ();

    // Set methods
    void Reset   (Array<int> const & NodesPairs, bool Directed=false); ///< Re-initialize the tree with a new set of edges
    void DelEdge (int LeftNodeID, int RightNodeID);                    ///< Remove an edge from current tree

    // Get methods
    size_t nEdges    () const { return igraph_ecount (&_graph); }             ///< Return the number of edges (size) of this tree (graph)
    size_t GetEdge   (int LeftNodeID, int RightNodeID) const;                 ///< Return the edge id corresponding to nodes LeftNodeID and RightNodeID
    void   GetNodes  (size_t i, int & LeftNodeID, int & RightNodeID) const;   ///< Return the left and right nodes of an edge (i)
    void   ShortPath (size_t FromNodeID, size_t ToNodeID, Array<int> & Path); ///< Find the shortest path from FromNodeID to ToNodeID

    // Methods
    void GetClusters (Array< Array<int> > & Clusters, bool Strong=true) const;
    void WriteDOT    (char const * FileKey) const;

    // Operators
    void operator= (Tree const & Other); ///< Assignment operator

private:
    // Data
    bool                _directed;
    igraph_vector_t     _edges; ///< All edges pairs
    igraph_t            _graph; ///< The graph
    igraph_vs_t         _anode; ///< A node
    igraph_es_t         _epair; ///< An edge pair
    igraph_vector_ptr_t _path;  ///< Array with path
    igraph_vector_ptr_t _segs;  ///< Array with edges as segments

    // Methods
    void _init_graph (Array<int> const & NodesPairs, bool Directed);
    void _init_path  ();

}; // class Tree


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Tree::~Tree()
{
    for (int i=0; i<igraph_vector_ptr_size(&_path); ++i)
    {
        igraph_vector_destroy (static_cast<igraph_vector_t*>(VECTOR(_path)[i]));
        free                  (VECTOR(_path)[i]);
    }
    igraph_vs_destroy         (&_anode);
    igraph_es_destroy         (&_epair);
    igraph_vector_ptr_destroy (&_path);
    igraph_vector_ptr_destroy (&_segs);
    igraph_vector_destroy     (&_edges);
    igraph_destroy            (&_graph);
}

inline void Tree::Reset (Array<int> const & NodesPairs, bool Directed)
{
    igraph_vector_destroy (&_edges);
    igraph_destroy        (&_graph);
    _init_graph           (NodesPairs, Directed);
}

inline void Tree::DelEdge (int LeftNodeID, int RightNodeID)
{
    igraph_es_pairs_small (&_epair, IGRAPH_DIRECTED, LeftNodeID, RightNodeID, /*flag*/-1);
    igraph_delete_edges   (&_graph, _epair);
}

inline void Tree::GetNodes (size_t i, int & LeftNodeID, int & RightNodeID) const
{
    igraph_integer_t left;
    igraph_integer_t right;
    igraph_edge (&_graph, i, &left, &right);
    LeftNodeID  = left;
    RightNodeID = right;
}

inline size_t Tree::GetEdge (int LeftNodeID, int RightNodeID) const
{
    igraph_integer_t eid;
    igraph_get_eid (&_graph, &eid, LeftNodeID, RightNodeID, /*directed*/1, false);
    return eid;
}

inline void Tree::ShortPath (size_t FromNodeID, size_t ToNodeID, Array<int> & Path)
{
    igraph_vs_vector_small    (&_anode, ToNodeID, /*flag*/-1);
    igraph_get_shortest_paths (&_graph, &_path, &_segs, FromNodeID, _anode, IGRAPH_ALL, NULL, NULL);
    size_t npaths = igraph_vector_ptr_size(&_path);
    if (npaths>1) throw new Fatal("Tree::ShortPath: There are (%d) more than one shortest path",npaths);
    size_t path_sz = igraph_vector_size(static_cast<igraph_vector_t*>(VECTOR(_path)[0]));
    Path.Resize(path_sz);
    for (size_t i=0; i<path_sz; ++i)
        Path[i] = VECTOR(*static_cast<igraph_vector_t*>(VECTOR(_path)[0]))[i];
}

inline void Tree::GetClusters (Array< Array<int> > & Clusters, bool Strong) const
{
    igraph_vector_t membership, csize;
    igraph_vector_init (&membership, 0);
    igraph_vector_init (&csize,      0);
    igraph_integer_t no;
    if (_directed && Strong) igraph_clusters (&_graph, &membership, &csize, &no, IGRAPH_STRONG);
    else                     igraph_clusters (&_graph, &membership, &csize, &no, IGRAPH_WEAK);

    Clusters.Resize (no);
    int sz = igraph_vector_size(&membership);
    for (int i=0; i<sz; ++i)
    {
        int cluster_id = VECTOR(membership)[i];
        Clusters[cluster_id].Push (i);
    }

    //std::cout << "no = " << no << std::endl;
    //std::cout << "sz = " << sz << std::endl;
}

inline void Tree::WriteDOT (char const * FileKey) const
{
    String fkey(FileKey);
    fkey.append(".dot");
    FILE * file = fopen(fkey.CStr(),"w");
    if (file!=NULL)
    {
        igraph_write_graph_dot (&_graph, file);
        //igraph_write_graph_graphml (&_graph, file);
        fclose(file);
    }
    else throw new Fatal("Tree::WriteDOT: Could not open file <%s>",fkey.CStr());
}

// Operators

inline void Tree::operator= (Tree const & Other)
{
    Array<int> nodes_pairs;
    nodes_pairs.Resize(Other.nEdges()*2);
    size_t k = 0;
    for (size_t i=0; i<Other.nEdges(); ++i)
    {
        int lef, rig;
        Other.GetNodes (i, lef, rig);
        nodes_pairs[k  ] = lef;
        nodes_pairs[k+1] = rig;
        k += 2;
    }
    Reset (nodes_pairs);
}


/* private */

inline void Tree::_init_graph (Array<int> const & NodesPairs, bool Directed)
{
    igraph_vector_init (&_edges, NodesPairs.Size());
    for (size_t i=0; i<NodesPairs.Size(); ++i)
        VECTOR(_edges)[i] = NodesPairs[i];
    igraph_create (&_graph, &_edges, /*nVerts=auto*/0, /*Directed*/Directed);
}

inline void Tree::_init_path ()
{
    igraph_vector_ptr_init (&_path, 1);
    for (int i=0; i<igraph_vector_ptr_size(&_path); ++i)
    {
        VECTOR(_path)[i] = calloc (1, sizeof(igraph_vector_t));
        igraph_vector_init (static_cast<igraph_vector_t*>(VECTOR(_path)[i]), 0);
    }
}


/* Output */

std::ostream & operator<< (std::ostream & os, Util::Tree const & T)
{
    for (size_t i=0; i<T.nEdges(); ++i)
    {
        int lef, rig;
        T.GetNodes (i, lef, rig);
        os << i << ":(" << lef << ", " << rig << ") ";
    }
    return os;
}


}; // namespace Util

#endif // MECHSYS_TREE_H
