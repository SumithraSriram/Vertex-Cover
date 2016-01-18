#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "utils.h"

using namespace std;
using Id = unsigned long;
using Edge = pair<Id, Id>;

class Vertex{
public:
	unordered_set<Id> neighs;
};

namespace std {
	template <>
	class hash<Edge>{
	public:
		size_t operator()( const Edge &name ) const;
	};
};

class Graph{
public:
	string filename;
	vector<Vertex> vertices;
	unordered_map<Edge, bool> edges;

	Graph( string infile );

	void check_coverage( vector<Id> &VC );
};

Edge edge( Id e1, Id e2 );

#endif
