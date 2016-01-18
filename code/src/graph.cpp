#include "graph.h"
#include <algorithm>
#include <cassert>
#include <string>

// Hash function for Edge pair. This only works if vertex id bits fit into half of unsigned long long
namespace std {
	size_t hash<Edge>::operator()( const Edge &name ) const
	{
		unsigned long long e = name.first;
		// Shift it to left half of e. sizeof = num of bytes. 1 byte = 8 bit. Shift by half the amount
		e |= static_cast<unsigned long long>( name.second ) << ( 4 * sizeof( unsigned long long ) );
		return hash<unsigned long long>()( e );
	}
};

Graph::Graph( string infile ){
	char delim = find( infile.begin(), infile.end(), '\\' ) != infile.end() ? '\\' : '/';
	auto path = split( infile, delim );
	filename = path.back();
	ifstream ifs( infile );
	if ( !ifs.good() ){
		cout << "File could not be opened!\n";
		exit(1);
	}

	Id N, M, W;
	ifs >> N >> M >> W;
	// Check if our edge data structure stupports this number of vertices
	// unsigned -1 overflows, so we get max value
	assert( N <= ( 1ull << ( sizeof( Id ) * 4 ) ) );
	ifs.ignore();
	vertices.resize( N );
	for ( Id i = 0; i < N; ++i ){
		string line;
		getline( ifs, line );
		auto neighs = split( line, ' ' );
		for ( string &s : neighs ){
			Id j = stoul( s );
			vertices[i].neighs.insert( --j );
			edges[edge(i, j)] = false;
		}
	}
	cout << "Graph initialized! |V| = " << N << ", |E| = " << M << '\n';
}

void Graph::check_coverage( vector<Id> &VC ){
	cout << "cover size: " << VC.size() << '\n';
	int coveredcount = 0;
	for ( auto &p : edges )
		p.second = false;
	for ( size_t i = 0; i < VC.size(); ++i ){
		for ( Id j : vertices[VC[i]].neighs ){
			auto e = edge( VC[i], j );
			if ( !edges[e] ){
				edges[e] = true;
				++coveredcount;
			}
		}
	}
	cout << "covered: " << coveredcount << '\n';
}

Edge edge( Id e1, Id e2 ){
	return Edge( min( e1, e2 ), max( e1, e2 ) );
}