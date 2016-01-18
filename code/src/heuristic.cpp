#include "heuristic.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <algorithm>

vector<Id> greedybad( Graph &G ) {
	// Init
	vector<Id> VC;
	vector<Vertex> vertices = G.vertices;

	unordered_set<Edge> uncovered;
	for ( auto &p : G.edges )
		uncovered.insert( p.first );
	while ( !uncovered.empty() )
	{
		// Find vertex vertices minimum degree
		Id v = vertices.size();
		for ( Id i = 0; i < vertices.size(); ++i )
			if ( vertices[i].neighs.size() > 0
				 && ( v == vertices.size() || vertices[v].neighs.size() > vertices[i].neighs.size() ) )
				v = i;
		if ( v == vertices.size() || vertices[v].neighs.empty() )
			break;

		//Finding v's neighbour with minimum degree
		Id minneigh = *min_element( vertices[v].neighs.begin(), vertices[v].neighs.end(),
									[&]( Id a, Id b ){
			return vertices[a].neighs.size() < vertices[b].neighs.size();
		} );

		Id e[2] = { v, minneigh };
		for ( Id i : e ){
			VC.push_back( i );
			for ( Id j : vertices[i].neighs )
			{
				// Remove edges from queues
				uncovered.erase( edge( i, j ) );
				vertices[j].neighs.erase( i );
			}
			vertices[i].neighs.clear();
		}
	}

	return VC;
}

// The same one that is outlined in the slides of Lecture 19
// TODO: Version where orphan edge only adds one vertex to VC.
vector<Id> heuristic1( Graph &G ){
	// Init
	vector<Id> VC;
	unordered_set<Edge> uncovered;
	for ( auto &p : G.edges )
		uncovered.insert( p.first );
	while ( !uncovered.empty() ){
		Edge e = *uncovered.begin();

		// Add to VC
		uncovered.erase( uncovered.begin() );
		VC.push_back( e.first );
		VC.push_back( e.second );

		// Remove neighbor edges
		for ( Id j : G.vertices[e.first].neighs )
			uncovered.erase( edge( e.first, j ) );
		for ( Id j : G.vertices[e.second].neighs )
			uncovered.erase( edge( e.second, j ) );
	}
	return VC;
}

// An optimization of the original heuristic
// It prioritizes the neighbor of edges that has a vertex with only one uncovered edge left
vector<Id> heuristic2( Graph &G ){
	// Init
	vector<Id> VC;
	unordered_set<Edge> uncovered;
	unordered_set<Edge> priority_uncovered;
	vector<Vertex> vertices = G.vertices;
	for ( auto &p : G.edges ){
		if ( G.vertices[p.first.first].neighs.size() == 1 || G.vertices[p.first.second].neighs.size() == 1 )
			priority_uncovered.insert( p.first );
		else
			uncovered.insert( p.first );
	}
	while ( !uncovered.empty() || !priority_uncovered.empty() ){
		Edge e;

		// Get next uncovered edge
		if ( !priority_uncovered.empty() ){
			e = *priority_uncovered.begin();
			priority_uncovered.erase( priority_uncovered.begin() );
			// Get the edges that touches our priority edge
			// And select the one which has a vertex with max number of uncovered neighbor
			if ( vertices[e.first].neighs.size() == 1 ){
				int maxneigh = *vertices[e.second].neighs.begin();
				for ( int j : vertices[e.second].neighs )
					if ( vertices[j].neighs.size() > vertices[maxneigh].neighs.size() )
						maxneigh = j;
				e = edge( e.second, maxneigh );
			}
			else if ( vertices[e.second].neighs.size() == 1 ){
				int maxneigh = *vertices[e.first].neighs.begin();
				for ( int j : vertices[e.first].neighs )
					if ( vertices[j].neighs.size() > vertices[maxneigh].neighs.size() )
						maxneigh = j;
				e = edge( e.first, maxneigh );
			}
		}
		else{
			e = *uncovered.begin();
			uncovered.erase( uncovered.begin() );
		}

		// Push to vertex cover list
		// If edge is isolated, only one vertex is needed.
		if ( vertices[e.first].neighs.size() > 1 )
			VC.push_back( e.first );
		if ( vertices[e.second].neighs.size() > 1 || vertices[e.first].neighs.size() == 1 )
			VC.push_back( e.second );

		// Cleanup
		for ( Id j : G.vertices[e.first].neighs ){
			Edge neigh = edge( e.first, j );
			// Remove edges from queues
			uncovered.erase( neigh );
			priority_uncovered.erase( neigh );
			// Remove edge from neighbor's edgelist
			vertices[j].neighs.erase( e.first );
			// If neighbor has one uncovered edge left, it becomes priority
			if ( vertices[j].neighs.size() == 1 ){
				Edge p_edge = edge( j, *vertices[j].neighs.begin() );
				if ( uncovered.find( p_edge ) != uncovered.end() ){
					priority_uncovered.insert( p_edge );
					uncovered.erase( p_edge );
				}
			}
		}
		for ( Id j : G.vertices[e.second].neighs ){
			Edge neigh = edge( e.second, j );
			// Remove edges from queues
			uncovered.erase( neigh );
			priority_uncovered.erase( neigh );
			// Remove edge from neighbor's edgelist
			vertices[j].neighs.erase( e.second );
			// If neighbor has one uncovered edge left, it becomes priority
			if ( vertices[j].neighs.size() == 1 ){
				Edge p_edge = edge( j, *vertices[j].neighs.begin() );
				if ( uncovered.find( p_edge ) != uncovered.end() ){
					priority_uncovered.insert( p_edge );
					uncovered.erase( p_edge );
				}
			}
		}
	}
	return VC;
}

// An optimization of the original heuristic
// It prioritizes edges that has a vertex with only one uncovered edge left, but only adds the other vertex
vector<Id> heuristic3( Graph &G ){
	// Init
	vector<Id> VC;
	unordered_set<Edge> uncovered;
	unordered_set<Edge> priority_uncovered;
	vector<Vertex> vertices = G.vertices;
	for ( auto &p : G.edges ){
		if ( G.vertices[p.first.first].neighs.size() == 1 || G.vertices[p.first.second].neighs.size() == 1 )
			priority_uncovered.insert( p.first );
		else
			uncovered.insert( p.first );
	}
	while ( !uncovered.empty() || !priority_uncovered.empty() ){
		Edge e;

		// Get next uncovered edge
		if ( !priority_uncovered.empty() ){
			e = *priority_uncovered.begin();
			priority_uncovered.erase( priority_uncovered.begin() );
		}
		else{
			e = *uncovered.begin();
			uncovered.erase( uncovered.begin() );
		}

		// Push to vertex cover list
		// If edge is isolated, only one vertex is needed.
		if ( vertices[e.first].neighs.size() > 1 )
			VC.push_back( e.first );
		if ( vertices[e.second].neighs.size() > 1 || vertices[e.first].neighs.size() == 1 )
			VC.push_back( e.second );

		// Cleanup
		for ( Id j : G.vertices[e.first].neighs ){
			Edge neigh = edge( e.first, j );
			// Remove edges from queues
			uncovered.erase( neigh );
			priority_uncovered.erase( neigh );
			// Remove edge from neighbor's edgelist
			vertices[j].neighs.erase( e.first );
			// If neighbor has one uncovered edge left, it becomes priority
			if ( vertices[j].neighs.size() == 1 ){
				Edge p_edge = edge( j, *vertices[j].neighs.begin() );
				if ( uncovered.find( p_edge ) != uncovered.end() ){
					priority_uncovered.insert( p_edge );
					uncovered.erase( p_edge );
				}
			}
		}
		for ( Id j : G.vertices[e.second].neighs ){
			Edge neigh = edge( e.second, j );
			// Remove edges from queues
			uncovered.erase( neigh );
			priority_uncovered.erase( neigh );
			// Remove edge from neighbor's edgelist
			vertices[j].neighs.erase( e.second );
			// If neighbor has one uncovered edge left, it becomes priority
			if ( vertices[j].neighs.size() == 1 ){
				Edge p_edge = edge( j, *vertices[j].neighs.begin() );
				if ( uncovered.find( p_edge ) != uncovered.end() ){
					priority_uncovered.insert( p_edge );
					uncovered.erase( p_edge );
				}
			}
		}
	}
	return VC;
}

void heuristic( Graph &G ){
	ostringstream oss;
	oss << "output/" << G.filename.substr( 0, G.filename.size() - 6u ) << "_Approx";

	chrono::time_point<std::chrono::system_clock> start, end;
	start = chrono::system_clock::now();
	vector<Id> VC = heuristic3( G );
	end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	ofstream ofs( oss.str() + ".trace" );
	ofs << elapsed_seconds.count() << ',' << VC.size();
	ofs.close();

	ofs.open( oss.str() + ".sol" );
	ofs << VC.size() << '\n';
	for ( size_t i = 0; i < VC.size(); ++i ){
		if ( i != 0 )
			ofs << ',';
		ofs << ( VC[i] + 1 );
	}
	G.check_coverage( VC );
}

vector<Id> getBestHeuristic( Graph &G ){
	vector<Id> VC = heuristic1( G );
	{
		vector<Id> temp = heuristic2( G );
		if ( temp.size() < VC.size() )
			VC = move( temp );
	}
	{
		vector<Id> temp = heuristic3( G );
		if ( temp.size() < VC.size() )
			VC = move( temp );
	}
	return VC;
}
