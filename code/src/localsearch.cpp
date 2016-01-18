#include "localsearch.h"
#include "heuristic.h"
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <iomanip>
#include <random>
#include <chrono>

using namespace std;

using TimePoint = chrono::time_point<std::chrono::system_clock>;
using SecondsDouble = chrono::duration<double>;

// G					- The graph
// VC					- Current solution
// MinVC				- best solution
// nodeW				- Weight for nodes
// edgeW				- Weight for edges
// taboo_swap			- Taboo edge
// gen					- random generator
// cutoff				- cutoff time in seconds
// start				- starting time point of solver
// ofs					- ofstream to write to trace
class StochasticSolver{
private:
	Graph &G;
	unordered_set<Id> VC, MinVC;
	vector<double> nodeW;
	unordered_map<Edge, double> edgeW;
	Edge taboo_swap;

	mt19937 gen;

	double cutoff;
	TimePoint start;
	ofstream ofs;
	string outfile;

public:
	StochasticSolver( Graph &G_, double cutoff_, int seed ) : G( G_ ), gen( seed ), cutoff( cutoff_ ){
		ostringstream oss;
		oss << "output/" << G.filename.substr( 0, G.filename.size() - 6u ) << "_LS1_" << cutoff << '_' << seed;
		outfile = oss.str();
	}

	pair<Id, Id> vertexPairToExchange( const Edge &e ){
		vector<Id> elist = { e.first, e.second };
		// Adjust weight for vertices in VC and e. More edge -> more weight
		for ( Id u : VC ){
			nodeW[u] = 0.0;
			for ( Id v : G.vertices[u].neighs )
				if ( VC.find( v ) == VC.end() )
					nodeW[u] += edgeW[edge( u, v )];
		}
		for ( Id u : elist ){
			nodeW[u] = 0.0;
			for ( Id v : G.vertices[u].neighs )
				if ( VC.find( v ) == VC.end() )
					nodeW[u] += edgeW[edge( u, v )];
		}

		// Calculate weighted gain for neighbors
		double maxGain = -1.0;
		Id a = 0, b = 0;
		for ( Id i : elist ){
			double bWt = nodeW[i];
			for( Id j : VC ){
				Edge cand = edge( i, j );
				if ( taboo_swap == cand )
					continue;
				double gain = bWt - nodeW[j];
				if ( edgeW.find( cand ) != edgeW.end() )
					gain += edgeW[cand];
				if ( gain > maxGain ){
					maxGain = gain;
					a = i;
					b = j;
				}
			}
		}
		return {a, b};
	}

	void ILS(){
		unordered_set<Edge> UncoveredEdges;

		TimePoint end = chrono::system_clock::now();
		for ( SecondsDouble elapsed_seconds = end - start; elapsed_seconds.count() <= cutoff;
			end = chrono::system_clock::now(), elapsed_seconds = end - start ){
			// While there are no uncovered edges
			while ( UncoveredEdges.empty() ){
				if ( VC.size() < MinVC.size() ){
					MinVC = VC;
					end = chrono::system_clock::now();
					elapsed_seconds = end - start;
					ofs << elapsed_seconds.count() << ',' << ( MinVC.size() ) << '\n';
				}
				// Remove a random vertex
				auto it = VC.begin();
				uniform_int_distribution<size_t> idis( 0u, VC.size() - 1u );
				size_t adv = idis( gen );
				advance( it, adv );
				// Add to uncovered edges if needed
				for ( Id i : G.vertices[*it].neighs )
					if ( VC.find( i ) == VC.end() )
						UncoveredEdges.insert( edge(*it, i) );
				// Reset these if we finally uncovered an edge
				if ( !UncoveredEdges.empty() ){
					taboo_swap = { 0, 0 };
					for ( auto &p : edgeW )
						p.second = 0.05;
				}
				VC.erase( it );
			}
			// Get a random uncovered edge
			auto it = UncoveredEdges.begin();
			uniform_int_distribution<size_t> idis( 0u, UncoveredEdges.size() - 1u );
			size_t adv = idis( gen );
			advance( it, adv );
			auto p = vertexPairToExchange( *it );
			// Erase second
			VC.erase( p.second );
			for ( Id i : G.vertices[p.second].neighs )
				if ( VC.find( i ) == VC.end() )
					UncoveredEdges.insert( edge( p.second, i ) );
			// Insert first
			VC.insert( p.first );
			for ( Id i : G.vertices[p.first].neighs )
				UncoveredEdges.erase( edge( p.first, i ) );

			// Add this edge to taboo
			taboo_swap = edge( p.first, p.second );

			// Increment weights
			for ( const Edge &e : UncoveredEdges )
				++edgeW[e];			
		}
	}

	vector<Id> solve(){
		vector<Id> res;
		
		// Solve
		start = chrono::system_clock::now();

		// Init sol
		auto sol = getBestHeuristic( G );
		MinVC.insert( sol.begin(), sol.end() );
		VC = MinVC;

		// Init weights
		nodeW.resize( G.vertices.size(), 0.0 );
		for ( auto &p : G.edges )
			edgeW[p.first] = 0.05;

		// Open trace ofstream
		ofs.open( outfile + ".trace" );
		TimePoint end = chrono::system_clock::now();
		SecondsDouble time_elapsed = end - start;
		ofs << time_elapsed.count() << ',' << MinVC.size() << '\n';

		ILS();

		// Close trace ofstream
		ofs.close();

		// Write best solution found
		ofstream solfs( outfile + ".sol" );
		solfs << MinVC.size() << '\n';
		res.insert( res.end(), MinVC.begin(), MinVC.end() );
		for ( size_t i = 0; i < res.size(); ++i ){
			if ( i != 0 )
				solfs << ',';
			solfs << ( res[i] + 1 );
		}

		return res;
	}
};

void localsearch1( Graph &G, int cutoff, int seed ){
	StochasticSolver solver( G, cutoff, seed );
	vector<Id> VC = solver.solve();

	G.check_coverage( VC );
}

// G					- The graph
// S					- Current solution
// tigthness			- Number of neighbors that are in S
// free					- Set of vertices with tightness 0, that are not in S
// newS/tightness/free	- Candidate solution and its data
// opt					- best solution
// gen					- random generator
// protection			- S cannot be replaced by worse solution for this many rounds
// cutoff				- cutoff time in seconds
// start				- starting time point of solver
// ofs					- ofstream to write to trace
class MISSolver{
private:
	Graph &G;
	unordered_set<Id> S, opt, newS, free, newfree;

	vector<int> tightness, newtightness;

	mt19937 gen;
	size_t protection;

	double cutoff;
	TimePoint start;
	ofstream ofs;
	string outfile;

public:
	MISSolver( Graph &G_, double cutoff_, int seed ) : G( G_ ), gen( seed ), cutoff( cutoff_ ){
		ostringstream oss;
		oss << "output/" << G.filename.substr( 0, G.filename.size() - 6u ) << "_LS2_" << cutoff << '_' << seed;
		outfile = oss.str();
	}

	void addToSol( Id i, unordered_set<Id> &s, unordered_set<Id> &f, vector<int> &t ){
		s.insert( i );
		f.erase( i );
		for ( Id j : G.vertices[i].neighs ){
			++t[j];
			if ( t[j] == 1 )
				f.erase( j );
		}
	}

	void remFromSol( Id i, unordered_set<Id> &s, unordered_set<Id> &f, vector<int> &t ){
		s.erase( i );
		f.insert( i );
		for ( Id j : G.vertices[i].neighs ){
			--t[j];
			if ( t[j] == 0 )
				f.insert( j );
		}
	}

	void perturb(){
		newS = S;

		// Determine k
		size_t k = 1;
		uniform_real_distribution<double> dis( 0.0, 1.0 );
		if ( protection == 0 && dis( gen ) <= 0.5/S.size() ){
			double roll = dis( gen );
			for ( double chance = 1.0; roll <= chance; chance *= 0.5 ) ++k;		
			k = min( k, S.size() );
		}

		newtightness = tightness;
		newfree = free;
		// S' <- S - { k random elements in S }
		vector<Id> temp( S.begin(), S.end() );
		shuffle( temp.begin(), temp.end(), gen );
		for ( size_t i = 0; i < k; ++i ){
			remFromSol( temp[i], newS, newfree, newtightness );
		}
		
		// Insert at most k new free vertices
		auto it = newfree.begin();
		uniform_int_distribution<size_t> idis( 0u, newfree.size()-1u );
		size_t adv = idis( gen );
		advance( it, adv );
		addToSol( *it, newS, newfree, newtightness );
		while ( --k && !newfree.empty() ){
			// Vertex in free has distance 2 from solution vertices if it has a neighbor with non-zero tightness
			for ( Id u : newfree ){
				bool good = false;
				for ( Id v : G.vertices[u].neighs ){
					good |= newtightness[v] != 0;
					if ( good ) break;
				}
				if ( good ){
					addToSol( u, newS, newfree, newtightness );
					break;
				}
			}
		}
	}

	void swapSols(){
		newS.swap( S );
		newfree.swap( free );
		newtightness.swap( tightness );
	}

	void two_improv(){
		// Initially all members of the solution are candidates
		unordered_set<Id> cand = newS;

		// While there are candidates
		while ( !cand.empty() ){
			Id x = *cand.begin();
			cand.erase( x );
			// Find replacements among x's neighbors that are 1-tight, not neighbors of each other
			for ( auto it1 = G.vertices[x].neighs.begin(), et = G.vertices[x].neighs.end(); it1 != et; ++it1 ){
				if ( newtightness[*it1] != 1 )
					continue;
				bool found = false;
				for ( auto it2 = it1; it2 != et; ++it2 ){
					if ( it1 == it2 || newtightness[*it2] != 1 
						 || G.vertices[*it1].neighs.find( *it2 ) != G.vertices[*it1].neighs.end() )
						continue;
					// Found
					remFromSol( x, newS, newfree, newtightness );
					addToSol( *it1, newS, newfree, newtightness );
					addToSol( *it2, newS, newfree, newtightness );

					// Add to candidates
					cand.insert( *it1 );
					cand.insert( *it2 );
					// If a neighbor of x became 1-tight due to x's removal
					for ( Id xn : G.vertices[x].neighs )
						if ( newtightness[xn] == 1 )
							// Then add the one neighbor that is in the solution to cand
							for ( Id xnn : G.vertices[xn].neighs )
								if ( xnn != *it1 && xnn != *it2 && newS.find( xnn ) != newS.end() )
									cand.insert( xnn );
					// Terminate loops
					found = true;
					break;
				}
				if ( found )
					break;
			}			
		}
	}

	void ILS(){
		// Try to improve initial solution
		swapSols();
		two_improv();
		swapSols();
		protection = S.size();

		TimePoint end = chrono::system_clock::now();
		for ( SecondsDouble elapsed_seconds = end - start; elapsed_seconds.count() <= cutoff; 
			end = chrono::system_clock::now(), elapsed_seconds = end - start ){			
			
			perturb();

			two_improv();

			// Decide to swap it
			if ( newS.size() > S.size() ){
				swapSols();
				protection = S.size();
				if ( S.size() > opt.size() ){
					opt = S;
					end = chrono::system_clock::now();
					elapsed_seconds = end - start;
					ofs << elapsed_seconds.count() << ',' << (G.vertices.size() - opt.size()) << '\n';
				}
			}
			else{
				if ( protection-- != 0 )
					continue;
				uniform_real_distribution<double> dis( 0.0, 1.0 );
				if ( dis( gen ) <= 1.0 / ( 1.0 + ( S.size() - newS.size() )*( opt.size() - newS.size() ) ) ){
					swapSols();
					protection = S.size();
				}
			}
		}		
	}

	vector<Id> solve(){
		vector<Id> res;

		// Solve
		start = chrono::system_clock::now();

		vector<Id> sol = getBestHeuristic( G );
		opt.reserve( G.vertices.size() );
		// Best solution
		for ( size_t i = 0; i < G.vertices.size(); ++i )
			opt.insert( i );
		for ( Id i : sol )
			opt.erase( i );
		// Current solution
		S = opt;

		// Initialize tightness
		tightness.resize( G.vertices.size(), 0 );
		for ( size_t i = 0; i < G.vertices.size(); ++i ){
			if ( S.find( i ) != S.end() )
				for ( Id j : G.vertices[i].neighs )
					++tightness[j];
		}
		// Initialize free
		for ( size_t i = 0; i < G.vertices.size(); ++i )
			if ( tightness[i] == 0 && S.find( i ) == S.end() )
				free.insert( i );

		// Open trace ofstream
		ofs.open( outfile + ".trace" );
		TimePoint end = chrono::system_clock::now();
		SecondsDouble time_elapsed = end - start;
		ofs << time_elapsed.count() << ',' << ( G.vertices.size() - opt.size() ) << '\n';
		ILS();

		// Close trace ofstream
		ofs.close();

		// Write best solution found
		ofstream solfs( outfile + ".sol" );
		solfs << ( G.vertices.size() - opt.size() ) << '\n';
		res.reserve( G.vertices.size() - opt.size() );
		for ( size_t i = 0; i < G.vertices.size(); ++i ){
			if ( opt.find( i ) != opt.end() )
				continue;
			if ( !res.empty() )
				solfs << ',';
			solfs << ( i + 1 );
			res.push_back( i );
		}

		return res;
	}
};

void localsearch2( Graph &G, int cutoff, int seed ){
	MISSolver solver( G, cutoff, seed );
	vector<Id> VC = solver.solve();

	G.check_coverage( VC );
}