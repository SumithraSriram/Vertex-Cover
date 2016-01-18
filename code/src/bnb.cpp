#include "bnb.h"
#include "simplex.h"
#include "heuristic.h"
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <chrono>

using namespace std;

using TimePoint = chrono::time_point<std::chrono::system_clock>;
using SecondsDouble = chrono::duration<double>;

// G		- The graph
// S		- Set of vertices that are to be explored
// uncov	- number of edges still uncovered
// soln		- holds the local solution
// opt		- holds the optimal solution
// cutoff	- cutoff time in seconds
// start	- starting time point of solver
// ofs		- ofstream to write to trace
class BnBSolver{
private:
	Graph &G;
	vector<Id> soln, opt;
	unordered_set<Id> S;
	Id uncov;
	double cutoff;
	TimePoint start;
	ofstream ofs;

public:
	BnBSolver( Graph &G_, double cutoff_ ): G(G_), uncov( G.edges.size() ), cutoff( cutoff_ ){
	}

	void branch(){
		TimePoint end = chrono::system_clock::now();
		SecondsDouble elapsed_seconds = end - start;
		if ( elapsed_seconds.count() > cutoff )
			return;
		// Recursion exit condition. All covered
		if ( uncov == 0 ){
			if ( soln.size() < opt.size() ){
				opt = soln;
				ofs << elapsed_seconds.count() << ',' << opt.size() << '\n';
			}
			return;
		}

		// Is it worth going further? Check lower bound
		// Use worse algorithm, because it gives better lower bound. (We want higher number)
		// THis is why ---> H/2 <= OPT <= H <= 2OPT
		vector<Id> x = greedybad( G );
		size_t low = x.size() / 2;
		if ( soln.size() + low >= opt.size() || low > S.size() )
			return;

		// Get next considered vertex
		Id u = *max_element( S.begin(), S.end(),
							 [&]( Id a, Id b ){
			return G.vertices[a].neighs.size() < G.vertices[b].neighs.size();
		} );
		S.erase( u );

		////////////////////////////////////////////////////////////////////////////////////
		// Case 1: Add it to the solution
		// Skip this case if u has degree 0 or 1, but its neighbor has more than 1 degree
		if ( !G.vertices[u].neighs.empty() && !( G.vertices[u].neighs.size() == 1 && G.vertices[*G.vertices[u].neighs.begin()].neighs.size() > 1 ) ){
			soln.push_back( u );
			// Erase edges from graph, no need to delete the vertex itself though
			for ( Id v : G.vertices[u].neighs ){
				G.edges.erase( edge( u, v ) );
				G.vertices[v].neighs.erase( u );
				--uncov;
			}
			auto backup = move( G.vertices[u].neighs );
			G.vertices[u].neighs.clear();
			// Branch
			branch();
			// Restore graph and solution
			G.vertices[u].neighs = move( backup );
			for ( Id v : G.vertices[u].neighs ){
				G.edges[edge( u, v )] = false;
				G.vertices[v].neighs.insert( u );
				++uncov;
			}
			soln.pop_back();
		}
		////////////////////////////////////////////////////////////////////////////////////
		// Case 2: Don't add it to the solution
		branch();

		// Insert it back to considered vertices
		S.insert( u );
	}

	vector<Id> solve(){
		ostringstream oss;
		oss << "output/" << G.filename.substr( 0, G.filename.size() - 6u ) << "_BnB_" << cutoff;
		
		// Solve
		start = chrono::system_clock::now();

		// Initial solution
		opt.reserve( G.vertices.size() );
		for ( Id i = 0; i < G.vertices.size(); ++i )
			opt.push_back( i );

		// Open trace ofstream
		ofs.open( oss.str() + ".trace" );
		TimePoint end = chrono::system_clock::now();
		SecondsDouble time_elapsed = end - start;
		ofs << time_elapsed.count() << ',' << opt.size() << '\n';

		// Initialize vertices to be explored
		S.reserve( G.vertices.size() );
		for ( Id i = 0; i < G.vertices.size(); ++i )
			if ( G.vertices[i].neighs.size() > 1 )
				S.insert( i );
		branch();

		// Close trace ofstream
		ofs.close();

		// Write best solution found
		ofstream solfs( oss.str() + ".sol" );
		solfs << opt.size() << '\n';
		for ( size_t i = 0; i < opt.size(); ++i ){
			if ( i != 0 )
				solfs << ',';
			solfs << ( opt[i] + 1 );
		}

		return opt;
	}
};

void branch_and_bound( Graph &G, int cutoff ){
	BnBSolver solver( G, cutoff );
	vector<Id> VC = solver.solve();

	G.check_coverage( VC );
}