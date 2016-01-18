#include <string>
#include <iostream>
#include <random>
#include "graph.h"
#include "utils.h"
#include "bnb.h"
#include "heuristic.h"
#include "localsearch.h"

using namespace std;

class CmdlineParser{
private:
	int _argc;
	char** _argv;
public:
	CmdlineParser( int argc, char* argv[] ): _argc( argc ), _argv( argv ) {}
	
	int get_opt_intarg( string s, int default_val ){
		for ( int i = 0; i < _argc; ++i )
			if ( s == _argv[i] && i + 1 < _argc )
				return atoi( _argv[i + 1] );
		return default_val;
	}

	string get_opt_strarg( string s, string default_val ){
		for ( int i = 0; i < _argc; ++i )
			if ( s == _argv[i] && i + 1 < _argc )
				return string( _argv[i + 1] );
		return default_val;
	}
};


int main( int argc, char* argv[] ){
	CmdlineParser parser( argc, argv );

	Graph G = Graph( parser.get_opt_strarg( "-inst", "input.txt" ) );

    int cutoff = parser.get_opt_intarg( "-time", 600 );

    string method = parser.get_opt_strarg( "-alg", "BnB" );

    int seed = parser.get_opt_intarg( "-seed", 0 );

    if ( method == "BnB" )
		branch_and_bound( G, cutoff );
	else if ( method == "Approx" )
		heuristic( G );
	else if ( method == "LS1" )
		localsearch1( G, cutoff, seed );
	else if ( method == "LS2" )
		localsearch2( G, cutoff, seed );
    
    return 0;
}