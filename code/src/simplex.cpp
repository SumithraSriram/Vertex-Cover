/*
Simplex in C++.
A modification of moshahmed/at/gmail's code that allows minimization.
http://stackoverflow.com/a/15080146
Turning maximization simplex into minimization:
http://college.cengage.com/mathematics/larson/elementary_linear/4e/shared/downloads/c09s4.pdf
*/

#include <cassert>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include "simplex.h"

static const double epsilon = 1.0e-8;
inline bool equal( double a, double b ) { return fabs( a - b ) < epsilon; }

using namespace std;

Tableau::Tableau() {}

Tableau::Tableau( Mode m, vector<double> vars ) : obj( m ){
    varnum = vars.size();
    if ( m == Mode::MAX ){
        mat.emplace_back( vars.size() + 1 );
        auto it = mat[0].begin();
        *it = 0;
        copy( vars.begin(), vars.end(), ++it );
    }
    else {
        for ( size_t i = 0; i <= vars.size(); ++i )
            mat.emplace_back( 1, 0.0 );
        for ( size_t i = 0; i < vars.size(); ++i )
            mat[i + 1][0] = vars[i];
    }
}

void Tableau::add_constraint( vector<double> vars ){
    if ( obj == Mode::MAX ){
        if ( vars.size() != col_size() ) return;
        mat.push_back( vars );
    }
    else{
        if ( vars.size() != row_size() ) return;
        for ( size_t i = 0; i < vars.size(); ++i )
            mat[i].push_back( vars[i] );
    }
}

size_t Tableau::row_size(){ return mat.size(); }

size_t Tableau::col_size(){ return mat[0].size(); }

void nl( int k ){
    cout << string( k, '-' ) << '\n';
}

void print_tableau( Tableau &tab, const string &mes ) {
    static int counter = 0;
    cout << '\n' << ++counter << ". Tableau " << mes << ":\n";
    nl( 70 );

    cout << "col:\tb[i]\t";
    for ( size_t j = 1; j < tab.col_size(); j++ )
        cout << 'x' << j << '\t';
    cout << '\n';

    for ( size_t i = 0; i < tab.row_size(); i++ ) {
        if ( i == 0 )
            cout << "max:";
        else
            cout << 'b' << i << ": ";
        cout << '\t';
        for ( size_t j = 0; j < tab.col_size(); j++ )
            cout << setprecision( 2 ) << tab.mat[i][j] << '\t';
        cout << '\n';
    }
    nl( 70 );
}

Tableau read_tableau( const string &filename ) {
    Tableau tab;
    ifstream ifs( filename );

    size_t m, n;
    ifs >> m >> n;


    for ( size_t i = 0; i < m; i++ ) {
        tab.mat.emplace_back( n );
        for ( size_t j = 0; j < n; j++ ) {
            ifs >> tab.mat[i][j];
        }
    }
    //cout << "Read tableau [" << m << " rows x " << n << "columns] from file '" << filename << "'.\n";
    return tab;
}

void pivot_on( Tableau &tab, const size_t row, const size_t col ) {
    double pivot = tab.mat[row][col];
    assert( pivot > 0 );
    for ( size_t j = 0; j < tab.col_size(); j++ )
        tab.mat[row][j] /= pivot;
    assert( equal( tab.mat[row][col], 1. ) );

    for ( size_t i = 0; i < tab.row_size(); i++ ) { // foreach remaining row i do
        double multiplier = tab.mat[i][col];
        if ( i == row ) continue;
        for ( size_t j = 0; j < tab.col_size(); j++ ) { // r[i] = r[i] - z * r[row];
            tab.mat[i][j] -= multiplier * tab.mat[row][j];
        }
    }
}

// Find pivot_col = most negative column in mat[0][1..n]
int find_pivot_column( Tableau &tab ) {
    int pivot_col = 1;
    double lowest = tab.mat[0][pivot_col];
    for ( size_t j = 1; j < tab.col_size(); j++ ) {
        if ( tab.mat[0][j] < lowest ) {
            lowest = tab.mat[0][j];
            pivot_col = j;
        }
    }
    //cout << "Most negative column in row[0] is col " << pivot_col << " = " << setprecision( 2 ) << lowest << ".\n";
    if ( lowest >= 0 ) {
        return -1; // All positive columns in row[0], this is optimal.
    }
    return pivot_col;
}

// Find the pivot_row, with smallest positive ratio = col[0] / col[pivot]
int find_pivot_row( Tableau &tab, const int pivot_col ) {
    int pivot_row = 0;
    bool first = true;
    double min_ratio = 0.0;
    //cout << "Ratios A[row_i,0]/A[row_i," << pivot_col << "] = [";
    for ( size_t i = 1; i < tab.row_size(); i++ ){
        if ( tab.mat[i][pivot_col] <= epsilon )
            continue;
        double ratio = tab.mat[i][0] / tab.mat[i][pivot_col];
        //cout << ratio << ", ";
        if ( ( ratio > 0 && ratio < min_ratio ) || first ) {
            min_ratio = ratio;
            pivot_row = i;
            first = false;
        }
    }
    //cout << "].\n";
    if ( first )
        return -1; // Unbounded.
    //cout << "Found pivot A[" << pivot_row << ", " << pivot_col << "], min positive ratio="
        //<< setprecision( 2 ) << min_ratio << " in row=" << pivot_row << ".\n";
    return pivot_row;
}

void add_slack_variables( Tableau &tab ) {
    size_t targetsize = tab.col_size() + tab.row_size() - 1;
    for ( size_t i = 0; i < tab.row_size(); i++ ){
        tab.mat[i].resize( targetsize, 0.0 );
        if ( i != 0 )
            tab.mat[i][i + targetsize - tab.row_size()] = 1.0;
    }
}

void check_b_positive( Tableau &tab ) {
    for ( size_t i = 1; i < tab.row_size(); i++ )
        assert( tab.mat[i][0] >= 0 );
}

// Given a column of identity matrix, find the row containing 1.
// return -1, if the column as not from an identity matrix.
int find_basis_variable( Tableau &tab, const int col ) {
    int xi = -1;
    for ( size_t i = 1; i < tab.row_size(); i++ ) {
        if ( equal( tab.mat[i][col], 1 ) ) {
            if ( xi == -1 )
                xi = i;   // found first '1', save this row number.
            else
                return -1; // found second '1', not an identity matrix.

        }
        else if ( !equal( tab.mat[i][col], 0 ) ) {
            return -1; // not an identity matrix column.
        }
    }
    return xi;
}

vector<double> get_solution( Tableau &tab ){
    vector<double> sol( tab.varnum, 0.0 );
    if ( tab.obj == Mode::MAX ){
        for ( size_t j = 1; j < tab.varnum; j++ ) { // for each column.
            int xi = find_basis_variable( tab, j );
            if ( xi != -1 )
                sol[j - 1] = tab.mat[xi][0];
        }
    }
    else{
        for ( size_t i = tab.col_size() - tab.varnum; i < tab.col_size(); ++i )
            sol[i + tab.varnum - tab.col_size()] = tab.mat[0][i];
    }
    return sol;
}

void print_optimal_vector( Tableau &tab, const string &message ) {
    if ( tab.obj == Mode::MAX ){
        cout << message << " at ";
        for ( size_t j = 1; j <= tab.varnum; j++ ) { // for each column.
            int xi = find_basis_variable( tab, j );
            if ( xi != -1 )
                cout << 'x' << j << '=' << setprecision( 2 ) << tab.mat[xi][0] << ", ";
            else
                cout << 'x' << j << "=0, ";
        }
        cout << '\n';
    }
    else{
        cout << message << " at ";
        for ( size_t i = tab.col_size() - tab.varnum; i < tab.col_size(); ++i )
            cout << 'x' << i + tab.varnum - tab.col_size() + 1 << '=' << setprecision( 2 ) << tab.mat[0][i] << ", ";
        cout << '\n';
    }
}

void transpose_tableau( Tableau &tab ){
    vector<vector<double>> newmat( tab.col_size(), vector<double>( tab.row_size() ) );
    for ( size_t r = 0; r < tab.row_size(); ++r )
        for ( size_t c = 0; c < tab.col_size(); ++c )
            newmat[c][r] = tab.mat[r][c];
    tab.mat = move( newmat );
}

bool simplex( Tableau &tab ) {
    int maxiter = max( tab.col_size(), tab.row_size() );
    for ( size_t i = 0; i < tab.col_size(); ++i )
        tab.mat[0][i] *= -1.0;
    add_slack_variables( tab );
    check_b_positive( tab );
    //print_tableau( tab, "Padded with slack variables" );
    while ( maxiter-- ) {
        int pivot_col, pivot_row;

        pivot_col = find_pivot_column( tab );
        if ( pivot_col < 0 ) {
            //cout << "Found optimal value=A[0,0]=" << setprecision( 2 ) << tab.mat[0][0] << " (no negatives in row 0).\n";
            //print_optimal_vector( tab, "Optimal vector" );
            return true;
        }
        //cout << "Entering variable x" << pivot_col << " to be made basic, so pivot_col=" << pivot_col << ".\n";

        pivot_row = find_pivot_row( tab, pivot_col );
        if ( pivot_row < 0 ) {
            //cout << "unbounded (no pivot_row).\n";
            return false;
        }
        //cout << "Leaving variable x" << pivot_row << ", so pivot_row=" << pivot_row << "\n";

        pivot_on( tab, pivot_row, pivot_col );
        //print_tableau( tab, "After pivoting" );
        //print_optimal_vector( tab, "Basic feasible solution" );
    }
    //cout << "Too many iterations > " << max( tab.col_size(), tab.row_size() ) << ".\n";
    return false;
}

//
//
//int main( int argc, char *argv[] ){
//    Tableau tab;
//    if ( argc > 8 )  // usage: cmd datafile
//        tab = read_tableau( argv[1] );
//    else{
//        tab.obj = Mode::MIN;
//        tab.varnum = 2;
//        tab.mat = {                    
//            { 0.0 , 3.0 , 2.0, },    // Min: 3x + 2y,
//            { 6.0 , 2.0 , 1.0, },    //      2x + y >= 6 .. b1
//            { 4.0 , 1.0 , 1.0, } };  //      x + y >= 4 .. b2
//                                     // Sol: x1 = 2, x2 = 2
//		transpose_tableau( tab );
//
//        tab.obj = Mode::MAX;
//        tab.varnum = 2;
//        tab.mat = {                     
//            { 0.0 , 7.0 , 10.0, },    // Max: 7x + 10y,
//            { 9.0 , 1.0 , 9.0, },     //      x + 9y <= 9 .. b1
//            { 12.0 , 1.0 , 2.0, },    //      x + 2y <= 10 .. b2
//            { 20.0 , 1.0 , 4.0, } };  //      x + 4y <= 12 .. b3
//                                      // Sol: x1 = 9, x2 = 0
//
//        tab.obj = Mode::MIN;
//        tab.varnum = 3;
//        tab.mat = {                     
//        { 0.0 , 1.0 , 1.0, 1.0 },  // Min: x + y + z,
//        { 1.0 , 0.0 , 1.0, 1.0 },  //      x + y >= 1 .. b1
//        { 1.0 , 1.0 , 0.0, 1.0 },  //      x + z >= 1 .. b2
//        { 1.0 , 1.0 , 1.0, 0.0 }   //      y + z >= 1 .. b3
//        };
//		transpose_tableau( tab );
//    }
//    print_tableau( tab, "Initial" );
//    simplex( tab );
//    return 0;
//}

