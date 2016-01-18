#ifndef HEURISTIC_H
#define HEURISTIC_H

#include "graph.h"

vector<Id> greedybad( Graph &G );

vector<Id> heuristic1( Graph &G );

vector<Id> heuristic2( Graph &G );

void heuristic( Graph &G );

vector<Id> getBestHeuristic( Graph &G );

#endif
