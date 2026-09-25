#ifndef GRAPH_GENERATORS_HPP
#define GRAPH_GENERATORS_HPP

#include <string>
#include "WPLC_Solver.hpp"

void setup_dimacs_graph(DiagnosticBiObjectiveBellmanFord& bd, const std::string& filenameD, const std::string& filenameT);
void setup_dimacs_graph_with_elevation(DiagnosticBiObjectiveBellmanFord& bd, const std::string& fD, const std::string& fT, int num_nodes);
void setup_grid_graph_paper(DiagnosticBiObjectiveBellmanFord& bd, int height, int width, int s, int t, int i);

#endif // GRAPH_GENERATORS_HPP