#ifndef GRAPH_GENERATORS_HPP
#define GRAPH_GENERATORS_HPP

#include <string>
#include "WPLC_Solver.hpp"

void setup_dimacs_graph(DiagnosticBiObjectiveBellmanFord& bd, const std::string& filenameD, const std::string& filenameT);
void setup_test_graph_with_2_cycles(DiagnosticBiObjectiveBellmanFord& bd);
void setup_research_based_dense_graph(DiagnosticBiObjectiveBellmanFord& bd, int n, int density_m_per_n, int max_chord_length, double correlation, bool allow_negative);
void setup_dimacs_graph_with_elevation(DiagnosticBiObjectiveBellmanFord& bd, const std::string& fD, const std::string& fT, int num_nodes);
void setup_missing_solution_graph(DiagnosticBiObjectiveBellmanFord& bd);
void setup_random_graph_negative_no_cycle(DiagnosticBiObjectiveBellmanFord& bd, int num_nodes, int num_edges, int base_min, int base_max, int potential_scale);
void setup_hidden_cycles_graph(DiagnosticBiObjectiveBellmanFord& bd);
void setup_grid_graph1(DiagnosticBiObjectiveBellmanFord& bd, int height, int width, int lower_limit, int upper_limit);
void setup_grid_graph(DiagnosticBiObjectiveBellmanFord& bd, int width, int height);
void setup_grid_graph_paper(DiagnosticBiObjectiveBellmanFord& bd, int height, int width, int s, int t, int i);
void setup_nwmoa_killer_graph(DiagnosticBiObjectiveBellmanFord& bd);

#endif // GRAPH_GENERATORS_HPP