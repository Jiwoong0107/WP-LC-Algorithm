#ifndef BDIJKSTRA_GRAPH_GENERATORS_HPP
#define BDIJKSTRA_GRAPH_GENERATORS_HPP

#include <string>
#include "BDijkstra_Solver.hpp"

// DIMACS 形式の距離ファイル・時間ファイルから BDijkstra 用のグラフを構築する
void setup_dimacs_graph(BDijkstra& bd, const std::string& filenameD, const std::string& filenameT);

// 論文仕様のグリッドグラフを構築する（BDijkstra は Dijkstra ベースのため、
// 辺重みは 1〜10 の正の範囲に限定している点に注意）
void setup_grid_graph_paper(BDijkstra& bd, int height, int width, int s, int t, int i);

#endif // BDIJKSTRA_GRAPH_GENERATORS_HPP#pragma once
