#ifndef BDIJKSTRA_HPP
#define BDIJKSTRA_HPP

#include "GraphTypes.hpp" // Cost, Arc, ReverseArc 等を利用
#include <vector>
#include <queue>
#include <functional>
#include <limits>
#include <cmath>
#include <cstddef>

// 永続ラベル（Permanent Label）を表す構造体
struct Label {
    Cost cost;
    int pred_node;
    size_t pred_label_idx;
};

// 優先度付きキュー（Heap）で管理する状態
struct State {
    Cost cost;
    int node;
    int pred_node;
    size_t pred_label_idx;

    bool operator>(const State& other) const {
        return other.cost < cost;
    }
};

class BDijkstra {
public:
    BDijkstra(int num_nodes);

    void add_edge(int u, int v, double c1, double c2);

    std::vector<Label> solve(int s, int t);

    std::vector<int> reconstruct_path(int target, size_t label_idx);

private:
    int n;
    int start_node, target_node;
    std::vector <std::vector<Arc>> adj;
    std::vector <std::vector<ReverseArc>> rev_adj;
    std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
    std::vector<std::vector<Label>> L;

    // 枝刈り用データ
    std::vector<Cost> dist_t;
    std::vector<Cost> wist_t;
    Cost nadir_s_t;

    void new_candidate_label_with_pruning(int u);
    void relaxation_process_with_pruning(int u);
    bool is_pruned(int node_idx, const Cost& cost);
    std::vector<Cost> run_lex_min_dijkstra(int source, bool swap_costs, bool use_reverse_graph);
};

#endif // BDIJKSTRA_HPP