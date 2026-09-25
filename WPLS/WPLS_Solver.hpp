#ifndef WPLS_SOLVER_HPP
#define WPLS_SOLVER_HPP

#include "GraphTypes.hpp"
#include <vector>
#include <queue>
#include <functional>
#include <limits>
#include <cmath>
#include <cstddef>

// Cost 構造体は GraphTypes.hpp で定義されているためここからは削除

struct Label {
    Cost cost;
    int pred_node;
    size_t pred_label_idx;
};

// 優先度付きキュー（Heap）で管理する候補
struct State {
    Cost cost;
    int node;
    int pred_node;
    size_t pred_label_idx;

    bool operator>(const State& other) const {
        return other.cost < cost;
    }
};

// solve() 実行中に観測されたキューの最大サイズ（診断用、定義は .cpp 側）
extern size_t max_queue_size;

// --- ソルバークラス ---
class WBDijkstra {
public:
    explicit WBDijkstra(int num_nodes);

    void add_edge(int u, int v, double c1, double c2);

    // L を解放して次の探索に備える
    void reset();

    // s から t への弱パレート解（確定ラベル集合）を求める
    std::vector<Label> solve(int s, int t);
    
    // target ノードの label_idx 番目のラベルから始点までの経路を復元する
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

    // 対象コストがノードのいずれかの確定ラベルに「強支配」されているか判定
    bool is_strongly_dominated(const std::vector<Label>& L_node, const Cost& target);

    void relaxation_process_with_pruning(int u);
    bool is_pruned(int node_idx, const Cost& cost);

    // source から（swap_costs / use_reverse_graph の指定に従い）辞書式最小コストを求める
    std::vector<Cost> run_lex_min_dijkstra(int source, bool swap_costs, bool use_reverse_graph);
};

#endif // WPLS_SOLVER_HPP