#include "BDijkstra_Solver.hpp"
#include <iostream>
#include <chrono>
#include <algorithm>
#include <utility>

BDijkstra::BDijkstra(int num_nodes)
    : n(num_nodes),
    adj(num_nodes),
    rev_adj(num_nodes)
{}

void BDijkstra::add_edge(int u, int v, double c1, double c2) {
    adj[u].push_back({ v, {c1, c2} });
    rev_adj[v].push_back({ u, {c1, c2} });
}

std::vector<Label> BDijkstra::solve(int s, int t) {
    start_node = s;
    target_node = t;

    // 枝刈りのための事前計算
    dist_t = run_lex_min_dijkstra(t, false, true);
    wist_t = run_lex_min_dijkstra(t, true, true);
    auto s_to_t_lex1 = run_lex_min_dijkstra(s, false, false)[t];
    auto s_to_t_lex2 = run_lex_min_dijkstra(s, true, false)[t];
    nadir_s_t = { s_to_t_lex2.c1, s_to_t_lex1.c2 };

    // 初期化
    L.assign(n, std::vector<Label>());
    pq = std::priority_queue<State, std::vector<State>, std::greater<State>>();

    // 開始点をヒープに追加
    pq.push({ {0.0, 0.0}, s, -1, 0 });

    while (!pq.empty()) {
        State current = pq.top();
        pq.pop();
        int u = current.node;

        // ===== 修正点 1：取り出したラベルが、既に確定したラベルに支配されていないかチェック =====
        if (!L[u].empty() && L[u].back().cost.dominates(current.cost)) {
            continue;
        }
        // 同じコストのラベルがすでにある場合もスキップ
        if (!L[u].empty() && L[u].back().cost == current.cost) {
            continue;
        }

        // ラベルを永続化
        L[u].push_back({ current.cost, current.pred_node, current.pred_label_idx });

        // 新しい候補の探索と緩和処理
        new_candidate_label_with_pruning(u);
        relaxation_process_with_pruning(u);
    }

    return L[t];
}

std::vector<int> BDijkstra::reconstruct_path(int target, size_t label_idx) {
    std::vector<int> path;
    if (label_idx >= L[target].size()) return path;

    int current_node = target;
    size_t current_label_idx = label_idx;

    while (current_node != -1) {
        path.push_back(current_node);
        const Label& current_label = L[current_node][current_label_idx];
        current_node = current_label.pred_node;
        current_label_idx = current_label.pred_label_idx;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

void BDijkstra::new_candidate_label_with_pruning(int u) {
    Cost best_new_cost;
    std::pair<bool,State> new_candidate_state = { false, State() };
    const auto& last_perm_label_cost = L[u].back().cost;

    for (auto& rev_arc : rev_adj[u]) {
        int v = rev_arc.from;
        if (L[v].empty()) continue;

        for (size_t r = 0; r < L[v].size(); ++r) {
            const auto& pred_label = L[v][r];
            Cost new_cost = pred_label.cost + rev_arc.cost;

            if (last_perm_label_cost < new_cost && !last_perm_label_cost.dominates(new_cost)) {
                if (is_pruned(u, new_cost)) continue;
                if (new_cost < best_new_cost) {
                    best_new_cost = new_cost;
                    new_candidate_state = { true, {new_cost, u, v, r} };
                }
            }
        }
    }
    if (new_candidate_state.first) {
        pq.push(new_candidate_state.second);
    }
}

void BDijkstra::relaxation_process_with_pruning(int u) {
    const auto& last_perm_label = L[u].back();
    size_t pred_label_idx = L[u].size() - 1;

    for (const auto& arc : adj[u]) {
        int v = arc.to;
        Cost new_cost = last_perm_label.cost + arc.cost;

        // 辞書式順序の制約を撤廃。
        // 永続ラベルに支配されていなければ、候補としてヒープに追加する。
        if (!L[v].empty() && L[v].back().cost.dominates(new_cost)) {
            continue;
        }
        if (is_pruned(v, new_cost)) {
            continue;
        }
        // そのままヒープに追加
        pq.push({ new_cost, v, u, pred_label_idx });
    }
}

bool BDijkstra::is_pruned(int node_idx, const Cost& cost) {
    Cost estimated_total_cost = cost + Cost{ dist_t[node_idx].c1, wist_t[node_idx].c2 };
    if (estimated_total_cost.c1 > nadir_s_t.c1 || estimated_total_cost.c2 > nadir_s_t.c2) {
        return true;
    }
    if (!L[target_node].empty() && L[target_node].back().cost.dominates(estimated_total_cost)) {
        return true;
    }
    return false;
}

std::vector<Cost> BDijkstra::run_lex_min_dijkstra(int source, bool swap_costs, bool use_reverse_graph) {
    std::vector<Cost> dists(n, { INF, INF });
    dists[source] = { 0.0, 0.0 };

    using PQ_State = std::pair<Cost, int>;
    std::priority_queue<PQ_State, std::vector<PQ_State>, std::greater<PQ_State>> local_pq;
    local_pq.push({ {0.0, 0.0}, source });

    while (!local_pq.empty()) {
        Cost cost_from_pq = local_pq.top().first;
        int u = local_pq.top().second;
        local_pq.pop();

        Cost stored_cost_for_comp = swap_costs ? Cost{ dists[u].c2, dists[u].c1 } : dists[u];
        if (stored_cost_for_comp < cost_from_pq) {
            continue;
        }
        if (use_reverse_graph) {
            for (const auto& edge : rev_adj[u]) {
                int v = edge.from;
                Cost new_cost = dists[u] + edge.cost;
                Cost new_cost_for_comp = swap_costs ? Cost{ new_cost.c2, new_cost.c1 } : new_cost;
                Cost stored_v_cost_for_comp = swap_costs ? Cost{ dists[v].c2, dists[v].c1 } : dists[v];
                if (new_cost_for_comp < stored_v_cost_for_comp) {
                    dists[v] = new_cost;
                    local_pq.push({ new_cost_for_comp, v });
                }
            }
        }
        else {
            for (const auto& edge : adj[u]) {
                int v = edge.to;
                Cost new_cost = dists[u] + edge.cost;
                Cost new_cost_for_comp = swap_costs ? Cost{ new_cost.c2, new_cost.c1 } : new_cost;
                Cost stored_v_cost_for_comp = swap_costs ? Cost{ dists[v].c2, dists[v].c1 } : dists[v];
                if (new_cost_for_comp < stored_v_cost_for_comp) {
                    dists[v] = new_cost;
                    local_pq.push({ new_cost_for_comp, v });
                }
            }
        }
    }
    return dists;
}