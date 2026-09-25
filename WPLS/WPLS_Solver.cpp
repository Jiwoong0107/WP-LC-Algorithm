#include "WPLS_Solver.hpp"
#include <iostream>
#include <chrono>
#include <algorithm>
#include <utility>

size_t max_queue_size = 0;

WBDijkstra::WBDijkstra(int num_nodes)
    : n(num_nodes),
    adj(num_nodes),
    rev_adj(num_nodes)
{}

void WBDijkstra::add_edge(int u, int v, double c1, double c2) {
    adj[u].push_back({ v, {c1, c2} });
    rev_adj[v].push_back({ u, {c1, c2} });
}

void WBDijkstra::reset() {
    // swap を使って、前回の探索で肥大化したメモリを OS に完全返却する
    std::vector<std::vector<Label>>().swap(L);
}

std::vector<Label> WBDijkstra::solve(int s, int t) {
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

    // 始点をヒープに追加
    pq.push({ {0.0, 0.0}, s, -1, 0 });

    long long iter = 0;
    // タイムアウト判定用の開始時刻を記録
    auto time_limit_start = std::chrono::steady_clock::now();

    while (!pq.empty()) {
        State current = pq.top();
        pq.pop();
        int u = current.node;

        if (iter % 1000 == 0) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed_minutes = std::chrono::duration_cast<std::chrono::minutes>(now - time_limit_start).count();
            if (elapsed_minutes >= 60) {
                std::cout << "\n[Timeout] Search exceeded 60 minutes. Terminating...\n";
                return L[t]; // タイムアウトによる強制終了
            }
        }

        if (max_queue_size < pq.size()) {
            max_queue_size = pq.size();
        }

        // 弱パレート判定：L[u] の確定ラベルすべてと比較し、強支配されていないか確認
        if (is_strongly_dominated(L[u], current.cost)) {
            continue;
        }

        // 完全に同じコスト・同じ経路由来の重複をスキップ
        // （別ルートで同じコストの場合は別の弱パレート解として残す）
        bool duplicate = false;
        for (const auto& l : L[u]) {
            if (l.cost == current.cost && l.pred_node == current.pred_node && l.pred_label_idx == current.pred_label_idx) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) continue;

        // ラベルを確定
        L[u].push_back({ current.cost, current.pred_node, current.pred_label_idx });

        // 緩和処理
        relaxation_process_with_pruning(u);
        iter++;
    }

    return L[t];
}

std::vector<int> WBDijkstra::reconstruct_path(int target, size_t label_idx) {
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

bool WBDijkstra::is_strongly_dominated(const std::vector<Label>& L_node, const Cost& target) {
    for (const auto& label : L_node) {
        if (label.cost.strong_dominates(target)) {
            return true;
        }
    }
    return false;
}

void WBDijkstra::relaxation_process_with_pruning(int u) {
    const auto& last_perm_label = L[u].back();
    size_t pred_label_idx = L[u].size() - 1;

    for (const auto& arc : adj[u]) {
        int v = arc.to;
        Cost new_cost = last_perm_label.cost + arc.cost;

        // 遷移先ノード v の確定ラベルのいずれかに「強支配」されていればスキップ
        if (is_strongly_dominated(L[v], new_cost)) {
            continue;
        }
        if (is_pruned(v, new_cost)) {
            continue;
        }

        // 候補としてヒープに追加
        pq.push({ new_cost, v, u, pred_label_idx });
    }
}

bool WBDijkstra::is_pruned(int node_idx, const Cost& cost) {
    Cost estimated_total_cost = cost + Cost{ dist_t[node_idx].c1, wist_t[node_idx].c2 };
    if (estimated_total_cost.c1 >= nadir_s_t.c1 && estimated_total_cost.c2 >= nadir_s_t.c2) {
        return true;
    }

    // ターゲットノードの確定ラベルに強支配されていれば枝刈り
    if (is_strongly_dominated(L[target_node], estimated_total_cost)) {
        return true;
    }

    return false;
}

std::vector<Cost> WBDijkstra::run_lex_min_dijkstra(int source, bool swap_costs, bool use_reverse_graph) {
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