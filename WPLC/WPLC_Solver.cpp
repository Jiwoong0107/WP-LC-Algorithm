#include "WPLC_Solver.hpp"
#include <iostream>
#include <algorithm>
#include <chrono>
#include <deque>

size_t max_queue_size = 0;

DiagnosticBiObjectiveBellmanFord::DiagnosticBiObjectiveBellmanFord(int num_nodes) : n(num_nodes), adj(num_nodes), rev_adj(num_nodes) {
    temp_visited.assign(num_nodes, false);
}

void DiagnosticBiObjectiveBellmanFord::add_edge(int u, int v, long long c1, long long c2) {
    if (u < 0 || u >= n || v < 0 || v >= n) return;
    adj[u].push_back({ v, {c1, c2} });
    rev_adj[v].push_back({ u, {c1, c2} });
}

bool DiagnosticBiObjectiveBellmanFord::solve(int s, int t) {
    label_pool.clear();
    free_list.clear();
    fast_node_id.clear();
    fast_parent_idx.clear();

    L_indices.assign(n, {});
    U_indices.assign(n, {});

    problematic_cycles.clear();
    target_node_id = t;
    pruned_count = 0;

    std::cout << "Calculating Heuristics & Danger Zones..." << std::endl;
    calculate_heuristics_and_danger_zones(s, t);

    if (heuristics_bwd[s].c1 >= INF_COST || heuristics_bwd[s].c2 >= INF_COST) {
        std::cout << "Target is unreachable from start.\n";
        return false;
    }

    std::deque<int> Q;
    std::vector<bool> in_queue(n, false);

    int start_idx = create_label(Cost{ 0, 0 }, s, -1);
    L_indices[s].push_back(start_idx);
    U_indices[s].push_back(start_idx);
    label_pool[start_idx].ref_count++;

    Q.push_back(s);
    in_queue[s] = true;

    long long iter = 0;
    std::cout << "Starting SPFA search (Hybrid: GC + Danger-Zone Cycle Check)..." << std::endl;

    auto time_limit_start = std::chrono::steady_clock::now();

    while (!Q.empty()) {
        if (iter % 1000 == 0) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed_minutes = std::chrono::duration_cast<std::chrono::minutes>(now - time_limit_start).count();
            if (elapsed_minutes >= 60) {
                std::cout << "\n[Timeout] Search exceeded 60 minutes. Terminating...\n";
                return false;
            }
        }

        if (max_queue_size < Q.size()) {
            max_queue_size = Q.size();
        }
        int u = Q.front();
        Q.pop_front();
        in_queue[u] = false;

        std::vector<int> labels_to_propagate = std::move(U_indices[u]);
        U_indices[u].clear();

        for (const auto& edge : adj[u]) {
            int v = edge.to;
            bool changed = false;

            bool is_non_positive_edge = (edge.cost.c1 <= 0 || edge.cost.c2 <= 0);
            long long min_f1_v = INF_COST;
            long long min_f2_v = INF_COST;

            for (int l_idx : labels_to_propagate) {
                if (!label_pool[l_idx].is_valid) continue;

                Cost next_c = label_pool[l_idx].cost + edge.cost;

                if (is_pruned(v, next_c)) {
                    pruned_count++;
                    continue;
                }

                // 辺コストが正の場合は閉路検知をコメントアウト
                //if (CheckAndRecordCycle(l_idx, v, edge)) continue;
                

                int new_idx = UpdateWeakParetoSet(v, next_c, l_idx);
                if (new_idx != -1) {
                    U_indices[v].push_back(new_idx);
                    label_pool[new_idx].ref_count++;
                    changed = true;

                    long long current_f1_v = (next_c.c1 + heuristics_bwd[v].c1);
                    long long current_f2_v = (next_c.c2 + heuristics_bwd[v].c2);
                    if (current_f1_v < min_f1_v) {
                        min_f1_v = current_f1_v;
                    }
                    if (current_f2_v < min_f2_v) {
                        min_f2_v = current_f2_v;
                    }
                }
            }

            if (changed && !in_queue[v]) {
                in_queue[v] = true;
                bool push_front = false;
                if (!Q.empty()) {
                    int f_node = Q.front();
                    long long f1_front = INF_COST;
                    long long f2_front = INF_COST;

                    for (int idx : L_indices[f_node]) {
                        if (label_pool[idx].is_valid) {
                            f1_front = (label_pool[idx].cost.c1 + heuristics_bwd[f_node].c1);
                            f2_front = (label_pool[idx].cost.c2 + heuristics_bwd[f_node].c2);
                            break;
                        }
                    }

                    if (f1_front == INF_COST || f2_front == INF_COST) {
                        f1_front = heuristics_bwd[f_node].c1;
                        f2_front = heuristics_bwd[f_node].c2;
                    }

                    if (min_f1_v < f1_front && min_f2_v < f2_front) {
                        push_front = true;
                    }
                }

                if (push_front) {
                    Q.push_front(v);
                }
                else {
                    Q.push_back(v);
                }
            }
        }

        for (int l_idx : labels_to_propagate) {
            release_ref(l_idx);
        }

        iter++;
    }
    std::cout << "\nSearch finished in " << iter << " queue pops.\n";
    std::cout << "Pruned operations during search: " << pruned_count << "\n";
    return true;
}

void DiagnosticBiObjectiveBellmanFord::print_diagnosis(int t, bool cyc_result, bool sol_result) {
    std::cout << "\n=== [Diagnosis] Problematic Cycles Found ===\n";
    if (problematic_cycles.empty()) {
        std::cout << "No problematic cycles detected.\n";
    }
    else {
        int cycle_count = 1;
        if (cyc_result) {
            for (const auto& cyc : problematic_cycles) {
                std::cout << "Cycle " << cycle_count++ << " [" << cyc.type << "] Delta("
                    << cyc.total_cost.c1 << ", " << cyc.total_cost.c2 << ") | Path: ";
                for (size_t i = 0; i < cyc.path.size(); ++i) {
                    std::cout << cyc.path[i] << (i + 1 == cyc.path.size() ? "" : "->");
                }
                std::cout << "\n";
            }
        }
    }

    std::cout << "\n=== [Result] Weakly Pareto Solutions for Target " << t << " ===\n";
    if (L_indices[t].empty()) {
        std::cout << "No path found.\n";
        return;
    }

    int count = 0;
    for (int idx : L_indices[t]) {
        if (!label_pool[idx].is_valid) continue;
        count++;

        std::vector<int> path;
        int current_idx = idx;
        while (current_idx != -1) {
            path.push_back(fast_node_id[current_idx]);
            current_idx = fast_parent_idx[current_idx];
        }
        std::reverse(path.begin(), path.end());

        if (sol_result) {
            std::cout << "Solution " << count << " | Cost(" << label_pool[idx].cost.c1 << ", " << label_pool[idx].cost.c2 << ") | Route: ";
            for (size_t i = 0; i < path.size(); ++i) {
                std::cout << path[i] << (i + 1 == path.size() ? "" : "->");
            }
            std::cout << "\n";
        }
    }

    std::cout << "\n( weak pareto solutions , cycles ) = (" << count << " , " << problematic_cycles.size() << " )\n";
}

int DiagnosticBiObjectiveBellmanFord::get_weak_pareto_count(int t) const {
    if (L_indices.empty() || L_indices[t].empty()) return 0;
    int count = 0;
    for (int idx : L_indices[t]) {
        if (label_pool[idx].is_valid) count++;
    }
    return count;
}

void DiagnosticBiObjectiveBellmanFord::print_graph() {
    std::cout << "\n--- Graph Edges ---\n";
    for (int u = 0; u < n; ++u) {
        for (const auto& e : adj[u]) {
            std::cout << u << "->" << e.to << " : (" << e.cost.c1 << ", " << e.cost.c2 << ")\n";
        }
    }
    std::cout << "-------------------\n";
}

int DiagnosticBiObjectiveBellmanFord::create_label(const Cost& c, int id, int parent_idx) {
    int new_idx;
    if (!free_list.empty()) {
        new_idx = free_list.back();
        free_list.pop_back();
        label_pool[new_idx] = Label(c, id, parent_idx);
    }
    else {
        new_idx = (int)label_pool.size();
        label_pool.emplace_back(c, id, parent_idx);
    }

    if (new_idx >= (int)fast_node_id.size()) {
        fast_node_id.resize(std::max(10000, (int)fast_node_id.size() * 2));
        fast_parent_idx.resize(std::max(10000, (int)fast_parent_idx.size() * 2));
    }
    fast_node_id[new_idx] = id;
    fast_parent_idx[new_idx] = parent_idx;

    if (parent_idx != -1) {
        label_pool[parent_idx].ref_count++;
    }
    return new_idx;
}

void DiagnosticBiObjectiveBellmanFord::release_ref(int idx) {
    int cur = idx;
    while (cur != -1) {
        label_pool[cur].ref_count--;
        if (label_pool[cur].ref_count == 0 && !label_pool[cur].is_valid) {
            int p_idx = label_pool[cur].parent_index;
            free_list.push_back(cur);
            cur = p_idx;
        }
        else break;
    }
}

bool DiagnosticBiObjectiveBellmanFord::can_exist_dominate_new(int exist_idx, int p_idx, int v) {
    if (!is_danger[v] || v == target_node_id) return true;

    int cur_new = p_idx;
    while (cur_new != -1) {
        int node = fast_node_id[cur_new];
        if (is_danger[node]) temp_visited[node] = true;
        cur_new = fast_parent_idx[cur_new];
    }
    if (is_danger[v]) temp_visited[v] = true;

    bool is_subset = true;
    int cur_exist = exist_idx;
    while (cur_exist != -1) {
        int node = fast_node_id[cur_exist];
        if (is_danger[node] && !temp_visited[node]) {
            is_subset = false;
            break;
        }
        cur_exist = fast_parent_idx[cur_exist];
    }

    cur_new = p_idx;
    while (cur_new != -1) {
        int node = fast_node_id[cur_new];
        if (is_danger[node]) temp_visited[node] = false;
        cur_new = fast_parent_idx[cur_new];
    }
    if (is_danger[v]) temp_visited[v] = false;

    return is_subset;
}

bool DiagnosticBiObjectiveBellmanFord::can_new_dominate_exist(int p_idx, int v, int exist_idx) {
    if (!is_danger[v] || v == target_node_id) return true;

    int cur_exist = exist_idx;
    while (cur_exist != -1) {
        int node = fast_node_id[cur_exist];
        if (is_danger[node]) temp_visited[node] = true;
        cur_exist = fast_parent_idx[cur_exist];
    }

    bool is_subset = true;
    int cur_new = p_idx;
    while (cur_new != -1) {
        int node = fast_node_id[cur_new];
        if (is_danger[node] && !temp_visited[node]) {
            is_subset = false;
            break;
        }
        cur_new = fast_parent_idx[cur_new];
    }
    if (is_danger[v] && !temp_visited[v]) is_subset = false;

    cur_exist = exist_idx;
    while (cur_exist != -1) {
        int node = fast_node_id[cur_exist];
        if (is_danger[node]) temp_visited[node] = false;
        cur_exist = fast_parent_idx[cur_exist];
    }

    return is_subset;
}

bool DiagnosticBiObjectiveBellmanFord::is_pruned(int node_idx, const Cost& current_cost) {
    Cost h = heuristics_bwd[node_idx];
    if (h.c1 >= INF_COST || h.c2 >= INF_COST) return true;

    if (h.c1 <= NEG_INF_COST || h.c2 <= NEG_INF_COST) return false;

    Cost est = current_cost + h;
    if (est.c1 > nadir_s_t.c1 && est.c2 > nadir_s_t.c2) return true;

    const auto& target_labels = L_indices[target_node_id];
    if (target_labels.empty()) return false;

    auto comp = [&](int idx, const Cost& val) { return label_pool[idx].cost.c1 < val.c1; };
    auto it = std::lower_bound(target_labels.begin(), target_labels.end(), est, comp);

    if (it != target_labels.begin()) {
        auto prev_it = std::prev(it);
        while (prev_it != target_labels.begin() && !label_pool[*prev_it].is_valid) --prev_it;
        if (label_pool[*prev_it].is_valid) {
            long long min_c1_seen = label_pool[*prev_it].cost.c1;
            Cost search_cost = { min_c1_seen, 0 };
            auto first_of_group_it = std::lower_bound(target_labels.begin(), it, search_cost, comp);

            while (first_of_group_it != it) {
                if (label_pool[*first_of_group_it].is_valid) {
                    if (label_pool[*first_of_group_it].cost.strictly_dominates(est)) return true;
                    break;
                }
                ++first_of_group_it;
            }
        }
    }
    return false;
}

int DiagnosticBiObjectiveBellmanFord::UpdateWeakParetoSet(int v, const Cost& new_cost, int p_idx) {
    auto& indices = L_indices[v];
    auto comp = [&](int idx, const Cost& val) {
        if (label_pool[idx].cost.c1 != val.c1) return label_pool[idx].cost.c1 < val.c1;
        return label_pool[idx].cost.c2 < val.c2;
    };
    auto it = std::lower_bound(indices.begin(), indices.end(), new_cost, comp);

    auto dup_it = it;
    while (dup_it != indices.end()) {
        Cost exist_cost = label_pool[*dup_it].cost;
        if (!(exist_cost == new_cost)) break;
        if (label_pool[*dup_it].parent_index == p_idx) return -1;
        dup_it++;
    }

    if (it != indices.begin()) {
        auto prev_it = std::prev(it);
        while (true) {
            if (label_pool[*prev_it].cost.strictly_dominates(new_cost)) {
                if (can_exist_dominate_new(*prev_it, p_idx, v)) {
                    return -1;
                }
            }
            if (prev_it == indices.begin()) break;
            --prev_it;
        }
    }

    auto new_end = std::remove_if(it, indices.end(), [&](int idx) {
        if (new_cost.strictly_dominates(label_pool[idx].cost)) {
            if (can_new_dominate_exist(p_idx, v, idx)) {
                label_pool[idx].is_valid = false;
                release_ref(idx);
                return true;
            }
        }
        return false;
        });
    indices.erase(new_end, indices.end());

    it = std::lower_bound(indices.begin(), indices.end(), new_cost, comp);

    int new_idx = create_label(new_cost, v, p_idx);
    indices.insert(it, new_idx);
    return new_idx;
}

bool DiagnosticBiObjectiveBellmanFord::CheckAndRecordCycle(int end_idx, int v, const Edge& edge) {
    int cur = end_idx;
    bool found = false;

    bool is_non_positive_edge = (edge.cost.c1 <= 0 || edge.cost.c2 <= 0);

    while (cur != -1) {
        int u = fast_node_id[cur];

        if (!is_danger[u] && !is_non_positive_edge) {
            break;
        }
        if (u == v) {
            found = true;
            break;
        }
        cur = fast_parent_idx[cur];
    }

    if (!found) return false;

    std::vector<int> path;
    path.push_back(v);
    int trace = end_idx;
    while (trace != -1) {
        path.push_back(fast_node_id[trace]);
        if (fast_node_id[trace] == v) break;
        trace = fast_parent_idx[trace];
    }
    std::reverse(path.begin(), path.end());

    Cost end_cost = label_pool[end_idx].cost + edge.cost;
    Cost start_cost = label_pool[cur].cost;
    Cost delta = { end_cost.c1 - start_cost.c1, end_cost.c2 - start_cost.c2 };
    std::string type = classify_cycle(delta);

    if (type != "Positive Cycle") problematic_cycles.insert({ path, delta, type });
    return true;
}

std::string DiagnosticBiObjectiveBellmanFord::classify_cycle(const Cost& d) {
    bool c1_pos = d.c1 > 0; bool c1_neg = d.c1 < 0; bool c1_zero = d.c1 == 0;
    bool c2_pos = d.c2 > 0; bool c2_neg = d.c2 < 0; bool c2_zero = d.c2 == 0;

    if (c1_pos && c2_pos) return "Positive Cycle";
    if (c1_zero && c2_pos) return "c2 Positive Cycle";
    if (c1_pos && c2_zero) return "c1 Positive Cycle";
    if (c1_neg && c2_pos) return "Trade-off Cycle (-c1, +c2)";
    if (c1_pos && c2_neg) return "Trade-off Cycle (+c1, -c2)";
    if (c1_neg && c2_zero) return "c1 Negative Cycle";
    if (c1_zero && c2_neg) return "c2 Negative Cycle";
    if (c1_neg && c2_neg) return "Strictly Negative Cycle";
    if (c1_zero && c2_zero) return "Zero Cycle";
    return "Unknown Cycle";
}

bool DiagnosticBiObjectiveBellmanFord::is_better(const Cost& val, const Cost& cur, bool min_c1) {
    if (val.c1 <= NEG_INF_COST) {
        if (cur.c1 <= NEG_INF_COST) return false;
        return true;
    }
    if (val.c1 >= INF_COST) return false;
    if (cur.c1 >= INF_COST) return true;
    return min_c1 ? (val.c1 < cur.c1 || (val.c1 == cur.c1 && val.c2 < cur.c2))
        : (val.c2 < cur.c2 || (val.c2 == cur.c2 && val.c1 < cur.c1));
}

std::vector<bool> DiagnosticBiObjectiveBellmanFord::run_forward_spfa_for_danger(int s) {
    std::vector<bool> is_danger_flags(n, false);
    std::vector<Cost> dist_c1(n, { INF_COST, INF_COST });
    std::vector<Cost> dist_c2(n, { INF_COST, INF_COST });

    auto run_spfa = [&](bool min_c1, std::vector<Cost>& dists) {
        dists[s] = { 0, 0 };
        std::deque<int> Q;
        std::vector<bool> in_queue(n, false);
        std::vector<int> update_count(n, 0);
        std::vector<int> parent(n, -1);

        Q.push_back(s);
        in_queue[s] = true;

        while (!Q.empty()) {
            int u = Q.front();
            Q.pop_front();
            in_queue[u] = false;

            for (const auto& e : adj[u]) {
                int v = e.to;

                if (dists[u].c1 <= NEG_INF_COST) {
                    if (dists[v].c1 > NEG_INF_COST) {
                        dists[v] = { NEG_INF_COST, NEG_INF_COST };
                        is_danger_flags[v] = true;
                        if (!in_queue[v]) { Q.push_front(v); in_queue[v] = true; }
                    }
                    continue;
                }

                Cost next = dists[u] + e.cost;

                if (is_better(next, dists[v], min_c1)) {
                    dists[v] = next;
                    parent[v] = u;
                    update_count[v]++;

                    if (update_count[v] >= 4 && (update_count[v] & (update_count[v] - 1)) == 0) {
                        int curr = u;
                        bool in_cycle = false;
                        int steps = 0;
                        while (curr != -1 && steps <= n) {
                            if (curr == v) { in_cycle = true; break; }
                            curr = parent[curr];
                            steps++;
                        }
                        if (in_cycle) {
                            dists[v] = { NEG_INF_COST, NEG_INF_COST };
                            is_danger_flags[v] = true;
                            if (!in_queue[v]) { Q.push_front(v); in_queue[v] = true; }
                            continue;
                        }
                    }

                    if (update_count[v] > n) {
                        dists[v] = { NEG_INF_COST, NEG_INF_COST };
                        is_danger_flags[v] = true;
                        if (!in_queue[v]) { Q.push_front(v); in_queue[v] = true; }
                        continue;
                    }

                    if (!in_queue[v]) {
                        if (!Q.empty() && (min_c1 ? next.c1 < dists[Q.front()].c1 : next.c2 < dists[Q.front()].c2)) {
                            Q.push_front(v);
                        }
                        else {
                            Q.push_back(v);
                        }
                        in_queue[v] = true;
                    }
                }
            }
        }
    };

    run_spfa(true, dist_c1);
    run_spfa(false, dist_c2);
    return is_danger_flags;
}

DiagnosticBiObjectiveBellmanFord::HeuristicResult DiagnosticBiObjectiveBellmanFord::run_backward_spfa(int t) {
    HeuristicResult res;
    res.dist_c1.assign(n, { INF_COST, INF_COST });
    res.dist_c2.assign(n, { INF_COST, INF_COST });

    auto run_spfa = [&](bool min_c1, std::vector<Cost>& dists) {
        dists[t] = { 0, 0 };
        std::deque<int> Q;
        std::vector<bool> in_queue(n, false);
        std::vector<int> parent(n, -1);
        std::vector<int> update_count(n, 0);

        Q.push_back(t);
        in_queue[t] = true;

        while (!Q.empty()) {
            int u = Q.front();
            Q.pop_front();
            in_queue[u] = false;

            for (const auto& e : rev_adj[u]) {
                int v = e.from;

                if (update_count[v] > n) continue;

                Cost next = dists[u] + e.cost;

                if (is_better(next, dists[v], min_c1)) {

                    if (is_danger[v]) {
                        int curr = u;
                        bool in_cycle = false;
                        int steps = 0;
                        while (curr != -1 && steps <= n) {
                            if (curr == v) { in_cycle = true; break; }
                            if (!is_danger[curr]) break;
                            curr = parent[curr];
                            steps++;
                        }
                        if (in_cycle) continue;
                    }

                    dists[v] = next;
                    parent[v] = u;

                    if (!in_queue[v]) {
                        update_count[v]++;

                        if (!Q.empty() && (min_c1 ? next.c1 < dists[Q.front()].c1 : next.c2 < dists[Q.front()].c2)) {
                            Q.push_front(v);
                        }
                        else {
                            Q.push_back(v);
                        }
                        in_queue[v] = true;
                    }
                }
            }
        }
    };

    run_spfa(true, res.dist_c1);
    run_spfa(false, res.dist_c2);
    return res;
}

/*
// 事前探索 |V|-1 ver
// 注意: こちらを使用する場合は WPLC_Solver.hpp の定義を run_backward_spfa(int s, int t) に変更してください．
DiagnosticBiObjectiveBellmanFord::HeuristicResult DiagnosticBiObjectiveBellmanFord::run_backward_spfa(int s, int t) {
    HeuristicResult res;
    res.dist_c1.assign(n, { INF_COST, INF_COST });
    res.dist_c2.assign(n, { INF_COST, INF_COST });
    res.first_c1.assign(n, { INF_COST, INF_COST });
    res.first_c2.assign(n, { INF_COST, INF_COST });

    auto run_hop_bounded_bf = [&](bool min_c1, std::vector<Cost>& dists, std::vector<Cost>& first_dists) {
        dists[t] = { 0, 0 };
        first_dists[t] = { 0, 0 };

        std::vector<Cost> curr_dists = dists;
        std::vector<int> active_nodes;
        active_nodes.push_back(t);

        for (int hop = 1; hop <= n - 1; ++hop) {
            std::vector<Cost> next_dists = curr_dists;
            std::vector<int> next_active_nodes;
            std::vector<bool> next_is_active(n, false);
            std::vector<bool> newly_reached(n, false);

            bool changed = false;

            for (int u : active_nodes) {
                for (const auto& e : rev_adj[u]) {
                    int v = e.from;
                    Cost next_c = curr_dists[u] + e.cost;

                    if (first_dists[v].c1 >= INF_COST) {
                        first_dists[v] = next_c;
                        newly_reached[v] = true;
                    }
                    else if (newly_reached[v]) {
                        if (is_better(next_c, first_dists[v], min_c1)) {
                            first_dists[v] = next_c;
                        }
                    }

                    if (is_better(next_c, next_dists[v], min_c1)) {
                        next_dists[v] = next_c;
                        changed = true;
                        if (!next_is_active[v]) {
                            next_is_active[v] = true;
                            next_active_nodes.push_back(v);
                        }
                    }
                }
            }

            curr_dists = next_dists;
            active_nodes = next_active_nodes;

            if (!changed) break;
        }
        dists = curr_dists;
    };

    run_hop_bounded_bf(true, res.dist_c1, res.first_c1);
    run_hop_bounded_bf(false, res.dist_c2, res.first_c2);
    return res;
}
*/

void DiagnosticBiObjectiveBellmanFord::calculate_heuristics_and_danger_zones(int s, int t) {
    std::cout << "forward spfa (Identifying Danger Zones)...\n";
    is_danger = run_forward_spfa_for_danger(s);
    int danger_count = std::count(is_danger.begin(), is_danger.end(), true);
    std::cout << "[Info] Identified " << danger_count << " danger nodes from start.\n";

    std::cout << "backward spfa (Calculating Heuristics)...\n";
    auto bwd_res = run_backward_spfa(t);
    heuristics_bwd.resize(n);
    for (int i = 0; i < n; ++i) {
        heuristics_bwd[i] = { bwd_res.dist_c1[i].c1, bwd_res.dist_c2[i].c2 };
    }

    Cost solA = bwd_res.dist_c1[s];
    Cost solB = bwd_res.dist_c2[s];

    std::cout << "A(min c1) : (" << solA.c1 << ", " << solA.c2 << ")\n";
    std::cout << "B(min c2) : (" << solB.c1 << ", " << solB.c2 << ")\n";

    if (solA.c1 >= INF_COST || solB.c1 >= INF_COST) {
        nadir_s_t = { INF_COST, INF_COST };
    }
    else {
        nadir_s_t = { std::max(solA.c1, solB.c1), std::max(solA.c2, solB.c2) };
        std::cout << "Nadir Point Estimate: (" << nadir_s_t.c1 << ", " << nadir_s_t.c2 << ")\n";
    }
}


// 理想点・最悪点 |V|-1 ver
// 注意: WPLC_Solver.hpp 側の関数シグネチャの変更が必要です．
/*
void DiagnosticBiObjectiveBellmanFord::calculate_heuristics_and_danger_zones(int s, int t) {
    std::cout << "forward spfa (Identifying Danger Zones)...\n";
    is_danger = run_forward_spfa_for_danger(s);
    int danger_count = std::count(is_danger.begin(), is_danger.end(), true);
    std::cout << "[Info] Identified " << danger_count << " danger nodes from start.\n";

    std::cout << "backward spfa (Calculating Heuristics - Hop Bounded BF)...\n";
    auto bwd_res = run_backward_spfa(s, t);
    heuristics_bwd.resize(n);
    for (int i = 0; i < n; ++i) {
        heuristics_bwd[i] = { bwd_res.dist_c1[i].c1, bwd_res.dist_c2[i].c2 };
    }

    Cost solA = bwd_res.first_c1[s];
    Cost solB = bwd_res.first_c2[s];

    std::cout << "A(min c1) : (" << solA.c1 << ", " << solA.c2 << ") [using first reached simple path]\n";
    std::cout << "B(min c2) : (" << solB.c1 << ", " << solB.c2 << ") [using first reached simple path]\n";

    if (solA.c1 >= INF_COST || solB.c1 >= INF_COST) {
        nadir_s_t = { INF_COST, INF_COST };
    }
    else {
        nadir_s_t = { std::max(solA.c1, solB.c1), std::max(solA.c2, solB.c2) };
        std::cout << "Nadir Point Estimate: (" << nadir_s_t.c1 << ", " << nadir_s_t.c2 << ")\n";
    }
}
*/
