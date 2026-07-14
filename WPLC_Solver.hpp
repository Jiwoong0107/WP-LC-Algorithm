#ifndef WPLC_SOLVER_HPP
#define WPLC_SOLVER_HPP

#include "GraphTypes.hpp"
#include <vector>
#include <set>
#include <string>
#include <deque>

extern size_t max_queue_size;

class DiagnosticBiObjectiveBellmanFord {
public:
    DiagnosticBiObjectiveBellmanFord(int num_nodes);

    void add_edge(int u, int v, long long c1, long long c2);
    bool solve(int s, int t);
    void print_diagnosis(int t, bool cyc_result, bool sol_result);
    int get_weak_pareto_count(int t) const;
    void print_graph();

private:
    int n;
    int target_node_id;
    std::vector<std::vector<Edge>> adj;
    std::vector<std::vector<RevEdge>> rev_adj;

    std::vector<Label> label_pool;
    std::vector<int> free_list;
    std::vector<int> fast_node_id;
    std::vector<int> fast_parent_idx;

    std::vector<std::vector<int>> L_indices;
    std::vector<std::vector<int>> U_indices;

    std::set<DetectedCycle> problematic_cycles;

    std::vector<Cost> heuristics_bwd;
    std::vector<bool> is_danger;
    std::vector<bool> temp_visited;

    Cost nadir_s_t;
    long long pruned_count = 0;

    struct HeuristicResult {
        std::vector<Cost> dist_c1, dist_c2;
    };

    int create_label(const Cost& c, int id, int parent_idx);
    void release_ref(int idx);
    bool can_exist_dominate_new(int exist_idx, int p_idx, int v);
    bool can_new_dominate_exist(int p_idx, int v, int exist_idx);
    bool is_pruned(int node_idx, const Cost& current_cost);
    int update_label_set_and_create(int v, const Cost& new_cost, int p_idx);
    bool check_and_record_cycle(int end_idx, int v, const Edge& edge);
    std::string classify_cycle(const Cost& d);
    bool is_better(const Cost& val, const Cost& cur, bool min_c1);
    std::vector<bool> run_forward_spfa_for_danger(int s);
    HeuristicResult run_backward_spfa(int t);
    void calculate_heuristics_and_danger_zones(int s, int t);
};

#endif // WPLC_SOLVER_HPP