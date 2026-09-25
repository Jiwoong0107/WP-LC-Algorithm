#ifndef GRAPH_TYPES_HPP
#define GRAPH_TYPES_HPP

#include <vector>
#include <string>
#include <algorithm>
#include <limits>

constexpr double INF = std::numeric_limits<double>::max();
constexpr double NEG_INF = -std::numeric_limits<double>::max();

constexpr long long INF_COST = 4000000000000000000LL;
constexpr long long NEG_INF_COST = -4000000000000000000LL;

struct Cost {
    long long c1 = 0, c2 = 0;

    bool operator<(const Cost& other) const {
        if (c1 != other.c1) return c1 < other.c1;
        return c2 < other.c2;
    }
    bool operator==(const Cost& other) const {
        return c1 == other.c1 && c2 == other.c2;
    }
    Cost operator+(const Cost& other) const {
        if (c1 <= NEG_INF_COST || other.c1 <= NEG_INF_COST) return { NEG_INF_COST, NEG_INF_COST };
        if (c2 <= NEG_INF_COST || other.c2 <= NEG_INF_COST) return { NEG_INF_COST, NEG_INF_COST };
        if (c1 >= INF_COST || other.c1 >= INF_COST || c2 >= INF_COST || other.c2 >= INF_COST) return { INF_COST, INF_COST };
        return { c1 + other.c1, c2 + other.c2 };
    }
    bool strictly_dominates(const Cost& other) const {
        return (c1 < other.c1&& c2 < other.c2);
    }
};

struct Edge {
    int to;
    Cost cost;
};

struct RevEdge {
    int from;
    Cost cost;
};

struct Label {
    Cost cost;
    int node_id;
    int parent_index;
    int ref_count;
    bool is_valid;

    Label() : cost({ 0,0 }), node_id(-1), parent_index(-1), ref_count(0), is_valid(false) {}
    Label(Cost c, int id, int p_idx)
        : cost(c), node_id(id), parent_index(p_idx), ref_count(1), is_valid(true) {}
};

struct DetectedCycle {
    std::vector<int> path;
    std::vector<int> canonical_path;
    Cost total_cost;
    std::string type;

    DetectedCycle(std::vector<int> p, Cost c, std::string t)
        : path(p), total_cost(c), type(t) {
        std::vector<int> temp = p;
        if (!temp.empty() && temp.front() == temp.back()) temp.pop_back();
        if (!temp.empty()) {
            auto min_it = std::min_element(temp.begin(), temp.end());
            std::rotate(temp.begin(), min_it, temp.end());
        }
        if (!temp.empty()) temp.push_back(temp.front());
        canonical_path = temp;
    }

    bool operator<(const DetectedCycle& other) const {
        if (total_cost.c1 != other.total_cost.c1) return total_cost.c1 < other.total_cost.c1;
        if (total_cost.c2 != other.total_cost.c2) return total_cost.c2 < other.total_cost.c2;
        return canonical_path < other.canonical_path;
    }
};

#endif // GRAPH_TYPES_HPP