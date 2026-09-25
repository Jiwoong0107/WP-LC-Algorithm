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
    double c1 = INF, c2 = INF;

    bool operator<(const Cost& other) const {
        if (std::abs(c1 - other.c1) > 1e-9) {
            return c1 < other.c1;
        }
        return c2 < other.c2;
    }

    bool operator==(const Cost& other) const {
        return std::abs(c1 - other.c1) < 1e-9 && std::abs(c2 - other.c2) < 1e-9;
    }

    bool dominates(const Cost& other) const {
        return  (c1 <= other.c1 && c2 <= other.c2) && !(*this == other);
    }

    // 弱パレート解用の「強い支配」判定
    bool strong_dominates(const Cost& other) const {
        return  (c1 < other.c1&& c2 < other.c2);
    }

    Cost operator+(const Cost& other) const {
        return { c1 + other.c1, c2 + other.c2 };
    }
};


struct Arc {
    int to;
    Cost cost;
};

struct ReverseArc {
    int from;
    Cost cost;
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