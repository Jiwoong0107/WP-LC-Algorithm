#include "GraphGenerators.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <numeric>
#include <algorithm>

void setup_dimacs_graph(DiagnosticBiObjectiveBellmanFord& bd, const std::string& filenameD, const std::string& filenameT) {
    std::cout << "Loading graph from " << filenameD << " & " << filenameT << "..." << std::endl;
    std::ifstream fileD(filenameD);
    std::ifstream fileT(filenameT);
    if (!fileD.is_open() || !fileT.is_open()) {
        std::cerr << "Error: Cannot open DIMACS files." << std::endl;
        return;
    }

    std::string lineD, lineT;
    long long edge_count = 0;

    while (std::getline(fileD, lineD) && std::getline(fileT, lineT)) {
        while (!lineD.empty() && (lineD[0] == 'c' || lineD[0] == 'p')) {
            if (!std::getline(fileD, lineD)) break;
        }
        while (!lineT.empty() && (lineT[0] == 'c' || lineT[0] == 'p')) {
            if (!std::getline(fileT, lineT)) break;
        }
        if (lineD.empty() || lineT.empty()) break;

        std::stringstream ssD(lineD);
        char typeD; int uD, vD, weightD;
        ssD >> typeD >> uD >> vD >> weightD;

        std::stringstream ssT(lineT);
        char typeT; int uT, vT, weightT;
        ssT >> typeT >> uT >> vT >> weightT;

        if (typeD == 'a' && typeT == 'a') {
            bd.add_edge(uD - 1, vD - 1, weightD, weightT);
            edge_count++;
        }
    }
    fileD.close();
    fileT.close();

    std::cout << "Graph loading finish. Added " << edge_count << " edges." << std::endl;
}

void setup_test_graph_with_2_cycles(DiagnosticBiObjectiveBellmanFord& bd) {
    std::cout << "Loading custom test graph with 2 trade-off cycles...\n";

    bd.add_edge(0, 1, 10, 10);
    bd.add_edge(1, 2, 10, 10);
    bd.add_edge(2, 3, 10, 10);
    bd.add_edge(3, 7, 10, 10);

    bd.add_edge(0, 2, 15, 25);
    bd.add_edge(2, 7, 25, 15);

    bd.add_edge(1, 4, 5, 5);
    bd.add_edge(4, 5, -10, 20);
    bd.add_edge(5, 4, 2, -5);

    bd.add_edge(3, 6, 20, -10);
    bd.add_edge(6, 3, -5, 5);
}

void setup_research_based_dense_graph(
    DiagnosticBiObjectiveBellmanFord& bd,
    int n,
    int density_m_per_n,
    int max_chord_length,
    double correlation,
    bool allow_negative
) {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist_01(0.0, 1.0);
    std::uniform_int_distribution<int> node_dist(0, n - 1);

    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), rng);

    long long target_edges = (long long)n * density_m_per_n;
    long long current_edges = 0;

    auto generate_correlated_costs = [&](double base_magnitude) -> std::pair<double, double> {
        double u = dist_01(rng);
        double c1 = u * base_magnitude;

        double noise = (dist_01(rng) - 0.5) * (base_magnitude * 0.1);
        double c2 = base_magnitude * (1.0 - correlation * u + (1.0 + correlation) * (1.0 - u)) + noise;

        c1 = std::max(1.0, c1);
        c2 = std::max(1.0, c2);
        return { c1, c2 };
    };

    std::cout << "Generating Dense Research Instance..." << std::endl;
    std::cout << " - Nodes: " << n << std::endl;
    std::cout << " - Density Target: " << density_m_per_n << "x (" << target_edges << " edges)" << std::endl;
    std::cout << " - Max Chord: " << max_chord_length << " (Controls Diameter)" << std::endl;
    std::cout << " - Correlation: " << correlation << std::endl;

    for (int i = 0; i < n - 1; ++i) {
        auto costs = generate_correlated_costs(100.0);
        bd.add_edge(p[i], p[i + 1], costs.first, costs.second);
        current_edges++;
    }

    for (int u = 0; u < n; ++u) {
        int edges_added_for_node = 0;
        int attempts = density_m_per_n + 10;

        for (int k = 0; k < attempts; ++k) {
            if (current_edges >= target_edges) break;

            int jump = 2 + (k % (max_chord_length - 1));
            int v = u + jump;

            if (v < n) {
                auto costs = generate_correlated_costs(100.0);
                bd.add_edge(p[u], p[v], costs.first, costs.second);
                current_edges++;
            }
            else {
                int back_v = u - jump;
                if (back_v >= 0) {
                    auto costs = generate_correlated_costs(200.0);
                    bd.add_edge(u, back_v, costs.first, costs.second);
                    current_edges++;
                }
            }
        }
    }

    if (allow_negative) {
        int num_adversarial_cycles = std::max(5, n / 1000);
        for (int i = 0; i < num_adversarial_cycles; ++i) {
            int u = node_dist(rng);
            int jump = std::uniform_int_distribution<int>(1, std::min(50, u))(rng);
            int v = u - jump;
            if (v < 0) continue;

            double c1 = -50.0;
            double c2 = 40.0;

            bd.add_edge(u, v, c1, c2);
            current_edges++;
        }
    }
}

void setup_dimacs_graph_with_elevation(DiagnosticBiObjectiveBellmanFord& bd,
    const std::string& fD,
    const std::string& fT,
    int num_nodes) {
    std::ifstream fileD(fD), fileT(fT);
    if (!fileD.is_open() || !fileT.is_open()) {
        std::cerr << "Error: Cannot open files." << std::endl;
        return;
    }
    std::mt19937 rng(100);
    std::uniform_int_distribution<int> elevation_dist(0, 20000);
    std::vector<int> elevations(num_nodes);
    for (int i = 0; i < num_nodes; ++i) elevations[i] = elevation_dist(rng);

    std::string lD, lT;
    while (std::getline(fileD, lD) && std::getline(fileT, lT)) {
        if (lD[0] == 'a') {
            std::stringstream ssD(lD), ssT(lT);
            char t; int u, v, wD, wT_orig;
            ssD >> t >> u >> v >> wD; ssT >> t >> u >> v >> wT_orig;
            long long base_cost = std::max(1LL, (long long)wT_orig);
            long long potential_diff = (long long)elevations[v - 1] - (long long)elevations[u - 1];
            long long c2 = base_cost + potential_diff;
            bd.add_edge(u - 1, v - 1, (long long)wD, c2);
        }
    }
    std::cout << "Graph loaded with Potential-based costs.\n";
}

void setup_missing_solution_graph(DiagnosticBiObjectiveBellmanFord& bd) {
    std::cout << "Loading theoretical limitation trap graph...\n";
    bd.add_edge(0, 1, 5, 5);
    bd.add_edge(1, 2, 5, 5);
    bd.add_edge(0, 2, 20, 20);
    bd.add_edge(2, 1, -50, -50);
    bd.add_edge(1, 3, 0, 0);
}

void setup_random_graph_negative_no_cycle(DiagnosticBiObjectiveBellmanFord& bd, int num_nodes, int num_edges, int base_min, int base_max, int potential_scale) {
    std::cout << "making random graph with negative edges (no negative cycles)\n";

    if (base_min < 1) base_min = 1;

    std::mt19937 rng(123);

    std::vector<int> h1(num_nodes);
    std::vector<int> h2(num_nodes);
    std::uniform_int_distribution<int> pot_dist(-potential_scale, potential_scale);

    for (int i = 0; i < num_nodes; ++i) {
        h1[i] = pot_dist(rng);
        h2[i] = pot_dist(rng);
    }

    std::uniform_int_distribution<int> node_dist(0, num_nodes - 1);
    std::uniform_int_distribution<int> base_dist(base_min, base_max);

    std::set<std::pair<int, int>> existing_edges;

    int edges_count = 0;
    while (edges_count < num_edges) {
        int u = node_dist(rng);
        int v = node_dist(rng);

        if (u == v || existing_edges.count({ u, v })) {
            continue;
        }

        int base_cost1 = base_dist(rng);
        int base_cost2 = base_dist(rng);

        int cost1 = base_cost1 + h1[u] - h1[v];
        int cost2 = base_cost2 + h2[u] - h2[v];

        bd.add_edge(u, v, cost1, cost2);

        existing_edges.insert({ u, v });
        edges_count++;
    }
}

void setup_hidden_cycles_graph(DiagnosticBiObjectiveBellmanFord& bd) {
    std::cout << "Loading graph with hidden cycles (Pruning Test)...\n";

    bd.add_edge(0, 1, 10, 10);
    bd.add_edge(1, 2, -100, -100);
    bd.add_edge(2, 1, -100, -100);

    bd.add_edge(0, 3, 50, 50);
    bd.add_edge(3, 4, 0, 0);
    bd.add_edge(4, 3, 0, 0);
    bd.add_edge(3, 5, 10, 10);

    bd.add_edge(0, 5, 5, 5);
    bd.add_edge(3, 0, -50, -40);
}

void setup_grid_graph1(DiagnosticBiObjectiveBellmanFord& bd, int height, int width, int lower_limit, int upper_limit) {
    std::cout << "making grid\n";
    std::mt19937 rng(123);
    std::uniform_int_distribution<int> dist(lower_limit, upper_limit);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int u = y * width + x;
            if (y + 1 < height) {
                int a = dist(rng); int b = dist(rng);
                int v = (y + 1) * width + x;
                bd.add_edge(u, v, a, b);
                bd.add_edge(v, u, a, b);
            }
            if (x + 1 < width) {
                int c = dist(rng); int d = dist(rng);
                int v = y * width + (x + 1);
                bd.add_edge(u, v, c, d);
                bd.add_edge(v, u, c, d);
            }
        }
    }
    std::cout << "grid finish\n";
}

void setup_grid_graph(DiagnosticBiObjectiveBellmanFord& bd, int width, int height) {
    std::mt19937 rng(123);
    std::uniform_int_distribution<int> dist(-10, 10);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int u = y * width + x;
            if (y + 1 < height) {
                int v = (y + 1) * width + x;
                bd.add_edge(u, v, dist(rng), dist(rng));
                bd.add_edge(v, u, dist(rng), dist(rng));
            }
            if (x + 1 < width) {
                int v = y * width + (x + 1);
                bd.add_edge(u, v, dist(rng), dist(rng));
                bd.add_edge(v, u, dist(rng), dist(rng));
            }
        }
    }
}

void setup_grid_graph_paper(DiagnosticBiObjectiveBellmanFord& bd, int height, int width, int s, int t, int i) {
    std::cout << "making grid (paper specification)\n";
    std::mt19937 rng(123 + i);

    std::uniform_int_distribution<int> dist(-10, 10);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int u = y * width + x;

            if (y + 1 < height) {
                int v = (y + 1) * width + x;
                bd.add_edge(u, v, dist(rng), dist(rng));
                bd.add_edge(v, u, dist(rng), dist(rng));
            }

            if (x + 1 < width) {
                int v = y * width + (x + 1);
                bd.add_edge(u, v, dist(rng), dist(rng));
                bd.add_edge(v, u, dist(rng), dist(rng));
            }
        }
    }

    for (int y = 0; y < height; ++y) {
        int left_node = y * width + 0;
        bd.add_edge(s, left_node, dist(rng), dist(rng));

        int right_node = y * width + (width - 1);
        bd.add_edge(right_node, t, dist(rng), dist(rng));
    }

    std::cout << "grid finish\n";
}

void setup_nwmoa_killer_graph(DiagnosticBiObjectiveBellmanFord& bd) {
    std::cout << "Loading NWMOA* Killer Graph...\n";

    bd.add_edge(0, 1, 10, 10);
    bd.add_edge(1, 3, 10, 10);

    bd.add_edge(1, 5, 1, 0);
    bd.add_edge(5, 1, 0, 0);

    bd.add_edge(0, 2, 100, 100);
    bd.add_edge(2, 3, 100, 100);

    bd.add_edge(2, 4, 10, 10);
    bd.add_edge(4, 2, -10, -10);
}