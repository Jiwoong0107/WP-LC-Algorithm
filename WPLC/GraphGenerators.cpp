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
