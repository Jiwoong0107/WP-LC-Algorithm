#include "GraphGenerators.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

void setup_dimacs_graph(WBDijkstra& bd, const std::string& filenameD, const std::string& filenameT) {
    std::cout << "Loading graph from " << filenameD << " (as Cost1) & " << filenameT << " (as Cost2)..." << std::endl;

    std::ifstream fileD(filenameD);
    if (!fileD.is_open()) {
        std::cerr << "Error: Cannot open file " << filenameD << std::endl;
        return;
    }

    std::ifstream fileT(filenameT);
    if (!fileT.is_open()) {
        std::cerr << "Error: Cannot open file " << filenameT << std::endl;
        fileD.close();
        return;
    }

    std::string lineD, lineT;
    long long edge_count = 0;
    bool p_line_processed = false;

    while (std::getline(fileD, lineD) && std::getline(fileT, lineT)) {
        while (!lineD.empty() && (lineD[0] == 'c' || lineD[0] == 'p')) {
            if (lineD[0] == 'p' && !p_line_processed) {
                std::stringstream ssD(lineD);
                char typeD;
                int num_nodes, num_arcs;
                std::string problem_type;
                ssD >> typeD >> problem_type >> num_nodes >> num_arcs;
                std::cout << "Graph info (from " << filenameD << "): "
                    << num_nodes << " nodes, " << num_arcs << " arcs." << std::endl;
                p_line_processed = true;
            }
            if (!std::getline(fileD, lineD)) break;
        }

        while (!lineT.empty() && (lineT[0] == 'c' || lineT[0] == 'p')) {
            if (!std::getline(fileT, lineT)) break;
        }

        if (lineD.empty() || lineT.empty()) {
            if (lineD.empty() == lineT.empty()) {
                break;
            }
            else {
                std::cerr << "Warning: Files have different line counts." << std::endl;
                break;
            }
        }

        std::stringstream ssD(lineD);
        char typeD;
        int uD, vD, weightD;
        ssD >> typeD >> uD >> vD >> weightD;

        std::stringstream ssT(lineT);
        char typeT;
        int uT, vT, weightT;
        ssT >> typeT >> uT >> vT >> weightT;

        if (typeD == 'a' && typeT == 'a') {
            if (uD != uT || vD != vT) {
                std::cerr << "Error: Graph structure mismatch!" << std::endl;
                std::cerr << "  FileD: a " << uD << " " << vD << std::endl;
                std::cerr << "  FileT: a " << uT << " " << vT << std::endl;
                break;
            }

            int u_idx = uD - 1;
            int v_idx = vD - 1;

            bd.add_edge(u_idx, v_idx, (double)weightD, (double)weightT);
            edge_count++;
        }
        else if (typeD != 'a') {
            std::cerr << "Error: Expected 'a' line in " << filenameD << " but got: " << lineD << std::endl;
            break;
        }
        else {
            std::cerr << "Error: Expected 'a' line in " << filenameT << " but got: " << lineT << std::endl;
            break;
        }
    }
    fileD.close();
    fileT.close();
    std::cout << "Graph loading finish. Added " << edge_count << " edges." << std::endl;
}

void setup_grid_graph_paper(WBDijkstra& bd, int height, int width, int s, int t, int i) {
    std::cout << "making grid (paper specification)\n";
    std::mt19937 rng(123 + i); // 再現性のためにシードを固定

    // BDijkstra は Dijkstra ベースで負の辺を扱えないため、1〜10 の正の範囲でランダムに生成する
    std::uniform_int_distribution<int> dist(1, 10);

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