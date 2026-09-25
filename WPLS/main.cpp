#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <random>
#include "WPLS_Solver.hpp"
#include "GraphGenerators.hpp"

int main() {
    int height = 4;
    int width = 4;
    int n = height * width + 2;
    //int start = n - 2; //grid
    //int target = n - 1; //grid

    n = 264346; // DIMACS NY
    //n = 435666; // DIMACS COL
    //n = 1070376; // DIMACS FLA

    std::random_device rd;
    std::mt19937 gen(123);
    std::uniform_int_distribution<int> distrib(0, n - 1);

    int num_trials = 5;
    std::vector<long long> execution_times;
    long long total_time = 0;

    std::cout << "\n========================================\n";
    std::cout << "Starting " << num_trials << " random trials..." << std::endl;
    std::cout << "========================================\n";

    std::string csv_filename = "experiment_results.csv";
    std::ofstream csv_file(csv_filename);
    if (!csv_file.is_open()) {
        std::cerr << "Error: Cannot open " << csv_filename << " for writing.\n";
        return 1;
    }
    csv_file << "start,target,count,exec_ms,Q_max\n";

    for (int i = 0; i < num_trials; ++i) {

        WBDijkstra solver(n);

        //setup_grid_graph_paper(solver, height, width, start, target, distrib(gen));

        int start = distrib(gen);
        int target = distrib(gen);

        // 始点と終点が同じにならないように再抽選
        while (start == target) {
            target = distrib(gen);
        }
        setup_dimacs_graph(solver, "USA-road-d.NY.gr", "USA-road-t.NY.gr");
        //setup_dimacs_graph(solver, "USA-road-d.COL.gr", "USA-road-t.COL.gr");
        //setup_dimacs_graph(solver, "USA-road-d.FLA.gr", "USA-road-t.FLA.gr");


        std::cout << "\n[ Trial " << i + 1 << " / " << num_trials << " ]\n";
        std::cout << "Start searching from " << start << " to " << target << "\n";

        auto start_time = std::chrono::high_resolution_clock::now();

        // BDijkstra::solve は弱パレート解（確定ラベル）の一覧を返す
        std::vector<Label> results = solver.solve(start, target);

        auto end_time = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        long long exec_ms = duration.count();
        execution_times.push_back(exec_ms);
        total_time += exec_ms;
        std::cout << "Execution time: " << exec_ms << " ms\n";

        if (results.empty()) {
            std::cout << "No path found.\n";
        }
        else {
            std::cout << "Found " << results.size() << " weak Pareto solutions.\n";
        }

        csv_file << start << "," << target << "," << results.size() << "," << exec_ms << std::endl;
    }

    csv_file.close();

    std::cout << "\n========================================\n";
    std::cout << "          Summary of " << num_trials << " Trials\n";
    std::cout << "========================================\n";
    for (int i = 0; i < num_trials; ++i) {
        std::cout << "Trial " << i + 1 << ": " << execution_times[i] << " ms\n";
    }
    std::cout << "----------------------------------------\n";
    std::cout << "Average Execution time: " << (total_time / num_trials) << " ms\n";
    std::cout << "========================================\n";

    return 0;
}