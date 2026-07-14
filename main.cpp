#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <random>
#include "WPLC_Solver.hpp"
#include "GraphGenerators.hpp"

int main() {
    int height = 4;
    int width = 4;
    int n = height * width + 2;
    int start = n - 2; //grid
    int target = n - 1; //grid
    //int a = 3; n = a * a;
    //n = 264346; // DIMACS NY
    //n = 435666; // DIMACS COL
    //n = 1070376; // DIMACS FLA

    std::random_device rd;
    std::mt19937 gen(125);
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

        DiagnosticBiObjectiveBellmanFord solver(n);

        setup_grid_graph_paper(solver, height, width, start, target, distrib(gen));
        
        //setup_dimacs_graph(solver, "USA-road-d.NY.gr", "USA-road-t.NY.gr");
        //setup_dimacs_graph(solver, "USA-road-d.COL.gr", "USA-road-t.COL.gr");
        //setup_dimacs_graph(solver, "USA-road-d.FLA.gr", "USA-road-t.FLA.gr");

        //setup_dimacs_graph_with_elevation(solver, "USA-road-d.NY.gr", "USA-road-t.NY.gr", n);
        //setup_dimacs_graph_with_elevation(solver, "USA-road-d.COL.gr", "USA-road-t.COL.gr", n);


        //solver.print_graph();

        max_queue_size = 0;

        std::cout << "\n[ Trial " << i + 1 << " / " << num_trials << " ]\n";
        std::cout << "Start searching from " << start << " to " << target << "\n";

        auto start_time = std::chrono::high_resolution_clock::now();

        bool is_completed = solver.solve(start, target);

        auto end_time = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        long long exec_ms = duration.count();
        execution_times.push_back(exec_ms);
        total_time += exec_ms;
        std::cout << "Qmax : " << max_queue_size << "\n";
        std::cout << "Execution time: " << exec_ms << " ms\n";

        if (!is_completed) {
            std::cout << "Status: TIMEOUT (Over 60 minutes)\n";
        }

        int solution_count = is_completed ? solver.get_weak_pareto_count(target) : -1;

        csv_file << start << "," << target << "," << solution_count << "," << exec_ms << "," << max_queue_size << std::endl;
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