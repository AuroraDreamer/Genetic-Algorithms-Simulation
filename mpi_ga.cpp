#include <mpi.h>
#include <omp.h>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Widget.H>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <string>

// Define the Genetic Algorithm Parameters
struct GAParams {
    int population_size;
    int generations;
    double crossover_rate;
    double mutation_rate;
    int num_genes;
    int num_threads;
    int migration_frequency;
};

// Fitness Function (Example: Sum of genes)
double fitnessFunction(const std::vector<double>& individual) {
    return std::accumulate(individual.begin(), individual.end(), 0.0);
}

// Random Double Generator
double randomDouble(double min, double max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

// Standard Genetic Algorithm Implementation
void standardGA(const GAParams& params, int rank, int size) {
    std::vector<std::vector<double>> population(params.population_size, std::vector<double>(params.num_genes));

    // Initialize population
    for (auto& individual : population) {
        for (auto& gene : individual) {
            gene = randomDouble(0.0, 1.0); // Random gene initialization
        }
    }

    // Run for generations
    for (int gen = 0; gen < params.generations; ++gen) {
        // Evaluate fitness
        std::vector<double> fitness_values;
        for (const auto& individual : population) {
            fitness_values.push_back(fitnessFunction(individual));
        }

        // Selection, crossover, and mutation would go here (omitted for simplicity)
        std::cout << "Standard GA executed on rank " << rank << " at generation " << gen << "\n";
    }
}

// Island Model Genetic Algorithm Implementation
void islandGA(const GAParams& params, int rank, int size) {
    std::vector<std::vector<double>> population(params.population_size, std::vector<double>(params.num_genes));

    // Initialize population
    for (auto& individual : population) {
        for (auto& gene : individual) {
            gene = randomDouble(0.0, 1.0);
        }
    }

    // Run for generations
    for (int gen = 0; gen < params.generations; ++gen) {
        // Evaluate fitness
        std::vector<double> fitness_values;
        for (const auto& individual : population) {
            fitness_values.push_back(fitnessFunction(individual));
        }

        // Selection, crossover, mutation, and migration between islands (simplified)
        std::cout << "Island GA executed on rank " << rank << " at generation " << gen << "\n";
    }
}

// Stationary Genetic Algorithm Implementation
void stationaryGA(const GAParams& params, int rank, int size) {
    std::vector<std::vector<double>> population(params.population_size, std::vector<double>(params.num_genes));

    // Initialize population
    for (auto& individual : population) {
        for (auto& gene : individual) {
            gene = randomDouble(0.0, 1.0);
        }
    }

    // Run for generations
    for (int gen = 0; gen < params.generations; ++gen) {
        // Evaluate fitness
        std::vector<double> fitness_values;
        for (const auto& individual : population) {
            fitness_values.push_back(fitnessFunction(individual));
        }

        // Replace a portion of the population (stationary model)
        std::cout << "Stationary GA executed on rank " << rank << " at generation " << gen << "\n";
    }
}

// Elitist Genetic Algorithm Implementation
void elitistGA(const GAParams& params, int rank, int size) {
    std::vector<std::vector<double>> population(params.population_size, std::vector<double>(params.num_genes));

    // Initialize population
    for (auto& individual : population) {
        for (auto& gene : individual) {
            gene = randomDouble(0.0, 1.0);
        }
    }

    // Run for generations
    for (int gen = 0; gen < params.generations; ++gen) {
        // Evaluate fitness
        std::vector<double> fitness_values;
        for (const auto& individual : population) {
            fitness_values.push_back(fitnessFunction(individual));
        }

        // Preserve best individuals (elitist strategy)
        std::cout << "Elitist GA executed on rank " << rank << " at generation " << gen << "\n";
    }
}

// Plot Speedup and Efficiency using Gnuplot
void plotGraphs(const std::vector<int>& threads, const std::vector<double>& speedups, const std::vector<double>& efficiencies) {
    std::ofstream file("performance_data.dat");
    for (size_t i = 0; i < threads.size(); ++i) {
        file << threads[i] << " " << speedups[i] << " " << efficiencies[i] << "\n";
    }
    file.close();

    std::ofstream gnuplot_script("plot_speedup_efficiency.gnuplot");
    gnuplot_script << "set title 'Performance Metrics'\n";
    gnuplot_script << "set xlabel 'Number of Threads'\n";
    gnuplot_script << "set ylabel 'Speedup / Efficiency'\n";
    gnuplot_script << "set key left top\n";
    gnuplot_script << "plot 'performance_data.dat' using 1:2 title 'Speedup' with lines, \\\n";
    gnuplot_script << "     'performance_data.dat' using 1:3 title 'Efficiency' with lines\n";
    gnuplot_script.close();

    system("gnuplot -p plot_speedup_efficiency.gnuplot");
}

// FLTK GUI Integration
void runGA(Fl_Widget*, void* data) {
    GAParams* params = static_cast<GAParams*>(data);

    // MPI Initialization
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Execute Selected Algorithm
    if (params->num_genes == 0) {
        standardGA(*params, rank, size);
    } else if (params->num_genes == 1) {
        islandGA(*params, rank, size);
    } else if (params->num_genes == 2) {
        stationaryGA(*params, rank, size);
    } else if (params->num_genes == 3) {
        elitistGA(*params, rank, size);
    }

    // Finalize MPI
    MPI_Finalize();
}

int main(int argc, char** argv) {
    Fl_Window* window = new Fl_Window(400, 300, "Genetic Algorithms with MPI");
    GAParams params;

    Fl_Input* pop_size_input = new Fl_Input(150, 50, 200, 25, "Population Size:");
    Fl_Input* gen_input = new Fl_Input(150, 80, 200, 25, "Generations:");
    Fl_Input* threads_input = new Fl_Input(150, 110, 200, 25, "Number of Threads:");
    Fl_Input* crossover_input = new Fl_Input(150, 140, 200, 25, "Crossover Rate:");
    Fl_Input* mutation_input = new Fl_Input(150, 170, 200, 25, "Mutation Rate:");
    Fl_Choice* algo_choice = new Fl_Choice(150, 200, 200, 25, "Algorithm:");

    algo_choice->add("Standard");
    algo_choice->add("Island Model");
    algo_choice->add("Stationary");
    algo_choice->add("Elitist");

    // Set user data to pass parameters to callback
    Fl_Button* run_button = new Fl_Button(150, 250, 100, 30, "Run");
    run_button->user_data((void*)&params); // Set user data

    run_button->callback(runGA, &params); // Link callback function

    window->end();
    window->show(argc, argv);

    return Fl::run();
}

