#include <mpi.h>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Text_Display.H>
#include <FL/fl_draw.H>
#include <thread>
#include <atomic>
#include <vector>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <string>

using namespace std;

// Custom Line Chart Class (same as before)
class LineChart : public Fl_Widget {
    vector<pair<double, double>> data;
    string xLabel, yLabel, title;

public:
    LineChart(int x, int y, int w, int h, const char *l = 0)
        : Fl_Widget(x, y, w, h, l), xLabel("X"), yLabel("Y"), title("Line Chart") {}

    void setLabels(const string &x, const string &y, const string &t) {
        xLabel = x;
        yLabel = y;
        title = t;
    }

    void addDataPoint(double x, double y) {
        data.emplace_back(x, y);
        if (data.size() > 20) data.erase(data.begin()); // Limit data points
        redraw();
    }

    void draw() override {
        fl_color(FL_WHITE);
        fl_rectf(x(), y(), w(), h()); // Clear background

        fl_color(FL_BLACK);
        fl_draw(title.c_str(), x() + w() / 2 - 50, y() + 20); // Title

        fl_color(FL_BLACK);
        fl_line(x() + 50, y() + h() - 50, x() + w() - 20, y() + h() - 50); // X-axis
        fl_line(x() + 50, y() + h() - 50, x() + 50, y() + 20);             // Y-axis

        if (data.empty()) return;

        double xMin = data.front().first, xMax = data.back().first;
        double yMin = 0, yMax = 0;
        for (const auto &[x, y] : data) yMax = max(yMax, y);

        for (size_t i = 0; i < data.size() - 1; ++i) {
            double x1 = mapValue(data[i].first, xMin, xMax, x() + 50, x() + w() - 20);
            double y1 = mapValue(data[i].second, yMin, yMax, y() + h() - 50, y() + 20);
            double x2 = mapValue(data[i + 1].first, xMin, xMax, x() + 50, x() + w() - 20);
            double y2 = mapValue(data[i + 1].second, yMin, yMax, y() + h() - 50, y() + 20);

            fl_color(FL_BLUE);
            fl_line(x1, y1, x2, y2);
        }
    }

private:
    double mapValue(double value, double minVal, double maxVal, double outMin, double outMax) {
        return outMin + ((value - minVal) / (maxVal - minVal)) * (outMax - outMin);
    }
};

// Global Variables
atomic<bool> running(false);
Fl_Input *inputPopulationSize, *inputMutationRate, *inputCrossoverRate;
Fl_Choice *algorithmChoice;
Fl_Text_Display *outputDisplay;
Fl_Text_Buffer *outputBuffer;
LineChart *speedupChart, *efficiencyChart;

// MPI Logic: Algorithm simulation (with MPI)
void runSelectedAlgorithm(int algorithm, int threads, int populationSize, double mutationRate, double crossoverRate, int rank) {
    // Simulate computation delay for each process based on the number of threads (simulate genetic algorithm)
    this_thread::sleep_for(chrono::milliseconds(100 * threads)); // Simulate computation delay for each thread
}

// MPI Simulation Function
void simulateAlgorithm() {
    running = true;

    // Validate inputs
    int populationSize = atoi(inputPopulationSize->value());
    double mutationRate = atof(inputMutationRate->value());
    double crossoverRate = atof(inputCrossoverRate->value());
    if (populationSize <= 0 || mutationRate < 0 || mutationRate > 1 || crossoverRate < 0 || crossoverRate > 1) {
        outputBuffer->append("Invalid input. Please check parameters.\n");
        running = false;
        return;
    }

    int selectedAlgorithm = algorithmChoice->value();
    outputBuffer->append("Simulation started...\n");

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the rank of the current process
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Get the number of processes in MPI_COMM_WORLD

    auto startTime = chrono::high_resolution_clock::now();

    // Loop through different numbers of threads (using MPI rank to distribute work)
    for (int threads = 1; threads <= 16 && running; threads *= 2) {
        // Simulate algorithm execution
        runSelectedAlgorithm(selectedAlgorithm, threads, populationSize, mutationRate, crossoverRate, rank);

        double elapsed = 0.1 * threads; // Simulate time based on the number of threads
        double speedup = 1.0 / elapsed;
        double efficiency = speedup / threads;

        // Add data points for chart updates (gather data using MPI)
        if (rank == 0) {
            speedupChart->addDataPoint(threads, speedup);
            efficiencyChart->addDataPoint(threads, efficiency);
        }

        string result = "Rank: " + to_string(rank) + " Threads: " + to_string(threads) +
                        " Speedup: " + to_string(speedup) +
                        " Efficiency: " + to_string(efficiency) + "\n";

        // Send results to rank 0
        if (rank != 0) {
            MPI_Send(result.c_str(), result.size() + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }

        // Rank 0 collects results and appends to output buffer
        if (rank == 0) {
            outputBuffer->append(result.c_str());
            for (int i = 1; i < size; ++i) {
                char buffer[256];
                MPI_Recv(buffer, 256, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                outputBuffer->append(buffer);
            }
        }
    }

    auto endTime = chrono::high_resolution_clock::now();
    double totalElapsed = chrono::duration<double>(endTime - startTime).count();

    // Final output (only from rank 0)
    if (rank == 0) {
        outputBuffer->append(("Simulation complete in " + to_string(totalElapsed) + " seconds.\n").c_str());
    }

    running = false;
}

// Button Callbacks
void startCallback(Fl_Widget *, void *) {
    if (running) {
        outputBuffer->append("Simulation already running!\n");
        return;
    }
    thread simulationThread(simulateAlgorithm);
    simulationThread.detach();
}

void stopCallback(Fl_Widget *, void *) {
    running = false;
    outputBuffer->append("Simulation stopped.\n");
}

// Main Function
int main(int argc, char **argv) {
    srand(time(0));

    // Initialize MPI
    MPI_Init(&argc, &argv);

    Fl_Window *window = new Fl_Window(900, 600, "MPI Genetic Algorithm Simulator");

    // Inputs
    inputPopulationSize = new Fl_Input(150, 30, 150, 25, "Population Size:");
    inputPopulationSize->value("100");

    inputMutationRate = new Fl_Input(150, 70, 150, 25, "Mutation Rate:");
    inputMutationRate->value("0.01");

    inputCrossoverRate = new Fl_Input(150, 110, 150, 25, "Crossover Rate:");
    inputCrossoverRate->value("0.7");

    algorithmChoice = new Fl_Choice(150, 150, 150, 25, "Algorithm:");
    algorithmChoice->add("Standard GA|Stationary GA|Elite GA|Island Model GA");
    algorithmChoice->value(0);

    // Buttons
    Fl_Button *startButton = new Fl_Button(50, 200, 100, 30, "Start");
    startButton->callback(startCallback);

    Fl_Button *stopButton = new Fl_Button(200, 200, 100, 30, "Stop");
    stopButton->callback(stopCallback);

    // Output Display
    outputDisplay = new Fl_Text_Display(50, 250, 400, 300);
    outputBuffer = new Fl_Text_Buffer();
    outputDisplay->buffer(outputBuffer);

    // Line Charts
    speedupChart = new LineChart(500, 50, 350, 200);
    speedupChart->setLabels("Threads", "Speedup (x)", "Speedup vs. Threads");

    efficiencyChart = new LineChart(500, 300, 350, 200);
    efficiencyChart->setLabels("Threads", "Efficiency (%)", "Efficiency vs. Threads");

    window->end();
    window->show(argc, argv);

    // Run FLTK event loop
    int ret = Fl::run();

    // Finalize MPI
    MPI_Finalize();

    return ret;
}
