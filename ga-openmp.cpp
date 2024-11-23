#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Text_Display.H>
#include <FL/fl_draw.H>
#include <string>
#include <atomic>
#include <vector>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <thread>

using namespace std;

// Constants
const long long MAX_POPULATION = 1000000000000LL;

// Custom Line Chart Class
class LineChart : public Fl_Widget {
    vector<pair<double, double>> data;
    string xLabel, yLabel, title, key;
    double xMax = 1.0, yMax = 1.0;

public:
    LineChart(int x, int y, int w, int h, const char *l = 0)
        : Fl_Widget(x, y, w, h, l), xLabel("X"), yLabel("Y"), title("Line Chart"), key("") {}

    void setLabels(const string &x, const string &y, const string &t, const string &k) {
        xLabel = x;
        yLabel = y;
        title = t;
        key = k;
    }

    void addDataPoint(double x, double y) {
        data.emplace_back(x, y);
        xMax = max(xMax, x);
        yMax = max(yMax, y);
        if (data.size() > 20) data.erase(data.begin()); // Limit data points
        redraw();
    }

    void draw() override {
        fl_color(FL_WHITE);
        fl_rectf(x(), y(), w(), h()); // Clear background

        // Draw title
        fl_color(FL_BLACK);
        fl_font(FL_HELVETICA_BOLD, 14);
        fl_draw(title.c_str(), x() + w() / 2 - 50, y() + 20);

        // Draw axes
        fl_color(FL_BLACK);
        fl_line(x() + 70, y() + h() - 70, x() + w() - 40, y() + h() - 70); // X-axis
        fl_line(x() + 70, y() + h() - 70, x() + 70, y() + 40);             // Y-axis

        fl_font(FL_HELVETICA, 12);
        fl_draw(xLabel.c_str(), x() + w() / 2 - 20, y() + h() - 30);  // Align X axis label
        fl_draw(yLabel.c_str(), x() + 20, y() + h() / 2);  // Align Y axis label

        drawTicks(); // Draw ticks and numbers

        // Draw legend/key
        fl_color(FL_BLACK);
        fl_rect(x() + w() - 140, y() + 10, 120, 50);
        fl_font(FL_HELVETICA, 10);  // Legend font
        fl_draw(key.c_str(), x() + w() - 130, y() + 35);

        if (data.empty()) return;

        // Draw data points
        for (size_t i = 0; i < data.size() - 1; ++i) {
            double x1 = mapValue(data[i].first, 0, xMax, x() + 70, x() + w() - 40);
            double y1 = mapValue(data[i].second, 0, yMax, y() + h() - 70, y() + 40);
            double x2 = mapValue(data[i + 1].first, 0, xMax, x() + 70, x() + w() - 40);
            double y2 = mapValue(data[i + 1].second, 0, yMax, y() + h() - 70, y() + 40);

            fl_color(FL_BLUE);
            fl_line(x1, y1, x2, y2);
        }
    }

private:
    double mapValue(double value, double minVal, double maxVal, double outMin, double outMax) {
        return outMin + ((value - minVal) / (maxVal - minVal)) * (outMax - outMin);
    }

    void drawTicks() {
        fl_font(FL_HELVETICA, 10);
        fl_color(FL_BLACK);

        for (int i = 0; i <= 5; ++i) {
            double xTickVal = xMax * i / 5;
            double xTickPos = mapValue(xTickVal, 0, xMax, x() + 70, x() + w() - 40);
            fl_draw(to_string(xTickVal).c_str(), xTickPos - 10, y() + h() - 55);
            fl_line(xTickPos, y() + h() - 70, xTickPos, y() + h() - 65);

            double yTickVal = yMax * i / 5;
            double yTickPos = mapValue(yTickVal, 0, yMax, y() + h() - 70, y() + 40);
            fl_draw(to_string(yTickVal).c_str(), x() + 45, yTickPos + 5);
            fl_line(x() + 70, yTickPos, x() + 75, yTickPos);
        }
    }
};

// Global Variables
atomic<bool> running(false);
Fl_Input *inputPopulationSize, *inputMutationRate, *inputCrossoverRate, *inputThreads;
Fl_Choice *algorithmChoice;
Fl_Text_Display *outputDisplay;
Fl_Text_Buffer *outputBuffer;
LineChart *speedupChart, *efficiencyChart;

// Function Prototypes
void runSelectedAlgorithm(int algorithm, int threads, long long populationSize, double mutationRate, double crossoverRate);
void runStandardGA(int threads, long long populationSize, double mutationRate, double crossoverRate);
void runStationaryGA(int threads, long long populationSize, double mutationRate, double crossoverRate);
void runEliteGA(int threads, long long populationSize, double mutationRate, double crossoverRate);
void runIslandModelGA(int threads, long long populationSize, double mutationRate, double crossoverRate);

// Simulation Function
void simulateAlgorithm() {
    running = true;

    long long populationSize = atoll(inputPopulationSize->value());
    double mutationRate = atof(inputMutationRate->value());
    double crossoverRate = atof(inputCrossoverRate->value());
    int threads = atoi(inputThreads->value());

    // Validate Inputs
    if (populationSize <= 0 || populationSize > MAX_POPULATION || mutationRate < 0 || mutationRate > 1 || crossoverRate < 0 || crossoverRate > 1 || threads <= 0) {
        outputBuffer->append("Invalid input. Check parameters.\n");
        running = false;
        return;
    }

    int selectedAlgorithm = algorithmChoice->value();
    double baselineTime = 0.0;

    for (int currentThreads = 1; currentThreads <= threads && running; currentThreads *= 2) {
        omp_set_num_threads(currentThreads);

        auto startTime = chrono::high_resolution_clock::now();
        runSelectedAlgorithm(selectedAlgorithm, currentThreads, populationSize, mutationRate, crossoverRate);
        auto endTime = chrono::high_resolution_clock::now();

        double elapsed = chrono::duration<double>(endTime - startTime).count();
        if (currentThreads == 1) baselineTime = elapsed;
        double speedup = baselineTime / elapsed;
        double efficiency = speedup / currentThreads;

        speedupChart->addDataPoint(currentThreads, speedup);
        efficiencyChart->addDataPoint(currentThreads, efficiency);

        outputBuffer->append(("Threads: " + to_string(currentThreads) +
                              " Time: " + to_string(elapsed) + "s" +
                              " Speedup: " + to_string(speedup) +
                              " Efficiency: " + to_string(efficiency) + "\n").c_str());
    }

    running = false;
}

// Implement Genetic Algorithms
void runSelectedAlgorithm(int algorithm, int threads, long long populationSize, double mutationRate, double crossoverRate) {
    switch (algorithm) {
        case 0: runStandardGA(threads, populationSize, mutationRate, crossoverRate); break;
        case 1: runStationaryGA(threads, populationSize, mutationRate, crossoverRate); break;
        case 2: runEliteGA(threads, populationSize, mutationRate, crossoverRate); break;
        case 3: runIslandModelGA(threads, populationSize, mutationRate, crossoverRate); break;
        default: outputBuffer->append("Unknown algorithm selected.\n");
    }
}

void runStandardGA(int threads, long long populationSize, double mutationRate, double crossoverRate) {
    vector<double> fitness(populationSize);

#pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < populationSize; ++i) {
        fitness[i] = rand() / double(RAND_MAX); // Placeholder for actual Standard GA logic
    }
}

void runStationaryGA(int threads, long long populationSize, double mutationRate, double crossoverRate) {
    vector<double> fitness(populationSize);

#pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < populationSize; ++i) {
        fitness[i] = rand() / double(RAND_MAX);  // Placeholder: Implement actual Stationary GA logic
    }
}

void runEliteGA(int threads, long long populationSize, double mutationRate, double crossoverRate) {
    vector<double> fitness(populationSize);

#pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < populationSize; ++i) {
        fitness[i] = rand() / double(RAND_MAX);  // Placeholder: Implement actual Elite GA logic
    }
}

void runIslandModelGA(int threads, long long populationSize, double mutationRate, double crossoverRate) {
    vector<double> fitness(populationSize);

#pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < populationSize; ++i) {
        fitness[i] = rand() / double(RAND_MAX);  // Placeholder: Implement actual Island Model GA logic
    }
}

// Main Function
int main(int argc, char **argv) {
    srand(time(0));

    Fl_Window *window = new Fl_Window(950, 700, "Genetic Algorithm Simulator");

    inputPopulationSize = new Fl_Input(150, 30, 150, 25, "Population Size:");
    inputPopulationSize->value("100");

    inputMutationRate = new Fl_Input(150, 70, 150, 25, "Mutation Rate:");
    inputMutationRate->value("0.01");

    inputCrossoverRate = new Fl_Input(150, 110, 150, 25, "Crossover Rate:");
    inputCrossoverRate->value("0.7");

    inputThreads = new Fl_Input(150, 150, 150, 25, "Number of Threads:");
    inputThreads->value("4");

    algorithmChoice = new Fl_Choice(150, 190, 150, 25, "Algorithm:");
    algorithmChoice->add("Standard GA|Stationary GA|Elite GA|Island Model GA");
    algorithmChoice->value(0);

    // Buttons
    Fl_Button *startButton = new Fl_Button(50, 230, 100, 30, "Start");
    startButton->callback([](Fl_Widget *, void *) { thread(simulateAlgorithm).detach(); });

    Fl_Button *stopButton = new Fl_Button(200, 230, 100, 30, "Stop");
    stopButton->callback([](Fl_Widget *, void *) { running = false; });

    // Output Display
    outputDisplay = new Fl_Text_Display(50, 280, 400, 300);
    outputBuffer = new Fl_Text_Buffer();
    outputDisplay->buffer(outputBuffer);

    // Line Charts for Speedup and Efficiency
    speedupChart = new LineChart(500, 50, 400, 250);
    speedupChart->setLabels("Threads", "Speedup (x)", "Speedup vs. Threads", "Speedup");

    efficiencyChart = new LineChart(500, 350, 400, 250);
    efficiencyChart->setLabels("Threads", "Efficiency (%)", "Efficiency vs. Threads", "Efficiency");

    window->end();
    window->show(argc, argv);

    return Fl::run();
}

