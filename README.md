# Genetic-Algorithms-Simulation
Background:
Genetic Algorithms (GAs) are optimization techniques inspired by natural selection, useful
for complex problems. They require significant computational power, especially with large
datasets. This project addresses this by leveraging parallel computing to enhance GA
efficiency, providing a desktop application that simulates various GA models and showcases
parallelization benefits.
Project Overview:
This project includes the development of a desktop that presents a general-purpose genetic
algorithm (GA) in a friendly user interface designed for non-specialists. This application will
be aimed at demonstrating the advantages of GAs in solving optimization tasks and the
advantages of parallel computing in order to increase efficiency. The project is designed to be
developed in C++ using the FLTK GUI toolkit as the primary language.
Code Repository and Selection
Repository URL:
● MPI Distributed Genetic Algorithm

Code Description: The selected repositories for this project implement genetic algorithms
which are optimized using parallel
computing techniques (OpenMP and MPI).

Code Analysis: Upon examination of the code's complexity analysis reveals that it is
considered owing to the duration it takes to run on extensive datasets. Usually about O(n⋅g),
for a single thread operation. It encompasses parallelization commands like, over #pragmas
or MPI calls.


Parallelization Strategy
OpenMP and MPI Implementation: Identify sections of the code suitable for parallel
execution, such as fitness evaluation and selection processes.

Execution Plan
● Hardware Specifications: Use a multi-core processor with at least 8 cores and 16GB
RAM for testing.
● Baseline Execution: Record single-thread execution time.
● Parallel Execution: Test performance with varying numbers of threads/processes and
different input sizes.
Data and Performance Metrics
Data Size: Utilize large datasets to test the efficiency of the algorithms.
Performance Metrics:
● Metrics include execution time, speedup, scalability, and efficiency.
● Measure performance by comparing original and parallelized versions.
Testing and Validation
Expected Results: Anticipate significant speedup in execution time through parallelization.


Deliverables:
● A fully functional desktop application that demonstrates genetic algorithms with a GUI.
● Application source code and documentation.
● Performance analysis comparing standard and parallel implementations.
