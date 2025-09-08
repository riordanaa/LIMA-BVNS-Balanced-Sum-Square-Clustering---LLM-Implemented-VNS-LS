//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================
#include "Vns.h"
#include "Solution.h"
#include "LocalSearch.h"
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <ctime>
#include <cstdlib>
#include <random>
#include "tempsC++.h"
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "Pair.h"
#include <limits>
#include <algorithm>
using namespace std; 

struct c_unique {
  int current;
  c_unique() {current=0;}
  int operator()() {return current++;}
} UniqueNumber;

Vns::Vns(vector<Point>* _dataset, DistanceMatrix* _distances, int _nClusters, Random* _random, vector< vector<Pair> >* _rankedEntities) {
    dataset = _dataset;
    distances = _distances;
    nClusters = _nClusters;
    random = _random;
    rankedEntities = _rankedEntities;
    
    // Initialize k to 1 (smallest neighborhood)
    k = 1;
}

int Vns::execute(Solution& bestSolution, double timeMax, int kMin, int kStep, int kMax, string outputFileName) {
    int iteration = 0;
    
    // Start the timer
    timer.Start();
    double startTime = timer.GetTime(); // Record start time
    
    // Generate initial solution
    initialSolution(bestSolution);
    
    // Create local search object
    LocalSearch localSearch(dataset, random, rankedEntities);
    
    // Apply initial local search to improve starting solution
    localSearch.execute(bestSolution, &timer, timeMax, 100);
    
    // Current solution for exploration
    Solution currentSolution(bestSolution);
    
    // Best solution value found so far
    double bestValue = bestSolution.solutionValue;
    
    // Output file for tracking progress
    ofstream outputFile;
    if (!outputFileName.empty()) {
        outputFile.open(outputFileName.c_str(), ios::out);
        outputFile << "Iteration,Time,k,SolutionValue" << endl;
        outputFile << iteration << "," << timer.GetTime() << "," << k << "," << bestValue << endl;
    }
    
    // Log initial solution
    cout << "Initial solution: " << fixed << setprecision(6) << bestValue << endl;
    
    // Main VNS loop
    k = kMin; // Start with minimum neighborhood size
    while (timer.GetTime() < timeMax) {
        iteration++;
        
        // Copy best solution to current solution for this iteration
        currentSolution.copy(bestSolution);
        
        // SHAKING PHASE: perturb current solution with k random swaps
        if (!shaking(currentSolution)) {
            // If shaking fails (e.g., insufficient clusters), skip iteration
            continue;
        }
        
        // LOCAL SEARCH PHASE: improve the perturbed solution
        localSearch.execute(currentSolution, &timer, timeMax, 100);
        
        // Check if we've exceeded the time limit
        if (timer.GetTime() >= timeMax) {
            break;
        }
        
        // Update best solution if improvement found
        if (currentSolution.solutionValue < bestSolution.solutionValue) {
            // Copy improved solution
            bestSolution.copy(currentSolution);
            bestValue = bestSolution.solutionValue;
            
            // Reset k to minimum (key VNS principle)
            k = kMin;
            
            // Log progress
            cout << "Iter: " << iteration 
                 << ", Time: " << fixed << setprecision(2) << timer.GetTime() 
                 << "s, k: " << k 
                 << ", Value: " << fixed << setprecision(6) << bestValue << endl;
            
            if (outputFile.is_open()) {
                outputFile << iteration << "," << timer.GetTime() << "," << k << "," << bestValue << endl;
            }
        } else {
            // No improvement - move to larger neighborhood
            k += kStep;
            
            // If k exceeds maximum, reset to minimum (cycle through neighborhoods)
            if (k > kMax) {
                k = kMin;
            }
        }
        
        // Periodically validate solution integrity
        if (iteration % 100 == 0 && !checkSolution(&bestSolution)) {
            cerr << "Error: Invalid solution detected at iteration " << iteration << endl;
            break;
        }
    }
    
    // Close output file if open
    if (outputFile.is_open()) {
        outputFile.close();
    }
    
    // Final solution check
    if (!checkSolution(&bestSolution)) {
        cerr << "Warning: Final solution fails verification" << endl;
    }
    
    // Calculate total execution time
    double totalTime = timer.GetTime() - startTime;
    
    cout << "VNS completed " << iteration << " iterations in " 
         << timer.GetTime() << " seconds." << endl;
    cout << "Final solution value: " << fixed << setprecision(6) << bestSolution.solutionValue << endl;
    
    // Add explicit timing information that will be visible in the output
    cout << endl;
    cout << "************************** LIMA-VNS TIMING SUMMARY **************************" << endl;
    cout << "* Dataset size: " << bestSolution.nDataPoints << " points" << endl; 
    cout << "* Number of clusters: " << bestSolution.nClusters << endl;
    cout << "* Total execution time: " << fixed << setprecision(6) << totalTime << " seconds" << endl;
    cout << "* Best objective function value: " << fixed << setprecision(6) << bestSolution.solutionValue << endl;
    cout << "************************************************************************" << endl;
    cout << endl;
    
    // Store the execution time in the solution object to make it available to the caller
    bestSolution.time = totalTime;
    
    return iteration;
}

bool Vns::shaking(Solution& solution) {
    // Perform k random swaps - this is the core of VNS "shaking" phase
    for (int swapCount = 0; swapCount < k; swapCount++) {
        // Select two random different clusters
        int clusterI = static_cast<int>(random->randp() * nClusters);
        int clusterJ;
        do {
            clusterJ = static_cast<int>(random->randp() * nClusters);
        } while (clusterJ == clusterI);
        
        // Get points assigned to each cluster
        vector<int> pointsInClusterI;
        vector<int> pointsInClusterJ;
        
        for (int p = 0; p < solution.nDataPoints; p++) {
            if (solution.assignment[p] == clusterI) {
                pointsInClusterI.push_back(p);
            } else if (solution.assignment[p] == clusterJ) {
                pointsInClusterJ.push_back(p);
            }
        }
        
        // Safety check - clusters should have points for balanced clustering
        if (pointsInClusterI.empty() || pointsInClusterJ.empty()) {
            swapCount--; // Try again for this swap
            continue;
        }
        
        // Select random points from each cluster
        int pointV = pointsInClusterI[static_cast<int>(random->randp() * pointsInClusterI.size())];
        int pointU = pointsInClusterJ[static_cast<int>(random->randp() * pointsInClusterJ.size())];
        
        // Update assignments only during shaking
        solution.assignment[pointV] = clusterJ;
        solution.assignment[pointU] = clusterI;
        
        // Update SC matrix for all points
        for (int p = 0; p < solution.nDataPoints; p++) {
            // Update SC values for cluster I
            solution.sc[p][clusterI] = solution.sc[p][clusterI] 
                                     - solution.distances->getDistance(p, pointV) 
                                     + solution.distances->getDistance(p, pointU);
            
            // Update SC values for cluster J
            solution.sc[p][clusterJ] = solution.sc[p][clusterJ] 
                                     - solution.distances->getDistance(p, pointU) 
                                     + solution.distances->getDistance(p, pointV);
        }
    }
    
    // KEY FIX: After all swaps in shaking, recalculate the solution value from scratch
    // This prevents error accumulation while maintaining the efficiency of incremental SC updates
    solution.solutionValue = 0;
    for (int c = 0; c < solution.nClusters; c++) {
        double sumDistances = 0;
        int clusterSize = 0; // Count for verification
        
        for (int p = 0; p < solution.nDataPoints; p++) {
            if (solution.assignment[p] == c) {
                sumDistances += solution.sc[p][c];
                clusterSize++;
            }
        }
        
        // Add to solution value (divide by 2 because each distance is counted twice)
        solution.solutionValue += sumDistances / (2 * clusterSize);
    }
    
    return true;
}

void Vns::initialSolution(Solution& initial) {
    // Initialize a balanced clustering solution
    
    // Determine cluster sizes to ensure balance
    int pointsPerCluster = initial.nDataPoints / initial.nClusters;
    int remainder = initial.nDataPoints % initial.nClusters;
    
    // First, randomize point order
    vector<int> pointIndices(initial.nDataPoints);
    for (int i = 0; i < initial.nDataPoints; i++) {
        pointIndices[i] = i;
    }
    
    // Shuffle the indices (Fisher-Yates algorithm)
    for (int i = 0; i < initial.nDataPoints - 1; i++) {
        int range = initial.nDataPoints - i;
        int offset = static_cast<int>(random->randp() * range);
        int j = i + offset;
        std::swap(pointIndices[i], pointIndices[j]);
    }
    
    // Assign points to clusters ensuring balance
    int currentPoint = 0;
    for (int c = 0; c < initial.nClusters; c++) {
        // Calculate size for this cluster (handling remainder)
        int clusterSize = pointsPerCluster + (c < remainder ? 1 : 0);
        initial.clusterSizes[c] = clusterSize;
        
        // Assign points to this cluster
        for (int i = 0; i < clusterSize; i++) {
            initial.assignment[pointIndices[currentPoint]] = c;
            currentPoint++;
        }
    }
    
    // Initialize the SC matrix values
    initial.initializeSc();
    
    // Calculate initial solution value using Huygens' theorem
    initial.solutionValue = 0;
    
    for (int c = 0; c < initial.nClusters; c++) {
        double sumDistances = 0;
        
        for (int p = 0; p < initial.nDataPoints; p++) {
            if (initial.assignment[p] == c) {
                sumDistances += initial.sc[p][c];
            }
        }
        
        // Add to solution value (divide by 2 because each distance is counted twice)
        initial.solutionValue += sumDistances / (2 * initial.clusterSizes[c]);
    }
}

// Add this function to Vns.cpp
void Vns::loadInitialSolution(Solution& solution, const std::string& filename) {
    std::ifstream infile(filename, std::ios::binary);
    
    if (!infile) {
        std::cerr << "Error opening initial solution file: " << filename << std::endl;
        exit(1);
    }
    
    // Read version marker
    int version;
    infile.read(reinterpret_cast<char*>(&version), sizeof(int));
    
    if (version != 1) {
        std::cerr << "Unknown version format in initial solution file" << std::endl;
        exit(1);
    }
    
    // Read initialization time
    double initTime;
    infile.read(reinterpret_cast<char*>(&initTime), sizeof(double));
    
    // Read dimensions to verify
    int nDataPoints, nClusters;
    infile.read(reinterpret_cast<char*>(&nDataPoints), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));
    
    if (nDataPoints != solution.nDataPoints || nClusters != solution.nClusters) {
        std::cerr << "Error: Initial solution dimensions don't match" << std::endl;
        exit(1);
    }
    
    // Read cluster assignments
    for (int i = 0; i < solution.nDataPoints; i++) {
        infile.read(reinterpret_cast<char*>(&solution.assignment[i]), sizeof(int));
    }
    
    // Read cluster sizes
    for (int i = 0; i < solution.nClusters; i++) {
        infile.read(reinterpret_cast<char*>(&solution.clusterSizes[i]), sizeof(double));
    }
    
    // Initialize SC matrix
    solution.initializeSc();
    
    // Calculate solution value
    solution.solutionValue = 0.0;
    for (int i = 0; i < solution.nDataPoints - 1; i++) {
        for (int j = i + 1; j < solution.nDataPoints; j++) {
            if (solution.assignment[i] == solution.assignment[j]) {
                solution.solutionValue += solution.distances->getDistance(i, j) / 
                                          solution.clusterSizes[solution.assignment[i]];
            }
        }
    }
    
    // Account for the initialization time in the timer
    timer.Reset();
    timer.Start();
    
    // Simulate elapsed time to include initialization time
    // This is a bit of a hack but it works with the ChronoCPU class
    std::this_thread::sleep_for(std::chrono::milliseconds(static_cast<int>(initTime * 1000)));
    
    infile.close();
}

bool Vns::checkSolution(Solution* solution) {
    // Comprehensive solution verification
    
    // Verify all points have valid cluster assignments
    for (int p = 0; p < solution->nDataPoints; p++) {
        int c = solution->assignment[p];
        if (c < 0 || c >= solution->nClusters) {
            cerr << "Error: Point " << p << " assigned to invalid cluster " << c << endl;
            return false;
        }
    }
    
    // Count points in each cluster
    vector<int> clusterSizes(solution->nClusters, 0);
    for (int p = 0; p < solution->nDataPoints; p++) {
        clusterSizes[solution->assignment[p]]++;
    }
    
    // Verify balanced clustering constraint
    int expectedSize = solution->nDataPoints / solution->nClusters;
    int remainder = solution->nDataPoints % solution->nClusters;
    
    for (int c = 0; c < solution->nClusters; c++) {
        int expectedClusterSize = expectedSize + (c < remainder ? 1 : 0);
        if (clusterSizes[c] != expectedClusterSize) {
            cerr << "Balance constraint violated: Cluster " << c << " has " 
                 << clusterSizes[c] << " points, expected " << expectedClusterSize << endl;
            return false;
        }
    }
    
    // Verify SC matrix consistency
    for (int p = 0; p < solution->nDataPoints; p++) {
        for (int c = 0; c < solution->nClusters; c++) {
            double calculatedSC = 0.0;
            
            for (int q = 0; q < solution->nDataPoints; q++) {
                if (solution->assignment[q] == c) {
                    calculatedSC += solution->distances->getDistance(p, q);
                }
            }
            
            if (fabs(calculatedSC - solution->sc[p][c]) > 1e-6) {
                cerr << "SC matrix error at [" << p << "][" << c << "]: "
                     << "Calculated=" << calculatedSC << ", Stored=" << solution->sc[p][c] << endl;
                return false;
            }
        }
    }
    
    // Recalculate solution value from scratch
    double recalculatedValue = 0;
    for (int c = 0; c < solution->nClusters; c++) {
        double sumDistances = 0;
        
        for (int p = 0; p < solution->nDataPoints; p++) {
            if (solution->assignment[p] == c) {
                sumDistances += solution->sc[p][c];
            }
        }
        
        // Using Huygens' theorem for accurate calculation
        recalculatedValue += sumDistances / (2 * clusterSizes[c]);
    }
    
    // Verify solution value
    if (fabs(recalculatedValue - solution->solutionValue) > 1e-6) {
        cerr << "Solution value error: Calculated=" << recalculatedValue 
             << ", Stored=" << solution->solutionValue << endl;
        return false;
    }
    
    // All checks passed
    return true;
}