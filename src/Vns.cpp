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

// This private helper function is added to avoid creating new LocalSearch objects in the shaking phase.
void Vns::swap(Solution& solution, int pointI, int pointJ, double delta) {
    int clusterI = solution.assignment[pointI];
    int clusterJ = solution.assignment[pointJ];
    solution.solutionValue += delta;
    for (int k = 0; k < solution.nDataPoints; k++) {
        double dist_k_I = solution.distances->getDistance(k, pointI);
        double dist_k_J = solution.distances->getDistance(k, pointJ);
        solution.sc[k][clusterI] = solution.sc[k][clusterI] - dist_k_I + dist_k_J;
        solution.sc[k][clusterJ] = solution.sc[k][clusterJ] + dist_k_I - dist_k_J;
    }
    solution.assignment[pointI] = clusterJ;
    solution.assignment[pointJ] = clusterI;
}

// Constructor to initialize the VNS algorithm parameters
Vns::Vns(vector<Point>* _dataset, DistanceMatrix *_distances, int _nClusters, Random* _random, vector< vector<Pair> >* _rankedEntities) {
    dataset = _dataset;
    distances = _distances;
    nClusters = _nClusters;
    random = _random;
    rankedEntities = _rankedEntities;
    k = 1; // Initialize neighborhood size
    timer = ChronoCPU();
}

// Main execution method for the VNS algorithm
int Vns::execute(Solution& bestSolution, double tempMax, int kMin, int kStep, int kMax, string outputFileName) {
    timer.Start();
    int iter = 0;

    LocalSearch localSearch(dataset, random, rankedEntities);

    // 1. Generate a random, balanced initial solution
    initialSolution(bestSolution);
    // 2. Improve it with local search to find the first local optimum
    localSearch.execute(bestSolution, &timer, tempMax, iter);

    cout << "Initial Solution Value: " << fixed << setprecision(5) << bestSolution.solutionValue << endl;
    
    k = kMin; // Start with the smallest neighborhood size

    // Main VNS loop - CORRECTED to only check time limit
    while (timer.GetTime() < tempMax) {
        iter++;

        // 1. Create a working copy of the current best solution
        Solution currentSolution(bestSolution);

        // 2. Shaking: Perturb the solution by applying 'k' random swaps
        shaking(currentSolution);
        
        // 3. Local Search: Find the local optimum from the shaken solution
        localSearch.execute(currentSolution, &timer, tempMax, iter);

        // 4. Move or Stay: Compare the new local optimum with the best-so-far solution
        if (currentSolution.solutionValue < bestSolution.solutionValue - 1e-9) {
            bestSolution.copy(currentSolution);
            cout << "Iteration " << iter << ": Found new best solution = " << fixed << setprecision(5) << bestSolution.solutionValue << " (k=" << k << ")" << endl;
            k = kMin; // Improvement found, reset to the smallest neighborhood
        } else {
            k += kStep; // No improvement, increase the neighborhood size
            // CORRECTED: If k exceeds kMax, reset it to kMin to continue searching
            if (k > kMax) {
                k = kMin;
            }
        }
    }
    
    timer.Stop();
    cout << "VNS finished. Total iterations: " << iter << endl;
    cout << "Final best solution value: " << fixed << setprecision(5) << bestSolution.solutionValue << endl;
    cout << "Total time: " << timer.GetTime() << "s" << endl;
    
    return iter;
}


// Shaking function: Applies 'k' random swaps to the solution
bool Vns::shaking(Solution& solution) {
    // OPTIMIZED: Calls the internal swap method directly
    for (int i = 0; i < k; ++i) {
        int pointA, pointB;
        
        // Find two points in different clusters to ensure a valid swap
        do {
            pointA = random->get_rand(solution.nDataPoints);
            pointB = random->get_rand(solution.nDataPoints);
        } while (pointA == pointB || solution.assignment[pointA] == solution.assignment[pointB]);

        int clusterA = solution.assignment[pointA];
        int clusterB = solution.assignment[pointB];
        double dist_ab = solution.distances->getDistance(pointA, pointB);

        // Calculate the delta for this random swap
        double delta = ( (solution.sc[pointB][clusterA] - solution.sc[pointA][clusterA] - dist_ab) / solution.clusterSizes[clusterA] ) +
                       ( (solution.sc[pointA][clusterB] - solution.sc[pointB][clusterB] - dist_ab) / solution.clusterSizes[clusterB] );
        
        // Apply the swap using the efficient internal method
        this->swap(solution, pointA, pointB, delta);
    }
    return true;
}

// Generates a random, balanced initial solution
void Vns::initialSolution(Solution& initial) {
    int n = initial.nDataPoints;
    int k = initial.nClusters;

    vector<int> p(n);
    for(int i = 0; i < n; ++i) p[i] = i;
    random->random_shuffle(p.begin(), p.end());

    // Determine balanced cluster sizes
    int baseSize = n / k;
    int remainder = n % k;

    for(int i = 0; i < k; i++){
        initial.clusterSizes[i] = baseSize + (i < remainder ? 1 : 0);
    }

    // Assign points to clusters based on shuffled order
    int pointIdx = 0;
    for (int c = 0; c < k; c++) {
        for (int i = 0; i < initial.clusterSizes[c]; i++) {
            initial.assignment[p[pointIdx++]] = c;
        }
    }

    // Initialize the sc matrix based on the new assignments
    initial.initializeSc();

    // Calculate the initial solution value from scratch
    initial.solutionValue = 0;
    for (int c = 0; c < k; c++) {
        double intraClusterSum = 0;
        for (int i = 0; i < n; i++) {
            if (initial.assignment[i] == c) {
                // Using the sc matrix for efficient calculation
                // sc[i][c] is sum of distances from i to all points in cluster c
                intraClusterSum += initial.sc[i][c];
            }
        }
        // Each pair distance is counted twice in the sum, so divide by 2
        initial.solutionValue += (intraClusterSum / 2.0) / initial.clusterSizes[c];
    }
}


// Verification function to check solution integrity (for debugging)
bool Vns::checkSolution(Solution* solution) {
    double calculatedValue = 0;
    int n = solution->nDataPoints;
    int k = solution->nClusters;

    // 1. Recalculate solution value from scratch
    for (int c = 0; c < k; c++) {
        double intraClusterSum = 0;
        vector<int> pointsInCluster;
        for (int i = 0; i < n; i++) {
            if (solution->assignment[i] == c) {
                pointsInCluster.push_back(i);
            }
        }

        if (pointsInCluster.size() > 1) {
            for (size_t i = 0; i < pointsInCluster.size(); ++i) {
                for (size_t j = i + 1; j < pointsInCluster.size(); ++j) {
                    intraClusterSum += solution->distances->getDistance(pointsInCluster[i], pointsInCluster[j]);
                }
            }
            calculatedValue += intraClusterSum / pointsInCluster.size();
        }
    }

    // 2. Verify stored value matches recalculated value
    if (fabs(calculatedValue - solution->solutionValue) > 1e-6) {
        cerr << "ERROR: Stored solution value (" << solution->solutionValue 
             << ") does not match recalculated value (" << calculatedValue << ")." << endl;
        return false;
    }

    // 3. Confirm all clusters have the correct size (balance constraint)
    vector<int> currentSizes(k, 0);
    for(int i=0; i<n; ++i) {
        currentSizes[solution->assignment[i]]++;
    }

    int baseSize = n / k;
    int remainder = n % k;
    for(int c=0; c<k; ++c) {
        double expectedSize = baseSize + (c < remainder ? 1 : 0);
        if(currentSizes[c] != expectedSize) {
            cerr << "ERROR: Cluster " << c << " has incorrect size. Expected: " << expectedSize << ", Got: " << currentSizes[c] << endl;
            return false;
        }
    }
    
    return true; // Solution is valid
}

void Vns::loadInitialSolution(Solution& solution, const std::string& filename) {
    // This function is a stub and would need to be implemented if loading from a file is required.
    // It would involve reading a file that specifies the assignment of each point to a cluster.
    // After loading assignments, it must call initialSolution logic to set clusterSizes, sc matrix, and solutionValue.
    cerr << "Warning: loadInitialSolution is not implemented." << endl;
}