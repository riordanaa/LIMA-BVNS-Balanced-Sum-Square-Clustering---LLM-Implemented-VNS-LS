//============================================================================
// Author      : Leandro R. Costa, Daniel Aloise, Nenad Mladenovic.
// Description : Implementation of the LIMA-VNS published in the paper "Less is 
//               more: basic variable neighborhood search heuristic for 
//               balanced minimum sum-of-squares clustering". Please, check
//               https://doi.org/10.1016/j.ins.2017.06.019 for theoretical 
//               details. 
//============================================================================

#include "LocalSearch.h"
#include "Solution.h"
#include <vector>
#include "tempsC++.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include "Pair.h"
#include <map>
#include <random>
#include <cmath> // For fabs

using namespace std;

// Constructor to initialize the LocalSearch object
LocalSearch::LocalSearch(vector<Point>* _dataset, Random* _random, vector< vector<Pair> >* _rankedEntities) {
    dataset = _dataset;
    random = _random;
    rankedEntities = _rankedEntities;
}

// Main execution wrapper for the local search process.
// It repeatedly applies the first-improvement heuristic until a local minimum is reached.
void LocalSearch::execute(Solution& bestLocalSolution, ChronoCPU* timer, double maxTime, int nIteration) {
    // Continuously apply the first-improvement swap search until no more improvements can be found.
    while (swapLocalSearchFirstRand(bestLocalSolution, timer, maxTime));
}

// Performs a best-improvement search.
// It evaluates all possible swaps and executes the one that provides the maximum improvement.
// Note: The LIMA-VNS paper uses a first-improvement strategy, but this is included for completeness.
bool LocalSearch::swapLocalSearchBest(Solution& solution, ChronoCPU* timer, double maxTime) {
    double bestDelta = 1e-9; // Use a small positive epsilon to avoid floating point noise
    int bestI = -1, bestJ = -1;

    for (int i = 0; i < solution.nDataPoints; i++) {
        for (int j = i + 1; j < solution.nDataPoints; j++) {
            if (timer->GetTime() > maxTime) return false;

            int clusterI = solution.assignment[i];
            int clusterJ = solution.assignment[j];

            if (clusterI == clusterJ) continue;

            double dist_ij = solution.distances->getDistance(i, j);
            
            // Calculate the change in objective function (delta) for swapping points i and j.
            // This is the O(1) calculation derived from Huygens' theorem.
            double delta = ( (solution.sc[j][clusterI] - solution.sc[i][clusterI] - dist_ij) / solution.clusterSizes[clusterI] ) +
                           ( (solution.sc[i][clusterJ] - solution.sc[j][clusterJ] - dist_ij) / solution.clusterSizes[clusterJ] );


            if (delta < bestDelta) {
                bestDelta = delta;
                bestI = i;
                bestJ = j;
            }
        }
    }

    if (bestI != -1) {
        // An improving swap was found, execute it.
        swap(solution, bestI, bestJ, bestDelta);
        return true;
    }

    return false; // No improvement found
}

// Performs a first-improvement search with a randomized starting point.
// It iterates through all possible swaps and executes the *first* one that improves the solution.
bool LocalSearch::swapLocalSearchFirstRand(Solution& solution, ChronoCPU* timer, double maxTime) {
    // Create a shuffled list of indices to randomize the search starting point
    vector<int> indices(solution.nDataPoints);
    for(int i = 0; i < solution.nDataPoints; ++i) indices[i] = i;
    random->random_shuffle(indices.begin(), indices.end());

    for (int i_idx = 0; i_idx < solution.nDataPoints; ++i_idx) {
        int i = indices[i_idx];
        for (int j_idx = i_idx + 1; j_idx < solution.nDataPoints; ++j_idx) {
            int j = indices[j_idx];

            if (timer->GetTime() > maxTime) return false;

            int clusterI = solution.assignment[i];
            int clusterJ = solution.assignment[j];

            if (clusterI == clusterJ) continue;

            double dist_ij = solution.distances->getDistance(i, j);

            // Calculate the change in objective function (delta) for swapping points i and j.
            // This is the O(1) calculation derived from Huygens' theorem.
            double delta = ( (solution.sc[j][clusterI] - solution.sc[i][clusterI] - dist_ij) / solution.clusterSizes[clusterI] ) +
                           ( (solution.sc[i][clusterJ] - solution.sc[j][clusterJ] - dist_ij) / solution.clusterSizes[clusterJ] );
            
            // If the delta is negative (an improvement), perform the swap and exit immediately.
            if (delta < -1e-9) {
                swap(solution, i, j, delta);
                return true; // Improvement found and applied
            }
        }
    }

    return false; // No improvement found after checking all pairs
}

// The core swap operation.
// Updates the solution value in O(1) and the sc matrix in O(n).
void LocalSearch::swap(Solution& solution, int pointI, int pointJ, double delta) {
    int clusterI = solution.assignment[pointI];
    int clusterJ = solution.assignment[pointJ];

    // 1. Update the solution value in O(1) using the pre-calculated delta
    solution.solutionValue += delta;

    // 2. Update the sc matrix in O(n)
    // For every point k, adjust its summed-distance for clusterI and clusterJ
    // to reflect the swap of pointI and pointJ.
    for (int k = 0; k < solution.nDataPoints; k++) {
        double dist_k_I = solution.distances->getDistance(k, pointI);
        double dist_k_J = solution.distances->getDistance(k, pointJ);

        // For clusterI, remove pointI's contribution and add pointJ's
        solution.sc[k][clusterI] = solution.sc[k][clusterI] - dist_k_I + dist_k_J;
        
        // For clusterJ, add pointI's contribution and remove pointJ's
        solution.sc[k][clusterJ] = solution.sc[k][clusterJ] + dist_k_I - dist_k_J;
    }

    // 3. Update the point assignments *after* the sc matrix has been updated
    solution.assignment[pointI] = clusterJ;
    solution.assignment[pointJ] = clusterI;
}


// A validation/debugging function to check solution integrity.
// It re-calculates the objective function from scratch (slowly) and compares it
// to the incrementally updated solutionValue. Returns true if they match.
bool LocalSearch::checkSolution(Solution* solutionBefore, Solution* solutionAfter, double deltaSolutionValue) {
    double calculatedValue = 0;
    
    // Recalculate the cost of each cluster from scratch
    for (int j = 0; j < solutionAfter->nClusters; j++) {
        double clusterCost = 0;
        vector<int> pointsInCluster;
        for (int i = 0; i < solutionAfter->nDataPoints; i++) {
            if (solutionAfter->assignment[i] == j) {
                pointsInCluster.push_back(i);
            }
        }
        
        if (pointsInCluster.size() > 1) {
            double intraClusterSum = 0;
            for (size_t p1_idx = 0; p1_idx < pointsInCluster.size(); ++p1_idx) {
                for (size_t p2_idx = p1_idx + 1; p2_idx < pointsInCluster.size(); ++p2_idx) {
                    intraClusterSum += solutionAfter->distances->getDistance(pointsInCluster[p1_idx], pointsInCluster[p2_idx]);
                }
            }
            clusterCost = intraClusterSum / pointsInCluster.size();
        }
        calculatedValue += clusterCost;
    }

    // Compare the from-scratch calculation with the incrementally updated value
    if (fabs(calculatedValue - solutionAfter->solutionValue) > 1e-6) {
        cerr << "ERROR: Solution value mismatch after local search." << endl;
        cerr << "Expected: " << calculatedValue << ", Got: " << solutionAfter->solutionValue << endl;
        return false;
    }

    // Also check if the delta was calculated correctly
    if (fabs((solutionBefore->solutionValue + deltaSolutionValue) - solutionAfter->solutionValue) > 1e-6) {
		cerr << "ERROR: Delta calculation mismatch." << endl;
		cerr << "Before: " << solutionBefore->solutionValue << " + Delta: " << deltaSolutionValue << " != After: " << solutionAfter->solutionValue << endl;
		return false;
	}

    return true;
}