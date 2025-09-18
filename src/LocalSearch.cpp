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
using namespace std;

LocalSearch::LocalSearch(vector<Point>* _dataset, Random* _random, vector< vector<Pair> >* _rankedEntities) {
    dataset = _dataset;
    random = _random;
    rankedEntities = _rankedEntities;
}

void LocalSearch::execute(Solution& bestLocalSolution, ChronoCPU* timer, double maxTime, int nIteration) {
    bool improved = true;
    int iterations = 0;
    
    // Apply first-improvement local search until no more improvements or time limit reached
    while (improved && timer->GetTime() < maxTime && iterations < nIteration) {
        improved = swapLocalSearchFirstRand(bestLocalSolution, timer, maxTime);
        iterations++;
    }
}

bool LocalSearch::swapLocalSearchBest(Solution& solution, ChronoCPU* timer, double maxTime) {
    double bestDeltaF = 0; // No improvement initially
    int bestI = -1, bestPointV = -1, bestJ = -1, bestPointU = -1;
    bool improved = false;
    
    // Iterate through all possible point pairs
    for (int i = 0; i < solution.nClusters; i++) {
        // Check time limit
        if (timer->GetTime() >= maxTime) return false;
        
        for (int j = i + 1; j < solution.nClusters; j++) {
            // For each cluster pair, iterate through all point pairs
            for (int pointV = 0; pointV < solution.nDataPoints; pointV++) {
                // Check if point belongs to cluster i
                if (solution.assignment[pointV] != i) continue;
                
                for (int pointU = 0; pointU < solution.nDataPoints; pointU++) {
                    // Periodically check time limit
                    if ((pointU % 100) == 0 && timer->GetTime() >= maxTime) 
                        return improved;
                    
                    // Check if point belongs to cluster j
                    if (solution.assignment[pointU] != j) continue;
                    
                    // Calculate potential improvement using Huygens' theorem
                    double sizeI = solution.clusterSizes[i];
                    double sizeJ = solution.clusterSizes[j];
                    
                    // FIXED: Cost components for the swap - corrected formula
                    // When pointV leaves cluster i
                    double removeVCost = -solution.sc[pointV][i] / (2 * sizeI);
                    
                    // When pointU leaves cluster j
                    double removeUCost = -solution.sc[pointU][j] / (2 * sizeJ);
                    
                    // When pointV joins cluster j
                    double addVCost = solution.sc[pointV][j] / (2 * sizeJ);
                    
                    // When pointU joins cluster i
                    double addUCost = solution.sc[pointU][i] / (2 * sizeI);
                    
                    // Total change in solution value
                    double deltaF = removeVCost + removeUCost + addVCost + addUCost;
                    
                    // Update best if improvement found
                    if (deltaF < bestDeltaF) {
                        bestDeltaF = deltaF;
                        bestI = i;
                        bestJ = j;
                        bestPointV = pointV;
                        bestPointU = pointU;
                        improved = true;
                    }
                }
            }
        }
    }
    
    // If improvement found, perform swap
    if (improved) {
        swap(solution, bestI, bestPointV, bestJ, bestPointU, bestDeltaF);
    }
    
    return improved;
}

bool LocalSearch::swapLocalSearchFirstRand(Solution& solution, ChronoCPU* timer, double maxTime) {
    // Create randomized order of clusters to search
    vector<int> clusterIndices(solution.nClusters);
    for (int i = 0; i < solution.nClusters; i++) {
        clusterIndices[i] = i;
    }
    
    // Randomize order (using the Random class provided)
    for (int i = 0; i < solution.nClusters - 1; i++) {
        int range = solution.nClusters - i;
        // FIXED: Ensure offset is always < range to avoid buffer overflow
        int offset = static_cast<int>(random->randp() * range);
        if (offset >= range) offset = range - 1; // Safety check
        int j = i + offset;
        std::swap(clusterIndices[i], clusterIndices[j]);
    }
    
    // For each cluster pair in random order
    for (int ci = 0; ci < solution.nClusters - 1; ci++) {
        int i = clusterIndices[ci];
        
        for (int cj = ci + 1; cj < solution.nClusters; cj++) {
            int j = clusterIndices[cj];
            
            // Create vectors of points in each cluster
            vector<int> pointsInI;
            vector<int> pointsInJ;
            
            for (int p = 0; p < solution.nDataPoints; p++) {
                if (solution.assignment[p] == i) {
                    pointsInI.push_back(p);
                } else if (solution.assignment[p] == j) {
                    pointsInJ.push_back(p);
                }
            }
            
            // Skip if either cluster is empty (shouldn't happen in balanced clustering)
            if (pointsInI.empty() || pointsInJ.empty()) {
                continue;
            }
            
            // Randomize order of points with proper bounds checking
            if (pointsInI.size() > 1) {
                for (size_t pi = 0; pi < pointsInI.size() - 1; pi++) {
                    size_t range = pointsInI.size() - pi;
                    // FIXED: Ensure offset is always < range to avoid buffer overflow
                    size_t offset = static_cast<size_t>(random->randp() * range);
                    if (offset >= range) offset = range - 1; // Safety check
                    size_t pj = pi + offset;
                    std::swap(pointsInI[pi], pointsInI[pj]);
                }
            }
            
            if (pointsInJ.size() > 1) {
                for (size_t pj = 0; pj < pointsInJ.size() - 1; pj++) {
                    size_t range = pointsInJ.size() - pj;
                    // FIXED: Ensure offset is always < range to avoid buffer overflow
                    size_t offset = static_cast<size_t>(random->randp() * range);
                    if (offset >= range) offset = range - 1; // Safety check
                    size_t pk = pj + offset;
                    std::swap(pointsInJ[pj], pointsInJ[pk]);
                }
            }
            
            // Check all point pairs in random order
            for (size_t pointVIdx = 0; pointVIdx < pointsInI.size(); pointVIdx++) {
                int pointV = pointsInI[pointVIdx];
                
                for (size_t pointUIdx = 0; pointUIdx < pointsInJ.size(); pointUIdx++) {
                    // Check time limit periodically
                    if ((pointUIdx % 20) == 0 && timer->GetTime() >= maxTime) 
                        return false;
                    
                    int pointU = pointsInJ[pointUIdx];
                    
                    // Calculate potential improvement
                    double sizeI = solution.clusterSizes[i];
                    double sizeJ = solution.clusterSizes[j];
                    
                    // FIXED: Cost components for the swap - corrected formula
                    // When pointV leaves cluster i
                    double removeVCost = -solution.sc[pointV][i] / (2 * sizeI);
                    
                    // When pointU leaves cluster j
                    double removeUCost = -solution.sc[pointU][j] / (2 * sizeJ);
                    
                    // When pointV joins cluster j
                    double addVCost = solution.sc[pointV][j] / (2 * sizeJ);
                    
                    // When pointU joins cluster i
                    double addUCost = solution.sc[pointU][i] / (2 * sizeI);
                    
                    // Total change in solution value
                    double deltaF = removeVCost + removeUCost + addVCost + addUCost;
                    
                    // If improvement found, perform swap and return
                    if (deltaF < 0) {
                        swap(solution, i, pointV, j, pointU, deltaF);
                        return true;
                    }
                }
            }
        }
    }
    
    // No improvement found
    return false;
}

void LocalSearch::swap(Solution& solution, int i, int pointV, int j, int pointU, double df) {
    // Update point assignments
    solution.assignment[pointV] = j;
    solution.assignment[pointU] = i;
    
    // Keep track of if this is a single swap (LocalSearch) or part of multiple swaps (Shaking)
    static int swapCounter = 0;
    swapCounter++;
    
    // Update SC matrix for all points
    for (int p = 0; p < solution.nDataPoints; p++) {
        // Update SC values for cluster i
        // Remove pointV (distance from p to pointV) and add pointU (distance from p to pointU)
        solution.sc[p][i] = solution.sc[p][i] 
                           - solution.distances->getDistance(p, pointV) 
                           + solution.distances->getDistance(p, pointU);
        
        // Update SC values for cluster j
        // Remove pointU (distance from p to pointU) and add pointV (distance from p to pointV)
        solution.sc[p][j] = solution.sc[p][j] 
                           - solution.distances->getDistance(p, pointU) 
                           + solution.distances->getDistance(p, pointV);
    }
    
    // Update solution value with the precalculated delta or recalculate if we've done many swaps
    if (swapCounter <= 5) {
        // For typical LocalSearch use - just incremental update
        solution.solutionValue += df;
    } else {
        // If we've done many swaps, recalculate from scratch to avoid drift
        solution.solutionValue = 0;
        for (int c = 0; c < solution.nClusters; c++) {
            double sumDistances = 0;
            
            for (int p = 0; p < solution.nDataPoints; p++) {
                if (solution.assignment[p] == c) {
                    sumDistances += solution.sc[p][c];
                }
            }
            
            // Add to solution value (divide by 2 because each distance is counted twice)
            solution.solutionValue += sumDistances / (2 * solution.clusterSizes[c]);
        }
        
        // Reset counter after recalculation
        swapCounter = 0;
    }
}

bool LocalSearch::checkSolution(Solution* solutionBefore, Solution* solutionAfter, double deltaSolutionValue) {
    // Check that solution value has been correctly updated
    if (abs(solutionAfter->solutionValue - solutionBefore->solutionValue - deltaSolutionValue) > 1e-6) {
        cout << "Solution value incorrectly updated. Expected: " 
             << (solutionBefore->solutionValue + deltaSolutionValue) 
             << ", Actual: " << solutionAfter->solutionValue << endl;
        return false;
    }
    
    // Check that cluster sizes are maintained (for balanced clustering)
    for (int c = 0; c < solutionBefore->nClusters; c++) {
        int countBefore = 0, countAfter = 0;
        
        for (int p = 0; p < solutionBefore->nDataPoints; p++) {
            if (solutionBefore->assignment[p] == c) countBefore++;
            if (solutionAfter->assignment[p] == c) countAfter++;
        }
        
        if (countBefore != countAfter) {
            cout << "Cluster size balance violated for cluster " << c << endl;
            return false;
        }
    }
    
    // Verify SC matrix consistency
    for (int p = 0; p < solutionAfter->nDataPoints; p++) {
        for (int c = 0; c < solutionAfter->nClusters; c++) {
            double calculatedSC = 0.0;
            
            for (int q = 0; q < solutionAfter->nDataPoints; q++) {
                if (solutionAfter->assignment[q] == c) {
                    calculatedSC += solutionAfter->distances->getDistance(p, q);
                }
            }
            
            if (abs(calculatedSC - solutionAfter->sc[p][c]) > 1e-6) {
                cout << "SC matrix inconsistency at [" << p << "][" << c << "]" << endl;
                return false;
            }
        }
    }
    
    return true;
}