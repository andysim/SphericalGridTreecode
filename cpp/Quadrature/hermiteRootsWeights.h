#ifndef HERMITE_ROOTS_WEIGHTS_H
#define HERMITE_ROOTS_WEIGHTS_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

double T(double t);

void GetAiryRoots(std::vector<double> &airyRoots, const std::vector<double> &first11AiryRoots, std::size_t m);

void GetXinit_Gatteschi(std::vector<double> &xinitAiry, const std::vector<double> &airyRoots, double a, double nu);

void GetXinit_Tricomi(std::vector<double> &xinitSin, std::size_t m, double a, double nu);

bool GetInitialGuessesForHermiteRoots(std::vector<double> &initialGuesses, std::size_t n);

void GetHermitePolynomialValByRecursion(double &herm, double &derivHerm, std::size_t n, double x);

void GetRoot_and_WeightOfHermite(double &root, double &weight, std::size_t n, double guessRoot);

void GetRoots_and_WeightsOfHermite(std::vector<double> &roots, std::vector<double> &weights, std::size_t n,
                                   const std::vector<double> &initialGuesses);

double GetEwaldCoeff(double THRESH, double cutoff);

void GetNandNprime(std::size_t &N, std::size_t &Nprime, double ewaldCoeff, double THRESH, double cutoff);

bool GetTrimmedGaussHermiteRoots_and_Weights(std::vector<double> &trimmedRoots, std::vector<double> &trimmedWeights,
                                             double &ewaldCoeff, double THRESH, double cutoff);

unsigned GetDegreeOfSphericalDesignFromHermiteRootsWeights(double THRESH, double directSumCutoff, double ewaldCoeff,
                                                           const std::vector<double> &gaussHermiteRoots,
                                                           const std::vector<double> &gaussHermiteWeights);

#endif  // HERMITE_ROOTS_WEIGHTS_H
