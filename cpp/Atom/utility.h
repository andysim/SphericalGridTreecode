#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include "atom.h"

bool InputVecOfAtomsFromWaterBox(std::vector<Atom> &atoms, const std::string &fileName);

void GetLowerUpperBoundsOnAtomCoords(double &lowestXCoordOfMol, double &highestXCoordOfMol, double &lowestYCoordOfMol,
                                     double &highestYCoordOfMol, double &lowestZCoordOfMol, double &highestZCoordOfMol,
                                     const std::vector<Atom> &atoms);

void AdjustMaxMinOfMolForCubeWidth(double &lowestXCoordOfMol, double &highestXCoordOfMol, double &lowestYCoordOfMol,
                                   double &highestYCoordOfMol, double &lowestZCoordOfMol, double &highestZCoordOfMol,
                                   double cubeWidth);

double GetLargestDistanceInMolecule(const std::vector<Atom> &atoms);

bool GetSmallestDistanceBetweenWatersInMolecule(double &minDistance, const std::vector<Atom> &atoms);
bool OutputVecOfAtoms(const std::vector<Atom> &atoms, const std::string &fileName);

double GetMaxErrorInCoulombPotOverVecOfAtoms(std::vector<Atom> &vecOfAtoms);

double GetRMSErrorInCoulombPotOverVecOfAtoms(std::vector<Atom> &vecOfAtoms);

double GetRMSErrorInGradOfCoulombPotOverVecOfAtoms(std::vector<Atom> &vecOfAtoms);

double GetRelativeRMSErrorOfGradOfCoulombPotOverVecOfAtoms(std::vector<Atom> &vecOfAtoms);

#endif  // UTILITY_H
