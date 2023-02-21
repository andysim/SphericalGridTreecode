#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <math.h>

#include "atom.h"
#include "utility.h"

#include "../IO/readInputFiles.h"
#include "../Quadrature/hermiteRootsWeights.h"

bool InputVecOfAtomsFromWaterBox(std::vector<Atom> &atoms, const std::string &fileName) {
    atoms.clear();
    std::vector<std::string> names;
    std::vector<double> charges;
    std::vector<double> xCoords;
    std::vector<double> yCoords;
    std::vector<double> zCoords;
    std::vector<std::size_t> indicesInVecOfAtoms;
    std::vector<std::vector<std::size_t>> exclusionLists;

    if (!InputAtomInfoFromWaterBox(names, charges, xCoords, yCoords, zCoords, indicesInVecOfAtoms, exclusionLists,
                                   fileName))
        return false;

    atoms.reserve(names.size());
    std::vector<std::string>::const_iterator iterN = names.begin();
    std::vector<std::string>::const_iterator endN = names.end();
    std::vector<double>::const_iterator iterC = charges.begin();
    std::vector<double>::const_iterator endC = charges.end();
    std::vector<double>::const_iterator iterX = xCoords.begin();
    std::vector<double>::const_iterator endX = xCoords.end();
    std::vector<double>::const_iterator iterY = yCoords.begin();
    std::vector<double>::const_iterator endY = yCoords.end();
    std::vector<double>::const_iterator iterZ = zCoords.begin();
    std::vector<double>::const_iterator endZ = zCoords.end();
    std::vector<std::size_t>::const_iterator iterI = indicesInVecOfAtoms.begin();
    std::vector<std::size_t>::const_iterator endI = indicesInVecOfAtoms.end();
    std::vector<std::vector<std::size_t>>::const_iterator iterE = exclusionLists.begin();
    std::vector<std::vector<std::size_t>>::const_iterator endE = exclusionLists.end();

    for (; (iterN != endN) && (iterC != endC) && (iterX != endX) && (iterY != endY) && (iterZ != endZ) &&
           (iterI != endI) && (iterE != endE);
         ++iterN, ++iterC, ++iterX, ++iterY, ++iterZ, ++iterI, ++iterE) {
        Atom atom(*iterN, *iterC, *iterX, *iterY, *iterZ, *iterI, *iterE);
        atoms.push_back(atom);
    }

    return true;
}

bool OutputVecOfAtoms(const std::vector<Atom> &atoms, const std::string &fileName) {
    std::ofstream ofs;
    ofs.open(fileName);
    if (!ofs.is_open()) return false;

    ofs << atoms.size() << "  Waters by translation" << std::endl;
    for (std::size_t n = 0; n < atoms.size(); ++n) {
        ofs << n + 1 << "  " << atoms[n].GetName() << "  " << atoms[n].GetXCoord() << "  " << atoms[n].GetYCoord()
            << "  " << atoms[n].GetZCoord() << std::endl;
    }
    ofs.close();

    return true;
}

void GetLowerUpperBoundsOnAtomCoords(double &lowestXCoordOfMol, double &highestXCoordOfMol, double &lowestYCoordOfMol,
                                     double &highestYCoordOfMol, double &lowestZCoordOfMol, double &highestZCoordOfMol,
                                     const std::vector<Atom> &atoms) {
    lowestXCoordOfMol = std::numeric_limits<double>::max();
    lowestYCoordOfMol = std::numeric_limits<double>::max();
    lowestZCoordOfMol = std::numeric_limits<double>::max();
    highestXCoordOfMol = -std::numeric_limits<double>::max();
    highestYCoordOfMol = -std::numeric_limits<double>::max();
    highestZCoordOfMol = -std::numeric_limits<double>::max();

    std::vector<Atom>::const_iterator iterA = atoms.begin();
    std::vector<Atom>::const_iterator endA = atoms.end();
    for (; iterA != endA; ++iterA) {
        double xCoord = iterA->GetXCoord();
        double yCoord = iterA->GetYCoord();
        double zCoord = iterA->GetZCoord();

        if (xCoord > highestXCoordOfMol) highestXCoordOfMol = xCoord;
        if (yCoord > highestYCoordOfMol) highestYCoordOfMol = yCoord;
        if (zCoord > highestZCoordOfMol) highestZCoordOfMol = zCoord;
        if (xCoord < lowestXCoordOfMol) lowestXCoordOfMol = xCoord;
        if (yCoord < lowestYCoordOfMol) lowestYCoordOfMol = yCoord;
        if (zCoord < lowestZCoordOfMol) lowestZCoordOfMol = zCoord;
    }
}

void AdjustMaxMinOfMolForCubeWidth(double &lowestXCoordOfMol, double &highestXCoordOfMol, double &lowestYCoordOfMol,
                                   double &highestYCoordOfMol, double &lowestZCoordOfMol, double &highestZCoordOfMol,
                                   double cubeWidth) {
    double xWidthOfMol = highestXCoordOfMol - lowestXCoordOfMol;
    double numXCubes = floor(xWidthOfMol / cubeWidth);
    if ((numXCubes < 1.0) || (numXCubes * cubeWidth < xWidthOfMol)) numXCubes += 1.0;
    double delXWidth = 0.5 * (numXCubes * cubeWidth - xWidthOfMol);
    lowestXCoordOfMol -= delXWidth;
    highestXCoordOfMol += delXWidth;

    double yWidthOfMol = highestYCoordOfMol - lowestYCoordOfMol;
    double numYCubes = floor(yWidthOfMol / cubeWidth);
    if ((numYCubes < 1.0) || (numYCubes * cubeWidth < yWidthOfMol)) numYCubes += 1.0;
    double delYWidth = 0.5 * (numYCubes * cubeWidth - yWidthOfMol);
    lowestYCoordOfMol -= delYWidth;
    highestYCoordOfMol += delYWidth;

    double zWidthOfMol = highestZCoordOfMol - lowestZCoordOfMol;
    double numZCubes = floor(zWidthOfMol / cubeWidth);
    if ((numZCubes < 1.0) || (numZCubes * cubeWidth < zWidthOfMol)) numZCubes += 1.0;
    double delZWidth = 0.5 * (numZCubes * cubeWidth - zWidthOfMol);
    lowestZCoordOfMol -= delZWidth;
    highestZCoordOfMol += delZWidth;
}

double GetLargestDistanceInMolecule(const std::vector<Atom> &atoms) {
    double maxDistance = 0.0;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        double xCoordI = atoms[i].GetXCoord();
        double yCoordI = atoms[i].GetYCoord();
        double zCoordI = atoms[i].GetZCoord();
        for (std::size_t j = i + 1; j < atoms.size(); ++j) {
            double xCoordJ = atoms[j].GetXCoord();
            double yCoordJ = atoms[j].GetYCoord();
            double zCoordJ = atoms[j].GetZCoord();
            double distance =
                sqrt((xCoordJ - xCoordI) * (xCoordJ - xCoordI) + (yCoordJ - yCoordI) * (yCoordJ - yCoordI) +
                     (zCoordJ - zCoordI) * (zCoordJ - zCoordI));
            if (distance > maxDistance) maxDistance = distance;
        }
    }
    return maxDistance;
}

bool GetSmallestDistanceBetweenTwoWatersInMol(double &distance, const std::vector<Atom> &atoms, std::size_t n,
                                              std::size_t m) {
    distance = std::numeric_limits<double>::max();
    std::size_t numWaters = atoms.size() / 3;

    if ((n >= numWaters) || (m >= numWaters) || (n == m)) return false;

    for (std::size_t i = 0; i < 3; ++i) {
        std::size_t ind = 3 * n + i;
        double x1 = atoms[ind].GetXCoord();
        double y1 = atoms[ind].GetYCoord();
        double z1 = atoms[ind].GetZCoord();
        for (std::size_t j = 0; j < 3; ++j) {
            std::size_t jnd = 3 * m + j;
            double x2 = atoms[jnd].GetXCoord();
            double y2 = atoms[jnd].GetYCoord();
            double z2 = atoms[jnd].GetZCoord();
            double dist = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
            if (dist < distance) distance = dist;
        }
    }
    return true;
}

bool GetSmallestDistanceBetweenWatersInMolecule(double &minDistance, const std::vector<Atom> &atoms) {
    minDistance = std::numeric_limits<double>::max();

    std::size_t numWaters = atoms.size() / 3;

    for (std::size_t n = 0; n < numWaters; ++n) {
        for (std::size_t m = n + 1; m < numWaters; ++m) {
            double distance = std::numeric_limits<double>::max();
            if (!GetSmallestDistanceBetweenTwoWatersInMol(distance, atoms, n, m)) return false;
            if (distance < minDistance) minDistance = distance;
        }
    }
    return true;
}

double GetMaxErrorInCoulombPotOverVecOfAtoms(std::vector<Atom> &vecOfAtoms) {
    double maxError = 0.0;
    for (std::size_t i = 0; i < vecOfAtoms.size(); ++i) {
        double absErrorInApproxCoulombPotential = 0.0;
        double sumSqErrorInGradOfApproxCoulombPotential;
        vecOfAtoms[i].GetErrorInApproxCoulombPotential_and_Grad(absErrorInApproxCoulombPotential,
                                                                sumSqErrorInGradOfApproxCoulombPotential);
        if (absErrorInApproxCoulombPotential > maxError) maxError = absErrorInApproxCoulombPotential;
    }
    return maxError;
}

double GetRMSErrorInCoulombPotOverVecOfAtoms(std::vector<Atom> &vecOfAtoms) {
    double msqError = 0.0;
    for (std::size_t i = 0; i < vecOfAtoms.size(); ++i) {
        double absErrorInApproxCoulombPotential = 0.0;
        double sumSqErrorInGradOfApproxCoulombPotential;
        vecOfAtoms[i].GetErrorInApproxCoulombPotential_and_Grad(absErrorInApproxCoulombPotential,
                                                                sumSqErrorInGradOfApproxCoulombPotential);
        msqError += absErrorInApproxCoulombPotential * absErrorInApproxCoulombPotential;
    }
    msqError /= static_cast<double>(vecOfAtoms.size());
    return sqrt(msqError);
}

double GetRMSErrorInGradOfCoulombPotOverVecOfAtoms(std::vector<Atom> &vecOfAtoms) {
    double msqError = 0.0;
    for (std::size_t i = 0; i < vecOfAtoms.size(); ++i) {
        double absErrorInApproxCoulombPotential = 0.0;
        double sumSqErrorInGradOfApproxCoulombPotential;
        vecOfAtoms[i].GetErrorInApproxCoulombPotential_and_Grad(absErrorInApproxCoulombPotential,
                                                                sumSqErrorInGradOfApproxCoulombPotential);
        msqError += sumSqErrorInGradOfApproxCoulombPotential;
    }
    msqError /= (3.0 * static_cast<double>(vecOfAtoms.size()));
    return sqrt(msqError);
}

double GetRelativeRMSErrorOfGradOfCoulombPotOverVecOfAtoms(std::vector<Atom> &vecOfAtoms) {
    double rmsError = GetRMSErrorInGradOfCoulombPotOverVecOfAtoms(vecOfAtoms);
    double euclideanLengthVec = 0.0;
    for (std::size_t i = 0; i < vecOfAtoms.size(); ++i) {
        double exactCoulombPotential;
        std::vector<double> gradOfExactCoulombPotential;
        vecOfAtoms[i].GetExactCoulombPotential_and_Grad(exactCoulombPotential, gradOfExactCoulombPotential);
        for (std::size_t j = 0; j < 3; ++j)
            euclideanLengthVec += gradOfExactCoulombPotential[j] * gradOfExactCoulombPotential[j];
    }
    euclideanLengthVec /= (3.0 * static_cast<double>(vecOfAtoms.size()));
    return rmsError / sqrt(euclideanLengthVec);
}
