#include <cmath>
#include <vector>
#include <iostream>

#include "../ErfPointCharge/erfPointCharge.h"
#include "../Atom/atom.h"
#include "hierarchy.h"

// public procs for Hierarchy

Hierarchy::Hierarchy(double widthOfCubes, double minXCoordInMol, double minYCoordInMol, double minZCoordInMol,
                     double ewaldCoeff, std::size_t numSphericalQuadraturePoints,
                     std::vector<double> &gaussHermiteRoots, std::vector<double> &gaussHermiteWeights,
                     double maxXCoordInMol, double maxYCoordInMol, double maxZCoordInMol,
                     const std::vector<Atom> &vecOfAtoms)
    : m_widthOfCubes(widthOfCubes),
      m_minXCoordInMol(minXCoordInMol),
      m_minYCoordInMol(minYCoordInMol),
      m_minZCoordInMol(minZCoordInMol),
      m_ewaldCoeff(ewaldCoeff),
      m_numSphericalQuadraturePoints(numSphericalQuadraturePoints),
      m_gaussHermiteRoots(gaussHermiteRoots),
      m_gaussHermiteWeights(gaussHermiteWeights) {
    SetupTopIndices(widthOfCubes, minXCoordInMol, minYCoordInMol, minZCoordInMol, maxXCoordInMol, maxYCoordInMol,
                    maxZCoordInMol);

    SetupVecOfVecOfCubes();

    if (!SetupVecOfErfPointCharges(vecOfAtoms)) return;

    std::size_t numRadialQuadraturePts = gaussHermiteRoots.size();
    std::size_t numLevelsInHierarchyOfCubes = GetNumLevelsInHierarchyOfCubes();
    AllocatePointChargeArrays(numLevelsInHierarchyOfCubes, numRadialQuadraturePts);
    m_valid = true;
}

bool Hierarchy::IsValid() const { return m_valid; }

bool Hierarchy::AccumulatePotentialsForNewSphericalQuadraturePoint(

    double sphPtx, double sphPty, double sphPtz) {
    SetNewSphericalQuadraturePoint(sphPtx, sphPty, sphPtz);

    LoadCosineSineBasisFunctions();

    if (!AccumulatePotentialsUpwardFromChildToParent()) {
        std::cout << "AccumulateCubePotentials failed!" << std::endl;
        return false;
    }

    if (!AccumulatePotentialsDownwardFromParentToChild()) {
        std::cout << "AccumulatePotentialsDownwardFromParentToChild failed!" << std::endl;
        return false;
    }

    return true;
}

void Hierarchy::LoadApproxErfPotential_and_Grad() {
    IntegrateAccumPotentials();

    std::vector<ErfPointCharge>::iterator iterP = m_vecOfErfPointCharges.begin();
    std::vector<ErfPointCharge>::iterator endP = m_vecOfErfPointCharges.end();
    for (; iterP != endP; ++iterP) iterP->LoadApproxErfPotential_and_Grad();
}

bool Hierarchy::GetVecsOfApproxErfPot_and_Grads(std::vector<double> &approxErfPots,
                                                std::vector<double> &approxGradErfPot_wrt_X,
                                                std::vector<double> &approxGradErfPot_wrt_Y,
                                                std::vector<double> &approxGradErfPot_wrt_Z,
                                                std::size_t sizeOfListOfAtoms) {
    approxErfPots.resize(sizeOfListOfAtoms, 0.0);
    approxGradErfPot_wrt_X.resize(sizeOfListOfAtoms, 0.0);
    approxGradErfPot_wrt_Y.resize(sizeOfListOfAtoms, 0.0);
    approxGradErfPot_wrt_Z.resize(sizeOfListOfAtoms, 0.0);

    std::vector<unsigned> flag(sizeOfListOfAtoms, 0);

    std::vector<ErfPointCharge>::const_iterator iterP = m_vecOfErfPointCharges.begin();
    std::vector<ErfPointCharge>::const_iterator endP = m_vecOfErfPointCharges.end();
    for (; iterP != endP; ++iterP) {
        std::size_t index = iterP->GetIndexInListOfAtoms();

        if (index >= sizeOfListOfAtoms) {
            std::cout << "index in list of atoms out of range!" << std::endl;
            return false;
        }
        if (flag[index] != 0) {
            std::cout << "duplicate atom index in vec of erfPointCharges!" << std::endl;
            return false;
        }

        double approxErfPotential = 0.0;
        std::vector<double> gradOfApproxErfPotential;
        iterP->GetApproxErfPotential_and_Grad(approxErfPotential, gradOfApproxErfPotential);
        approxErfPots[index] = approxErfPotential;
        approxGradErfPot_wrt_X[index] = gradOfApproxErfPotential[0];
        approxGradErfPot_wrt_Y[index] = gradOfApproxErfPotential[1];
        approxGradErfPot_wrt_Z[index] = gradOfApproxErfPotential[2];
        flag[index] = 1;
    }

    return true;
}

// private procs for Hierarchy

void Hierarchy::IntegrateAccumPotentials() {
    std::size_t numLevelsInHierarchy = GetNumLevelsInHierarchyOfCubes();
    std::vector<ErfPointCharge>::iterator iterP = m_vecOfErfPointCharges.begin();
    std::vector<ErfPointCharge>::iterator endP = m_vecOfErfPointCharges.end();
    for (; iterP != endP; ++iterP) {
        iterP->IntegrateAccumPotentials(m_ewaldCoeff, numLevelsInHierarchy, m_gaussHermiteWeights);
    }
}

void Hierarchy::SetNewSphericalQuadraturePoint(double sphPtx, double sphPty, double sphPtz) {
    m_spherePointX = sphPtx;
    m_spherePointY = sphPty;
    m_spherePointZ = sphPtz;
    // zero fill Accum Potentials
    std::vector<VecOfCubes>::iterator iterV = m_vecOfVecOfCubes.begin();
    std::vector<VecOfCubes>::iterator endV = m_vecOfVecOfCubes.end();
    for (; iterV != endV; ++iterV) iterV->ZeroFillAccumulatedSourcePotentials();
}

void Hierarchy::LoadCosineSineBasisFunctions() {
    std::vector<ErfPointCharge>::iterator iterP = m_vecOfErfPointCharges.begin();
    std::vector<ErfPointCharge>::iterator endP = m_vecOfErfPointCharges.end();
    std::size_t numLevelsInHierarchyOfCubes = GetNumLevelsInHierarchyOfCubes();
    for (; iterP != endP; ++iterP)
        iterP->LoadCosineSineBasisFunctions(m_ewaldCoeff, m_numSphericalQuadraturePoints, numLevelsInHierarchyOfCubes,
                                            m_gaussHermiteRoots, m_spherePointX, m_spherePointY, m_spherePointZ);
}

bool Hierarchy::AccumulatePotentialsUpwardFromChildToParent() {
    std::size_t numLevels = GetNumLevelsInHierarchyOfCubes();
    if (numLevels == 0) {
        return false;
    }
    // first from erfCharges into level 0
    if (!m_vecOfVecOfCubes[0].AccumulatePotentialsFromErfPointChargesIntoCubes(m_numAtomsPerLevelZeroCube,
                                                                               m_vecOfErfPointCharges)) {
        std::cout << "AccumulatePotentialsFromErfPointChargesIntoCubes failed!" << std::endl;
        return false;
    }
    // next accumulate from child to parent
    if (numLevels > 1) {
        std::size_t topChildLevel = numLevels - 2;
        for (std::size_t level = 0; level <= topChildLevel; ++level) {
            std::size_t parentLevel = level + 1;
            std::vector<Cube> vecOfCubesAtParentLevel = m_vecOfVecOfCubes[parentLevel].GetVecOfCubes();
            std::vector<Cube> tmpVec = m_vecOfVecOfCubes[level].GetVecOfCubes();
            if (!m_vecOfVecOfCubes[level].AddUpwardAccumPotentialsToThoseOfParentCubes(vecOfCubesAtParentLevel))
                return false;
            m_vecOfVecOfCubes[parentLevel].SetVecOfCubes(vecOfCubesAtParentLevel);
        }
    }
    return true;
}

bool Hierarchy::AccumulatePotentialsDownwardFromParentToChild() {
    if (GetNumLevelsInHierarchyOfCubes() == 0) {
        return false;
    }

    std::size_t topLevel = GetNumLevelsInHierarchyOfCubes() - 1;
    if (!m_vecOfVecOfCubes[topLevel].SumUpwardAccumPotsAtTopLevelAndPutIntoDownwardAccumPots()) {
        std::cout << "SumAccumPotentialsOverCubesAndPutBackIntoCubes failed!" << std::endl;
        return false;
    }

    if (topLevel >= 1) {
        for (std::size_t j = 0; j < topLevel; ++j) {
            std::size_t parLevel = topLevel - j;
            std::size_t childLevel = parLevel - 1;
            if (!m_vecOfVecOfCubes[parLevel].AdjustDownwardAccumPotentialsByNearestNeighbors()) return false;
            std::vector<Cube> vecOfCubesAtParentLevel = m_vecOfVecOfCubes[parLevel].GetVecOfCubes();
            if (!m_vecOfVecOfCubes[childLevel].GetDownwardAccumPotentialsFromThoseOfParentCubes(
                    vecOfCubesAtParentLevel))
                return false;
        }
    }

    if (!m_vecOfVecOfCubes[0].AddDownwardAccumCubePotentialsToErfPointCharges(
            m_vecOfErfPointCharges, m_numAtomsPerLevelZeroCube, m_spherePointX, m_spherePointY, m_spherePointZ)) {
        std::cout << "AddAccumCubePotentialsToErfPointCharges fails!" << std::endl;
        return false;
    }

    return true;
}

std::size_t Hierarchy::GetNumLevelsInHierarchyOfCubes() const { return m_topIndexX.size(); }

std::size_t Hierarchy::SetupTopIndicesAtLevel(std::size_t &indexX, std::size_t &indexY, std::size_t &indexZ,
                                              double xWidthOfMol, double yWidthOfMol, double zWidthOfMol,
                                              double cubeWidthAtThisLevel) {
    double eps = 1.0e-14;

    if (cubeWidthAtThisLevel >= xWidthOfMol) {
        indexX = 1;
    } else {
        double numCubesX = floor(xWidthOfMol / cubeWidthAtThisLevel);
        if (xWidthOfMol - numCubesX * cubeWidthAtThisLevel > eps) numCubesX += 1.0;
        indexX = static_cast<unsigned>(numCubesX);
    }

    if (cubeWidthAtThisLevel >= yWidthOfMol) {
        indexY = 1;
    } else {
        double numCubesY = floor(yWidthOfMol / cubeWidthAtThisLevel);
        if (yWidthOfMol - numCubesY * cubeWidthAtThisLevel > eps) numCubesY += 1.0;
        indexY = static_cast<unsigned>(numCubesY);
    }

    if (cubeWidthAtThisLevel >= zWidthOfMol) {
        indexZ = 1;
    } else {
        double numCubesZ = floor(zWidthOfMol / cubeWidthAtThisLevel);
        if (zWidthOfMol - numCubesZ * cubeWidthAtThisLevel > eps) numCubesZ += 1.0;
        indexZ = static_cast<unsigned>(numCubesZ);
    }
    std::size_t maxIndex = 0;
    if (indexX > maxIndex) maxIndex = indexX;
    if (indexY > maxIndex) maxIndex = indexY;
    if (indexX > maxIndex) maxIndex = indexZ;

    return maxIndex;
}

void Hierarchy::SetupTopIndices(double widthOfCubes, double minXCoordInMol, double minYCoordInMol,
                                double minZCoordInMol, double maxXCoordInMol, double maxYCoordInMol,
                                double maxZCoordInMol) {
    double xWidthOfMol = maxXCoordInMol - minXCoordInMol;
    double yWidthOfMol = maxYCoordInMol - minYCoordInMol;
    double zWidthOfMol = maxZCoordInMol - minZCoordInMol;

    // start at level 0
    double cubeWidthAtThisLevel = widthOfCubes;
    std::size_t indexX = 0;
    std::size_t indexY = 0;
    std::size_t indexZ = 0;
    std::size_t maxIndex =
        SetupTopIndicesAtLevel(indexX, indexY, indexZ, xWidthOfMol, yWidthOfMol, zWidthOfMol, cubeWidthAtThisLevel);
    m_topIndexX.push_back(indexX);
    m_topIndexY.push_back(indexY);
    m_topIndexZ.push_back(indexZ);

    // continue as long as > 4x4x4 situation (quadrature is good to 4*sqrt(3)*widthOfCubes
    // at bottom level ewald ---
    while (maxIndex > 4) {
        cubeWidthAtThisLevel *= 2.0;
        maxIndex =
            SetupTopIndicesAtLevel(indexX, indexY, indexZ, xWidthOfMol, yWidthOfMol, zWidthOfMol, cubeWidthAtThisLevel);

        m_topIndexX.push_back(indexX);
        m_topIndexY.push_back(indexY);
        m_topIndexZ.push_back(indexZ);
    }
}

void Hierarchy::SetupVecOfVecOfCubes() {
    std::size_t numLevelsInHierarchyOfCubes = m_topIndexX.size();
    std::size_t numRadialQuadraturePts = m_gaussHermiteRoots.size();

    m_vecOfVecOfCubes.clear();
    m_vecOfVecOfCubes.reserve(numLevelsInHierarchyOfCubes);
    double widthOfCubes = m_widthOfCubes;
    for (std::size_t level = 0; level < numLevelsInHierarchyOfCubes; ++level) {
        VecOfCubes vecOfCubes(widthOfCubes, m_minXCoordInMol, m_minYCoordInMol, m_minZCoordInMol, m_topIndexX,
                              m_topIndexY, m_topIndexZ, level, numRadialQuadraturePts);
        m_vecOfVecOfCubes.push_back(vecOfCubes);
    }
}

std::size_t Hierarchy::GetIndexInZeroLevelVecOfCubes(std::size_t indexX, std::size_t indexY, std::size_t indexZ) const {
    return indexZ * m_topIndexX[0] * m_topIndexY[0] + indexY * m_topIndexX[0] + indexX;
}

bool Hierarchy::SetIndexInZeroLevelVecOfCubesToAtoms(std::vector<Atom> &tmpVec, const std::vector<Atom> &oldVec) const {
    tmpVec.clear();
    tmpVec.reserve(oldVec.size());
    std::size_t numCubes = m_topIndexX[0] * m_topIndexY[0] * m_topIndexZ[0];
    for (auto atom : oldVec) {
        double xCoord = atom.GetXCoord();
        double yCoord = atom.GetYCoord();
        double zCoord = atom.GetZCoord();
        double relX = floor((xCoord - m_minXCoordInMol) / m_widthOfCubes);
        double relY = floor((yCoord - m_minYCoordInMol) / m_widthOfCubes);
        double relZ = floor((zCoord - m_minZCoordInMol) / m_widthOfCubes);
        if ((relX < 0.0) || (relY < 0.0) || (relZ < 0.0)) {
            std::cout << "relX,Y,Z = " << relX << ", " << relY << ", " << relZ << std::endl;
            return false;
        }
        std::size_t indexX = static_cast<std::size_t>(relX);
        std::size_t indexY = static_cast<std::size_t>(relY);
        std::size_t indexZ = static_cast<std::size_t>(relZ);
        if ((indexX >= m_topIndexX[0]) || (indexY >= m_topIndexY[0]) || (indexZ >= m_topIndexZ[0])) {
            std::cout << "relX,Y,Z = " << relX << ", " << relY << ", " << relZ << std::endl;
            return false;
        }

        std::size_t indexInListOfCubes = GetIndexInZeroLevelVecOfCubes(indexX, indexY, indexZ);
        if (indexInListOfCubes >= numCubes) {
            std::cout << "some atom cube index out of range of cubes!" << std::endl;
            return false;
        }
        atom.SetIndexInListOfCubes(indexInListOfCubes);
        tmpVec.push_back(atom);
    }

    return true;
}

bool Hierarchy::SortVecOfAtomsByCubeIndex(std::vector<Atom> &newVec, const std::vector<Atom> &oldVec) {
    std::vector<std::size_t> numAtomsPerCube;
    std::size_t numCubes = m_topIndexX[0] * m_topIndexY[0] * m_topIndexZ[0];
    numAtomsPerCube.reserve(numCubes);

    newVec.clear();
    newVec.reserve(oldVec.size());

    std::vector<Atom> tmpVec;
    if (!SetIndexInZeroLevelVecOfCubesToAtoms(tmpVec, oldVec)) {
        std::cout << "AtomsToCubes failed!" << std::endl;
        return false;
    }

    std::vector<std::vector<Atom>> atomsByCubes;
    std::vector<Atom> emptyVec;
    for (std::size_t n = 0; n < numCubes; ++n) atomsByCubes.push_back(emptyVec);

    for (auto atom : tmpVec) {
        std::size_t index = atom.GetIndexInListOfCubes();
        atomsByCubes.at(index).push_back(atom);
    }

    for (auto atomVec : atomsByCubes) {
        numAtomsPerCube.push_back(atomVec.size());
        std::copy(atomVec.begin(), atomVec.end(), std::back_inserter(newVec));
    }
    // check ordering

    std::size_t offSet = 0;
    for (std::size_t n = 0; n < numCubes; ++n) {
        for (std::size_t i = 0; i < numAtomsPerCube[n]; ++i) {
            if (newVec[offSet + i].GetIndexInListOfCubes() != n) {
                std::cout << "sorted newVec order is wrong!" << std::endl;
                return false;
            }
        }
        offSet += numAtomsPerCube[n];
    }
    // check sizes and total
    if (numAtomsPerCube.size() != numCubes) {
        std::cout << "numAtomsPerCube.size() != numCubes" << std::endl;
        return false;
    }
    std::size_t sum = 0;
    for (std::size_t n = 0; n < numCubes; ++n) sum += numAtomsPerCube[n];

    if (sum != newVec.size()) {
        std::cout << "total num atoms in cubes doesn't match vecOfAtoms.size()" << std::endl;
        return false;
    }

    m_numAtomsPerLevelZeroCube = numAtomsPerCube;

    return true;
}

void Hierarchy::AtomsToErfPointCharges(std::vector<ErfPointCharge> &vecOfErfPointCharges,
                                       const std::vector<Atom> &vecOfAtoms) const {
    vecOfErfPointCharges.clear();
    vecOfErfPointCharges.reserve(vecOfAtoms.size());

    std::vector<Atom>::const_iterator iterA = vecOfAtoms.begin();
    std::vector<Atom>::const_iterator endA = vecOfAtoms.end();
    for (; iterA != endA; ++iterA) {
        ErfPointCharge q(iterA->GetCharge(), iterA->GetXCoord(), iterA->GetYCoord(), iterA->GetZCoord(),
                         iterA->GetIndexInListOfAtoms());
        vecOfErfPointCharges.push_back(q);
    }
}

bool Hierarchy::SetupVecOfErfPointCharges(const std::vector<Atom> &vecOfAtoms) {
    std::vector<Atom> sortedVecOfAtoms;
    if (!SortVecOfAtomsByCubeIndex(sortedVecOfAtoms, vecOfAtoms)) {
        std::cout << "SortVecOfAtomsByCubes failed!" << std::endl;
        return false;
    }
    std::vector<ErfPointCharge> vecOfErfPointCharges;
    AtomsToErfPointCharges(vecOfErfPointCharges, sortedVecOfAtoms);
    m_vecOfErfPointCharges = vecOfErfPointCharges;
    return true;
}

void Hierarchy::AllocatePointChargeArrays(std::size_t numLevelsInHierarchyOfCubes, std::size_t numRadialQuadPts) {
    std::vector<ErfPointCharge>::iterator iterP = m_vecOfErfPointCharges.begin();
    std::vector<ErfPointCharge>::iterator endP = m_vecOfErfPointCharges.end();
    for (; iterP != endP; ++iterP) iterP->AllocateArrays(numLevelsInHierarchyOfCubes, numRadialQuadPts);
}
