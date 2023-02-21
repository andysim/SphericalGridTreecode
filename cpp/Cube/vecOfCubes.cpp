#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "cube.h"
#include "vecOfCubes.h"
#include "../Atom/atom.h"

// procs for class VecOfCubes
//
VecOfCubes::VecOfCubes(double widthOfCubes, double minXCoordInVecOfCubes, double minYCoordInVecOfCubes,
                       double minZCoordInVecOfCubes, const std::vector<std::size_t> &topXIndexPerLevel,
                       const std::vector<std::size_t> &topYIndexPerLevel,
                       const std::vector<std::size_t> &topZIndexPerLevel, std::size_t levelInHeirarchyOfCubes,
                       std::size_t numRadialQuadraturePts)
    : m_widthOfCubes(widthOfCubes),
      m_topXIndexPerLevel(topXIndexPerLevel),
      m_topYIndexPerLevel(topYIndexPerLevel),
      m_topZIndexPerLevel(topZIndexPerLevel),
      m_minXCoordInVecOfCubes(minXCoordInVecOfCubes),
      m_minYCoordInVecOfCubes(minYCoordInVecOfCubes),
      m_minZCoordInVecOfCubes(minZCoordInVecOfCubes),
      m_levelInHeirarchyOfCubes(levelInHeirarchyOfCubes),
      m_numRadialQuadraturePts(numRadialQuadraturePts) {
    SetupVecOfCubes();
}

void VecOfCubes::ZeroFillAccumulatedSourcePotentials() {
    std::size_t numCubes = GetNumCubesAtMyLevel();
    for (std::size_t n = 0; n < numCubes; ++n) m_vecOfCubes[n].ZeroFillAccumulatedSourcePotentials();
}

bool VecOfCubes::AccumulatePotentialsFromErfPointChargesIntoCubes(
    const std::vector<std::size_t> &numAtomsPerLevelZeroCube, const std::vector<ErfPointCharge> &vecOfErfPointCharges) {
    if (GetMyLevelInHeirarchyOfCubes() != 0) return false;

    std::size_t numElements = GetNumLevelsInHeirarchyOfCubes() * GetNumRadialQuadraturePoints();
    std::size_t offSet = 0;
    std::size_t numCubes = GetNumCubesAtMyLevel();
    for (std::size_t n = 0; n < numCubes; ++n) {
        std::vector<double> cosPot(numElements, 0.0);
        std::vector<double> sinPot(numElements, 0.0);

        for (std::size_t j = offSet; j < offSet + numAtomsPerLevelZeroCube[n]; ++j) {
            if (!vecOfErfPointCharges[j].AddMySourcePotentialsToAccumulatedSourcePotentials(cosPot, sinPot))
                return false;
        }

        m_vecOfCubes[n].SetUpwardsAccumCosPotentials(cosPot);
        m_vecOfCubes[n].SetUpwardsAccumSinPotentials(sinPot);
        offSet += numAtomsPerLevelZeroCube[n];
    }
    return true;
}

bool VecOfCubes::AddUpwardAccumPotentialsToThoseOfParentCubes(std::vector<Cube> &vecOfCubesAtParentLevel) {
    std::size_t myLevel = GetMyLevelInHeirarchyOfCubes();
    std::size_t parentLevel = myLevel + 1;

    if (parentLevel >= GetNumLevelsInHeirarchyOfCubes()) return false;  // no parent level

    for (std::size_t indexZ = 0; indexZ < GetTopZIndexAtMyLevel(); ++indexZ) {
        for (std::size_t indexY = 0; indexY < GetTopYIndexAtMyLevel(); ++indexY) {
            for (std::size_t indexX = 0; indexX < GetTopXIndexAtMyLevel(); ++indexX) {
                std::size_t index = GetIndexInVecOfCubesAtLevel(indexX, indexY, indexZ, myLevel);
                std::size_t parentIndex = GetIndexOfParentOfCube(myLevel, indexX, indexY, indexZ);

                // accum cos arrays
                std::vector<double> fromCos = m_vecOfCubes[index].GetUpwardsAccumCosPotentials();
                std::vector<double>::const_iterator iterFC = fromCos.begin();
                std::vector<double>::const_iterator endFC = fromCos.end();

                std::vector<double> toCos = vecOfCubesAtParentLevel[parentIndex].GetUpwardsAccumCosPotentials();
                std::vector<double>::iterator iterTC = toCos.begin();
                std::vector<double>::iterator endTC = toCos.end();
                for (; (iterFC != endFC) && (iterTC != endTC); ++iterFC, ++iterTC) {
                    *iterTC += *iterFC;
                }

                vecOfCubesAtParentLevel[parentIndex].SetUpwardsAccumCosPotentials(toCos);

                // accum sin arrays
                std::vector<double> fromSin = m_vecOfCubes[index].GetUpwardsAccumSinPotentials();
                std::vector<double>::const_iterator iterFS = fromSin.begin();
                std::vector<double>::const_iterator endFS = fromSin.end();

                std::vector<double> toSin = vecOfCubesAtParentLevel[parentIndex].GetUpwardsAccumSinPotentials();
                std::vector<double>::iterator iterTS = toSin.begin();
                std::vector<double>::iterator endTS = toSin.end();
                for (; (iterFS != endFS) && (iterTS != endTS); ++iterFS, ++iterTS) {
                    *iterTS += *iterFS;
                }
                vecOfCubesAtParentLevel[parentIndex].SetUpwardsAccumSinPotentials(toSin);
            }
        }
    }

    return true;
}

bool VecOfCubes::SumUpwardAccumPotsAtTopLevelAndPutIntoDownwardAccumPots() {
    if (GetMyLevelInHeirarchyOfCubes() + 1 != GetNumLevelsInHeirarchyOfCubes()) return false;

    std::size_t numElements = GetNumLevelsInHeirarchyOfCubes() * GetNumRadialQuadraturePoints();

    std::vector<double> cosPot(numElements, 0.0);
    std::vector<double> sinPot(numElements, 0.0);

    std::size_t numElementsInRow = GetNumRadialQuadraturePoints();
    for (auto cube : m_vecOfCubes) {
        std::vector<double> cosPotCube = cube.GetUpwardsAccumCosPotentials();
        if (cosPotCube.size() != numElements) return false;
        std::vector<double> sinPotCube = cube.GetUpwardsAccumSinPotentials();
        if (sinPotCube.size() != numElements) return false;

        for (std::size_t j = 0; j < numElementsInRow; ++j) {
            cosPot[j] += cosPotCube[j];
            sinPot[j] += sinPotCube[j];
        }
    }

    std::vector<Cube>::iterator iterV = m_vecOfCubes.begin();
    std::vector<Cube>::iterator endV = m_vecOfCubes.end();
    for (; iterV != endV; ++iterV) {
        iterV->SetDownwardsAccumCosPotentials(cosPot);
        iterV->SetDownwardsAccumSinPotentials(sinPot);
    }

    return true;
}

bool VecOfCubes::AdjustDownwardAccumPotentialsByNearestNeighbors() {
    if (GetMyLevelInHeirarchyOfCubes() == 0) return false;  // level 0 doesn't do this adjustment

    std::size_t topLevel = GetNumLevelsInHeirarchyOfCubes() - 1;
    std::size_t myLevel = GetMyLevelInHeirarchyOfCubes();
    std::size_t rowNum = topLevel - myLevel;
    std::size_t nextRowNum = rowNum + 1;
    std::size_t numElementsInRow = GetNumRadialQuadraturePoints();

    for (std::size_t indexZ = 0; indexZ < GetTopZIndexAtMyLevel(); ++indexZ) {
        for (std::size_t indexY = 0; indexY < GetTopYIndexAtMyLevel(); ++indexY) {
            for (std::size_t indexX = 0; indexX < GetTopXIndexAtMyLevel(); ++indexX) {
                std::vector<double> sumNghbrCosPot(2 * numElementsInRow, 0.0);
                std::vector<double> sumNghbrSinPot(2 * numElementsInRow, 0.0);
                std::vector<std::size_t> nearestNeighborList;
                if (!LoadIndicesOfNearestNeighborsOfCube(nearestNeighborList, indexX, indexY, indexZ, myLevel))
                    return false;
                std::size_t offSet = rowNum * numElementsInRow;
                std::size_t nextOffSet = offSet + numElementsInRow;
                for (auto nghbr : nearestNeighborList) {
                    std::vector<double> cosPot = m_vecOfCubes[nghbr].GetUpwardsAccumCosPotentials();
                    std::vector<double> sinPot = m_vecOfCubes[nghbr].GetUpwardsAccumSinPotentials();
                    for (std::size_t j = 0; j < 2 * numElementsInRow; ++j) {
                        sumNghbrCosPot[j] += cosPot[offSet + j];
                        sumNghbrSinPot[j] += sinPot[offSet + j];
                    }
                }
                std::size_t index = GetIndexInVecOfCubesAtLevel(indexX, indexY, indexZ, myLevel);

                std::vector<double> dCosPot = m_vecOfCubes[index].GetDownwardsAccumCosPotentials();
                std::vector<double> dSinPot = m_vecOfCubes[index].GetDownwardsAccumSinPotentials();

                for (std::size_t j = 0; j < numElementsInRow; ++j) {
                    dCosPot[offSet + j] -= sumNghbrCosPot[j];  // remove nghbr interactions at this level
                    dCosPot[nextOffSet + j] += sumNghbrCosPot[numElementsInRow + j];  // add at next level
                    dSinPot[offSet + j] -= sumNghbrSinPot[j];
                    dSinPot[nextOffSet + j] += sumNghbrSinPot[numElementsInRow + j];
                }
                m_vecOfCubes[index].SetDownwardsAccumCosPotentials(dCosPot);
                m_vecOfCubes[index].SetDownwardsAccumSinPotentials(dSinPot);
            }
        }
    }

    return true;
}

bool VecOfCubes::GetDownwardAccumPotentialsFromThoseOfParentCubes(std::vector<Cube> &vecOfCubesAtParentLevel) {
    std::size_t myLevel = GetMyLevelInHeirarchyOfCubes();
    std::size_t parentLevel = myLevel + 1;

    if (parentLevel >= GetNumLevelsInHeirarchyOfCubes()) return false;  // no parent level

    for (std::size_t indexZ = 0; indexZ < GetTopZIndexAtMyLevel(); ++indexZ) {
        for (std::size_t indexY = 0; indexY < GetTopYIndexAtMyLevel(); ++indexY) {
            for (std::size_t indexX = 0; indexX < GetTopXIndexAtMyLevel(); ++indexX) {
                std::size_t index = GetIndexInVecOfCubesAtLevel(indexX, indexY, indexZ, myLevel);
                std::size_t parentIndex = GetIndexOfParentOfCube(myLevel, indexX, indexY, indexZ);
                std::vector<double> parentCosPot =
                    vecOfCubesAtParentLevel[parentIndex].GetDownwardsAccumCosPotentials();
                std::vector<double> parentSinPot =
                    vecOfCubesAtParentLevel[parentIndex].GetDownwardsAccumSinPotentials();
                m_vecOfCubes[index].SetDownwardsAccumCosPotentials(parentCosPot);
                m_vecOfCubes[index].SetDownwardsAccumSinPotentials(parentSinPot);
            }
        }
    }

    return true;
}

bool VecOfCubes::AddDownwardAccumCubePotentialsToErfPointCharges(
    std::vector<ErfPointCharge> &vecOfErfPointCharges, const std::vector<std::size_t> &numAtomsPerLevelZeroCube,
    double sphPtx, double sphPty, double sphPtz) {
    if (GetMyLevelInHeirarchyOfCubes() != 0) return false;

    std::size_t offSet = 0;
    std::size_t numCubes = GetNumCubesAtMyLevel();
    for (std::size_t n = 0; n < numCubes; ++n) {
        std::vector<double> cosPot = m_vecOfCubes[n].GetDownwardsAccumCosPotentials();
        std::vector<double> sinPot = m_vecOfCubes[n].GetDownwardsAccumSinPotentials();
        for (std::size_t j = offSet; j < offSet + numAtomsPerLevelZeroCube[n]; ++j) {
            if (!vecOfErfPointCharges[j].AddAccumCosSinPotentialsToMyAccumCosSinPotentials(cosPot, sinPot, sphPtx,
                                                                                           sphPty, sphPtz))
                return false;
        }
        offSet += numAtomsPerLevelZeroCube[n];
    }

    return true;
}

void VecOfCubes::SetVecOfCubes(const std::vector<Cube> &vecOfCubes) { m_vecOfCubes = vecOfCubes; }

const std::vector<Cube> &VecOfCubes::GetVecOfCubes() const { return m_vecOfCubes; }

//
// private procs for VecOfCubes
//

std::size_t VecOfCubes::GetTopXIndexAtMyLevel() const { return m_topXIndexPerLevel[m_levelInHeirarchyOfCubes]; }

std::size_t VecOfCubes::GetTopYIndexAtMyLevel() const { return m_topYIndexPerLevel[m_levelInHeirarchyOfCubes]; }

std::size_t VecOfCubes::GetTopZIndexAtMyLevel() const { return m_topZIndexPerLevel[m_levelInHeirarchyOfCubes]; }

std::size_t VecOfCubes::GetNumCubesAtMyLevel() const {
    return GetTopXIndexAtMyLevel() * GetTopYIndexAtMyLevel() * GetTopZIndexAtMyLevel();
}

std::size_t VecOfCubes::GetIndexInVecOfCubesAtLevel(std::size_t indexX, std::size_t indexY, std::size_t indexZ,
                                                    std::size_t level) const {
    return indexZ * m_topXIndexPerLevel[level] * m_topYIndexPerLevel[level] + indexY * m_topXIndexPerLevel[level] +
           indexX;
}

std::size_t VecOfCubes::GetNumLevelsInHeirarchyOfCubes() const { return m_topXIndexPerLevel.size(); }

std::size_t VecOfCubes::GetMyLevelInHeirarchyOfCubes() const { return m_levelInHeirarchyOfCubes; }

std::size_t VecOfCubes::GetNumRadialQuadraturePoints() const { return m_numRadialQuadraturePts; }

bool VecOfCubes::LoadIndicesOfNearestNeighborsOfCube(std::vector<std::size_t> &nearestNeighborList,
                                                     std::size_t indexXOfCube, std::size_t indexYOfCube,
                                                     std::size_t indexZOfCube, std::size_t levelInHeirarchy) {
    nearestNeighborList.clear();

    std::size_t topXIndex = GetTopXIndexAtMyLevel();
    std::size_t topYIndex = GetTopYIndexAtMyLevel();
    std::size_t topZIndex = GetTopZIndexAtMyLevel();
    if (levelInHeirarchy == 0) return false;  // no list for bottom level
    if (topXIndex == 0)
        return false;
    else if (topYIndex == 0)
        return false;
    else if (topZIndex == 0)
        return false;

    std::size_t lowerLimX;
    if (indexXOfCube > 1)
        lowerLimX = indexXOfCube - 1;
    else
        lowerLimX = 0;

    std::size_t upperLimX;
    if (indexXOfCube + 1 < topXIndex)
        upperLimX = indexXOfCube + 1;
    else
        upperLimX = topXIndex - 1;

    std::size_t lowerLimY;
    if (indexYOfCube > 1)
        lowerLimY = indexYOfCube - 1;
    else
        lowerLimY = 0;

    std::size_t upperLimY;
    if (indexYOfCube + 1 < topYIndex)
        upperLimY = indexYOfCube + 1;
    else
        upperLimY = topYIndex - 1;

    std::size_t lowerLimZ;
    if (indexZOfCube > 1)
        lowerLimZ = indexZOfCube - 1;
    else
        lowerLimZ = 0;

    std::size_t upperLimZ;
    if (indexZOfCube + 1 < topZIndex)
        upperLimZ = indexZOfCube + 1;
    else
        upperLimZ = topZIndex - 1;

    for (std::size_t k = lowerLimZ; k <= upperLimZ; ++k) {
        for (std::size_t j = lowerLimY; j <= upperLimY; ++j) {
            for (std::size_t i = lowerLimX; i <= upperLimX; ++i) {
                std::size_t index = GetIndexInVecOfCubesAtLevel(i, j, k, m_levelInHeirarchyOfCubes);
                nearestNeighborList.push_back(index);
            }
        }
    }
    return true;
}

std::size_t VecOfCubes::GetIndexOfParentOfCube(std::size_t levelInHeirarchy, std::size_t indexXOfCube,
                                               std::size_t indexYOfCube, std::size_t indexZOfCube) {
    if (levelInHeirarchy + 1 == GetNumLevelsInHeirarchyOfCubes()) return std::numeric_limits<std::size_t>::max();

    std::size_t indexXOfParent = indexXOfCube / 2;
    std::size_t indexYOfParent = indexYOfCube / 2;
    std::size_t indexZOfParent = indexZOfCube / 2;
    std::size_t indexOfParent =
        GetIndexInVecOfCubesAtLevel(indexXOfParent, indexYOfParent, indexZOfParent, levelInHeirarchy + 1);
    return indexOfParent;
}

void VecOfCubes::SetupVecOfCubes() {
    std::size_t numRadialQuadraturePts = GetNumRadialQuadraturePoints();

    for (std::size_t indexZ = 0; indexZ < GetTopZIndexAtMyLevel(); ++indexZ) {
        for (std::size_t indexY = 0; indexY < GetTopYIndexAtMyLevel(); ++indexY) {
            for (std::size_t indexX = 0; indexX < GetTopXIndexAtMyLevel(); ++indexX) {
                Cube cube(GetNumLevelsInHeirarchyOfCubes(), numRadialQuadraturePts);
                m_vecOfCubes.push_back(cube);
            }
        }
    }
}
