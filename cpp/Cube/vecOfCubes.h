#ifndef VEC_OF_CUBES_H
#define VEC_OF_CUBES_H

#include <vector>

#include "../ErfPointCharge/erfPointCharge.h"
#include "../Atom/atom.h"
#include "cube.h"

class VecOfCubes {
   public:
    VecOfCubes();
    VecOfCubes(double widthOfCubes, double minXCoordInVecOfCubes, double minYCoordInVecOfCubes,
               double minZCoordInVecOfCubes, const std::vector<std::size_t> &topXIndexPerLevel,
               const std::vector<std::size_t> &topYIndexPerLevel, const std::vector<std::size_t> &topZIndexPerLevel,
               std::size_t levelInHeirarchyOfCubes, std::size_t numRadialQuadraturePts);

    void ZeroFillAccumulatedSourcePotentials();

    bool AccumulatePotentialsFromErfPointChargesIntoCubes(const std::vector<std::size_t> &numAtomsPerLevelZeroCube,
                                                          const std::vector<ErfPointCharge> &vecOfErfPointCharges);

    bool AddUpwardAccumPotentialsToThoseOfParentCubes(std::vector<Cube> &vecOfCubesAtParentLevel);

    bool SumUpwardAccumPotsAtTopLevelAndPutIntoDownwardAccumPots();

    bool AdjustDownwardAccumPotentialsByNearestNeighbors();

    bool GetDownwardAccumPotentialsFromThoseOfParentCubes(std::vector<Cube> &vecOfCubesAtParentLevel);

    bool AddDownwardAccumCubePotentialsToErfPointCharges(std::vector<ErfPointCharge> &vecOfErfPointCharges,
                                                         const std::vector<std::size_t> &numAtomsPerLevelZeroCube,
                                                         double sphPtx, double sphPty, double sphPtz);

    void SetVecOfCubes(const std::vector<Cube> &vecOfCubes);
    const std::vector<Cube> &GetVecOfCubes() const;

   private:
    std::size_t GetTopXIndexAtMyLevel() const;
    std::size_t GetTopYIndexAtMyLevel() const;
    std::size_t GetTopZIndexAtMyLevel() const;

    std::size_t GetNumCubesAtMyLevel() const;

    std::size_t GetIndexInVecOfCubesAtLevel(std::size_t indexX, std::size_t indexY, std::size_t indexZ,
                                            std::size_t level) const;

    std::size_t GetNumLevelsInHeirarchyOfCubes() const;

    std::size_t GetMyLevelInHeirarchyOfCubes() const;

    std::size_t GetNumRadialQuadraturePoints() const;

    bool LoadIndicesOfNearestNeighborsOfCube(std::vector<std::size_t> &nearestNeighborList, std::size_t indexXOfCube,
                                             std::size_t indexYOfCube, std::size_t indexZOfCube,
                                             std::size_t levelInHeirarchy);

    std::size_t GetIndexOfParentOfCube(std::size_t levelInHeirarchy, std::size_t indexXOfCube, std::size_t indexYOfCube,
                                       std::size_t indexZOfCube);

    void SetupVecOfCubes();

    double m_widthOfCubes = 0.0;
    double m_minXCoordInVecOfCubes = std::numeric_limits<double>::max();
    double m_minYCoordInVecOfCubes = std::numeric_limits<double>::max();
    double m_minZCoordInVecOfCubes = std::numeric_limits<double>::max();
    std::vector<std::size_t> m_topXIndexPerLevel;
    std::vector<std::size_t> m_topYIndexPerLevel;
    std::vector<std::size_t> m_topZIndexPerLevel;
    std::size_t m_levelInHeirarchyOfCubes = std::numeric_limits<double>::max();
    std::size_t m_numRadialQuadraturePts = 0;
    std::vector<Cube> m_vecOfCubes;
};

#endif  // VEC_OF_CUBES_H
