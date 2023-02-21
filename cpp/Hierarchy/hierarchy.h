#ifndef HIERARCHY_H
#define HIERARCHY_H

#include <vector>

#include "../ErfPointCharge/erfPointCharge.h"
#include "../Atom/atom.h"
#include "../Cube/vecOfCubes.h"

class Hierarchy {
   public:
    Hierarchy(double widthOfCubes, double minXCoordInMol, double minYCoordInMol, double minZCoordInMol,
              double ewaldCoeff, std::size_t numSphericalQuadraturePoints, std::vector<double> &gaussHermiteRoots,
              std::vector<double> &gaussHermiteWeights, double maxXCoordInMol, double maxYCoordInMol,
              double maxZCoordInMol, const std::vector<Atom> &vecOfAtoms);

    bool IsValid() const;

    bool AccumulatePotentialsForNewSphericalQuadraturePoint(double sphPtx, double sphPty, double sphPtz);

    void LoadApproxErfPotential_and_Grad();

    bool GetVecsOfApproxErfPot_and_Grads(std::vector<double> &approxErfPots,
                                         std::vector<double> &approxGradErfPot_wrt_X,
                                         std::vector<double> &approxGradErfPot_wrt_Y,
                                         std::vector<double> &approxGradErfPot_wrt_Z, std::size_t sizeOfListOfAtoms);

   private:
    void IntegrateAccumPotentials();

    void SetNewSphericalQuadraturePoint(double sphPtx, double sphPty, double sphPtz);

    void LoadCosineSineBasisFunctions();

    bool AccumulatePotentialsUpwardFromChildToParent();

    bool AccumulatePotentialsDownwardFromParentToChild();

    std::size_t GetNumLevelsInHierarchyOfCubes() const;

    std::size_t SetupTopIndicesAtLevel(std::size_t &indexX, std::size_t &indexY, std::size_t &indexZ,
                                       double xWidthOfMol, double yWidthOfMol, double zWidthOfMol,
                                       double cubeWidthAtThisLevel);

    void SetupTopIndices(double widthOfCubes, double minXCoordInMol, double minYCoordInMol, double minZCoordInMol,
                         double maxXCoordInMol, double maxYCoordInMol, double maxZCoordInMol);

    void SetupVecOfVecOfCubes();

    std::size_t GetIndexInZeroLevelVecOfCubes(std::size_t indexX, std::size_t indexY, std::size_t indexZ) const;

    bool SetIndexInZeroLevelVecOfCubesToAtoms(std::vector<Atom> &tmpVec, const std::vector<Atom> &oldVec) const;

    bool SortVecOfAtomsByCubeIndex(std::vector<Atom> &newVec, const std::vector<Atom> &oldVec);

    void AtomsToErfPointCharges(std::vector<ErfPointCharge> &vecOfErfPointCharges,
                                const std::vector<Atom> &vecOfAtoms) const;

    bool SetupVecOfErfPointCharges(const std::vector<Atom> &vecOfAtoms);

    void AllocatePointChargeArrays(std::size_t numLevelsInHierarchyOfCubes, std::size_t numRadialQuadPts);

    double m_widthOfCubes = 0.0;
    double m_minXCoordInMol = std::numeric_limits<double>::max();
    double m_minYCoordInMol = std::numeric_limits<double>::max();
    double m_minZCoordInMol = std::numeric_limits<double>::max();
    double m_ewaldCoeff = 0.0;
    std::size_t m_numSphericalQuadraturePoints = 0;
    std::vector<double> m_gaussHermiteRoots;
    std::vector<double> m_gaussHermiteWeights;

    bool m_valid = false;

    std::vector<std::size_t> m_topIndexX;
    std::vector<std::size_t> m_topIndexY;
    std::vector<std::size_t> m_topIndexZ;

    double m_spherePointX = 0.0;
    double m_spherePointY = 0.0;
    double m_spherePointZ = 0.0;

    std::vector<VecOfCubes> m_vecOfVecOfCubes;

    std::vector<std::size_t> m_numAtomsPerLevelZeroCube;
    std::vector<ErfPointCharge> m_vecOfErfPointCharges;
};

#endif  // HIERARCHY_H
