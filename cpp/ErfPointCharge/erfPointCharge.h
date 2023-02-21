#ifndef ERF_POINT_CHARGE_H
#define ERF_POINT_CHARGE_H

#include <vector>
#include <string>
#include <sstream>

class ErfPointCharge {
   public:
    ErfPointCharge();
    ErfPointCharge(double charge, double xCoord, double yCoord, double zCoord, std::size_t indexInListOfAtoms);
    ErfPointCharge(const ErfPointCharge &rhs);
    ErfPointCharge &operator=(const ErfPointCharge &rhs);
    ~ErfPointCharge() {}

    double GetCharge() const;
    double GetXCoord() const;
    double GetYCoord() const;
    double GetZCoord() const;

    std::size_t GetIndexInListOfAtoms() const;

    void LoadApproxErfPotential_and_Grad();

    void GetApproxErfPotential_and_Grad(double &approxErfPotential,
                                        std::vector<double> &gradOfApproxErfPotential) const;

    void AllocateArrays(std::size_t numRows, std::size_t numInRow);
    void InitializePotentials();

    void LoadCosineSineBasisFunctions(double ewaldCoeff, std::size_t numSpherePts, std::size_t numRows,
                                      const std::vector<double> &gaussHermiteRoots, double sphPtx, double sphPty,
                                      double sphPtz);

    bool AddMySourcePotentialsToAccumulatedSourcePotentials(std::vector<double> &accumulatedSourceCosPotentials,
                                                            std::vector<double> &accumulatedSourceSinPotentials) const;

    bool AddAccumCosSinPotentialsToMyAccumCosSinPotentials(const std::vector<double> &accumulatedSourceCosPotentials,
                                                           const std::vector<double> &accumulatedSourceSinPotentials,
                                                           double sphPtx, double sphPty, double sphPtz);

    void IntegrateAccumPotentials(double ewaldCoeff, std::size_t numRows,
                                  const std::vector<double> &gaussHermiteWeights);

    std::string PrintApproxErfPotAndGrad() const;
    void Clear();

   private:
    double m_charge = 0.0;
    double m_xCoord = std::numeric_limits<double>::max();
    double m_yCoord = std::numeric_limits<double>::max();
    double m_zCoord = std::numeric_limits<double>::max();
    std::size_t m_indexInListOfAtoms = std::numeric_limits<std::size_t>::max();
    double m_approxErfPotential = 0.0;
    std::vector<double> m_gradOfApproxErfPotential;
    std::vector<double> m_integratedCosPotential;
    std::vector<double> m_integrated_dCosPotentialDx;
    std::vector<double> m_integrated_dCosPotentialDy;
    std::vector<double> m_integrated_dCosPotentialDz;

    std::vector<double> m_accumCosPotential;
    std::vector<double> m_dAccumCosPotentialDx;
    std::vector<double> m_dAccumCosPotentialDy;
    std::vector<double> m_dAccumCosPotentialDz;

    std::vector<double> m_cosBasisFunctions;
    std::vector<double> m_sinBasisFunctions;
    std::vector<double> m_dArg_dDot;
};

#endif  // ERF_POINT_CHARGE_H
