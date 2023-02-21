#include <vector>
#include <cmath>
#include <iostream>

#include "erfPointCharge.h"

// procs for ErfPointCharge

ErfPointCharge::ErfPointCharge() {}

ErfPointCharge::ErfPointCharge(double charge, double xCoord, double yCoord, double zCoord,
                               std::size_t indexInListOfAtoms)
    : m_charge(charge),
      m_xCoord(xCoord),
      m_yCoord(yCoord),
      m_zCoord(zCoord),
      m_indexInListOfAtoms(indexInListOfAtoms) {}

ErfPointCharge::ErfPointCharge(const ErfPointCharge &rhs)
    : m_charge(rhs.m_charge),
      m_xCoord(rhs.m_xCoord),
      m_yCoord(rhs.m_yCoord),
      m_zCoord(rhs.m_zCoord),
      m_indexInListOfAtoms(rhs.m_indexInListOfAtoms),

      m_approxErfPotential(rhs.m_approxErfPotential),
      m_gradOfApproxErfPotential(rhs.m_gradOfApproxErfPotential),

      m_integratedCosPotential(rhs.m_integratedCosPotential),
      m_integrated_dCosPotentialDx(rhs.m_integrated_dCosPotentialDx),
      m_integrated_dCosPotentialDy(rhs.m_integrated_dCosPotentialDy),
      m_integrated_dCosPotentialDz(rhs.m_integrated_dCosPotentialDz),

      m_accumCosPotential(rhs.m_accumCosPotential),
      m_dAccumCosPotentialDx(rhs.m_dAccumCosPotentialDx),
      m_dAccumCosPotentialDy(rhs.m_dAccumCosPotentialDy),
      m_dAccumCosPotentialDz(rhs.m_dAccumCosPotentialDz),

      m_cosBasisFunctions(rhs.m_cosBasisFunctions),
      m_sinBasisFunctions(rhs.m_sinBasisFunctions),
      m_dArg_dDot(rhs.m_dArg_dDot) {}

ErfPointCharge &ErfPointCharge::operator=(const ErfPointCharge &rhs) {
    if (this != &rhs) {
        m_charge = rhs.m_charge;
        m_xCoord = rhs.m_xCoord;
        m_yCoord = rhs.m_yCoord;
        m_zCoord = rhs.m_zCoord;
        m_indexInListOfAtoms = rhs.m_indexInListOfAtoms;

        m_approxErfPotential = rhs.m_approxErfPotential;
        m_gradOfApproxErfPotential = rhs.m_gradOfApproxErfPotential;

        m_integratedCosPotential = rhs.m_integratedCosPotential;
        m_integrated_dCosPotentialDx = rhs.m_integrated_dCosPotentialDx;
        m_integrated_dCosPotentialDy = rhs.m_integrated_dCosPotentialDz;
        m_integrated_dCosPotentialDz = rhs.m_integrated_dCosPotentialDz;

        m_accumCosPotential = rhs.m_accumCosPotential;
        m_dAccumCosPotentialDx = rhs.m_dAccumCosPotentialDx;
        m_dAccumCosPotentialDy = rhs.m_dAccumCosPotentialDy;
        m_dAccumCosPotentialDz = rhs.m_dAccumCosPotentialDz;

        m_cosBasisFunctions = rhs.m_cosBasisFunctions;
        m_sinBasisFunctions = rhs.m_sinBasisFunctions;
        m_dArg_dDot = rhs.m_dArg_dDot;
    }
    return *this;
}

double ErfPointCharge::GetCharge() const { return m_charge; }

double ErfPointCharge::GetXCoord() const { return m_xCoord; }

double ErfPointCharge::GetYCoord() const { return m_yCoord; }

double ErfPointCharge::GetZCoord() const { return m_zCoord; }

std::size_t ErfPointCharge::GetIndexInListOfAtoms() const { return m_indexInListOfAtoms; }

void ErfPointCharge::LoadApproxErfPotential_and_Grad() {
    m_approxErfPotential = 0.0;
    m_gradOfApproxErfPotential.resize(3, 0.0);
    std::size_t numRows = m_integratedCosPotential.size();

    for (std::size_t k = 0; k < numRows; ++k) {
        m_approxErfPotential += m_integratedCosPotential[k];
        m_gradOfApproxErfPotential[0] += m_integrated_dCosPotentialDx[k];
        m_gradOfApproxErfPotential[1] += m_integrated_dCosPotentialDy[k];
        m_gradOfApproxErfPotential[2] += m_integrated_dCosPotentialDz[k];
    }
}

void ErfPointCharge::GetApproxErfPotential_and_Grad(double &approxErfPotential,
                                                    std::vector<double> &gradOfApproxErfPotential) const {
    approxErfPotential = m_approxErfPotential;
    gradOfApproxErfPotential = m_gradOfApproxErfPotential;
}

void ErfPointCharge::AllocateArrays(std::size_t numRows, std::size_t numInRow) {
    std::size_t numElements = numRows * numInRow;

    m_integratedCosPotential.resize(numRows, 0.0);
    m_integrated_dCosPotentialDx.resize(numRows, 0.0);
    m_integrated_dCosPotentialDy.resize(numRows, 0.0);
    m_integrated_dCosPotentialDz.resize(numRows, 0.0);

    m_accumCosPotential.resize(numElements);
    m_dAccumCosPotentialDx.resize(numElements);
    m_dAccumCosPotentialDy.resize(numElements);
    m_dAccumCosPotentialDz.resize(numElements);

    m_cosBasisFunctions.resize(numElements);
    m_sinBasisFunctions.resize(numElements);
    m_dArg_dDot.resize(numElements);
}

void ErfPointCharge::InitializePotentials() {
    std::fill(m_integratedCosPotential.begin(), m_integratedCosPotential.end(), 0.0);
    std::fill(m_integrated_dCosPotentialDx.begin(), m_integrated_dCosPotentialDx.end(), 0.0);
    std::fill(m_integrated_dCosPotentialDy.begin(), m_integrated_dCosPotentialDy.end(), 0.0);
    std::fill(m_integrated_dCosPotentialDz.begin(), m_integrated_dCosPotentialDz.end(), 0.0);

    std::fill(m_accumCosPotential.begin(), m_accumCosPotential.end(), 0.0);
    std::fill(m_dAccumCosPotentialDx.begin(), m_dAccumCosPotentialDx.end(), 0.0);
    std::fill(m_dAccumCosPotentialDy.begin(), m_dAccumCosPotentialDy.end(), 0.0);
    std::fill(m_dAccumCosPotentialDz.begin(), m_dAccumCosPotentialDz.end(), 0.0);
}

void ErfPointCharge::LoadCosineSineBasisFunctions(double ewaldCoeff, std::size_t numSpherePts, std::size_t numRows,
                                                  const std::vector<double> &gaussHermiteRoots, double sphPtx,
                                                  double sphPty, double sphPtz) {
    // lowest indices correspond to highest level, smallest ewald coeff
    //  each level up we halve the ewald coeff
    //  in this data structure, for each row up we get cos, sin of double ewald coeff

    std::size_t numInRow = gaussHermiteRoots.size();

    double currEwaldCoeff = ewaldCoeff;
    for (std::size_t n = 1; n < numRows; ++n) currEwaldCoeff *= 0.5;

    double dotProd = sphPtx * GetXCoord() + sphPty * GetYCoord() + sphPtz * GetZCoord();

    double argPrefac = 2.0 * currEwaldCoeff;
    double multiplier = 1.0 / static_cast<double>(numSpherePts);
    double weightPrefac = sqrt(multiplier);
    for (std::size_t n = 0; n < numInRow; ++n) {
        m_dArg_dDot[n] = 2.0 * currEwaldCoeff * gaussHermiteRoots[n];
        double arg = dotProd * m_dArg_dDot[n];
        m_cosBasisFunctions[n] = cos(arg);
        m_sinBasisFunctions[n] = sin(arg);
    }

    if (numRows > 1) {
        for (std::size_t k = 1; k < numRows; ++k) {
            std::size_t startFrom = (k - 1) * numInRow;
            std::size_t startTo = k * numInRow;
            for (std::size_t n = 0; n < numInRow; ++n) {
                // arg doubles each row
                double dArg_dDotFrom = m_dArg_dDot[startFrom + n];
                double dArg_dDotTo = 2.0 * dArg_dDotFrom;
                m_dArg_dDot[startTo + n] = dArg_dDotTo;

                double cosFrom = m_cosBasisFunctions[startFrom + n];
                double sinFrom = m_sinBasisFunctions[startFrom + n];
                double cosTo = cosFrom * cosFrom - sinFrom * sinFrom;
                double sinTo = 2.0 * sinFrom * cosFrom;
                m_cosBasisFunctions[startTo + n] = cosTo;
                m_sinBasisFunctions[startTo + n] = sinTo;
            }
        }
    }

    std::size_t size = numRows * numInRow;
    for (std::size_t j = 0; j < size; ++j) {
        m_cosBasisFunctions[j] *= weightPrefac;
        m_sinBasisFunctions[j] *= weightPrefac;
    }
}

bool ErfPointCharge::AddMySourcePotentialsToAccumulatedSourcePotentials(
    std::vector<double> &accumulatedSourceCosPotentials, std::vector<double> &accumulatedSourceSinPotentials) const {
    if (accumulatedSourceCosPotentials.size() != m_cosBasisFunctions.size()) return false;
    if (accumulatedSourceSinPotentials.size() != m_sinBasisFunctions.size()) return false;

    double charge = GetCharge();
    std::size_t size = m_cosBasisFunctions.size();
    for (std::size_t n = 0; n < size; ++n) accumulatedSourceCosPotentials[n] += charge * m_cosBasisFunctions[n];

    for (std::size_t n = 0; n < size; ++n) accumulatedSourceSinPotentials[n] += charge * m_sinBasisFunctions[n];

    return true;
}

bool ErfPointCharge::AddAccumCosSinPotentialsToMyAccumCosSinPotentials(
    const std::vector<double> &accumulatedSourceCosPotentials,
    const std::vector<double> &accumulatedSourceSinPotentials, double sphPtx, double sphPty, double sphPtz) {
    if (accumulatedSourceCosPotentials.size() != m_cosBasisFunctions.size()) return false;
    if (accumulatedSourceSinPotentials.size() != m_sinBasisFunctions.size()) return false;

    std::size_t size = m_cosBasisFunctions.size();
    for (std::size_t n = 0; n < size; ++n) {
        // working from minus of angle in imaginary exp form
        // term is product of exp and exp complement or exp of diff of args
        // note we are accumulating; adding over the different sphere points
        double accumCosPot = accumulatedSourceCosPotentials[n];
        double accumSinPot = accumulatedSourceSinPotentials[n];

        double cosPot = m_cosBasisFunctions[n];
        double sinPot = m_sinBasisFunctions[n];
        double cosDiff = cosPot * accumCosPot + sinPot * accumSinPot;
        double sinDiff = sinPot * accumCosPot - cosPot * accumSinPot;
        m_accumCosPotential[n] += cosDiff;

        double dCos_dArg = -sinDiff;
        double dArg_dDot = m_dArg_dDot[n];
        double dDot_dX = sphPtx;
        double dDot_dY = sphPty;
        double dDot_dZ = sphPtz;
        m_dAccumCosPotentialDx[n] += dCos_dArg * dArg_dDot * dDot_dX;
        m_dAccumCosPotentialDy[n] += dCos_dArg * dArg_dDot * dDot_dY;
        m_dAccumCosPotentialDz[n] += dCos_dArg * dArg_dDot * dDot_dZ;
    }

    return true;
}

void ErfPointCharge::IntegrateAccumPotentials(double ewaldCoeff, std::size_t numRows,
                                              const std::vector<double> &gaussHermiteWeights) {
    // this is done after sum over sphere points
    double currEwaldCoeff = ewaldCoeff;
    for (std::size_t n = 1; n < numRows; ++n) currEwaldCoeff *= 0.5;

    std::size_t numPerRow = gaussHermiteWeights.size();

    double charge = GetCharge();
    for (std::size_t k = 0; k < numRows; ++k) {
        double integralCos = 0.0;
        double integral_dCosDx = 0.0;
        double integral_dCosDy = 0.0;
        double integral_dCosDz = 0.0;
        std::size_t startFrom = k * numPerRow;
        for (std::size_t n = 0; n < numPerRow; ++n) {
            integralCos += gaussHermiteWeights[n] * m_accumCosPotential[startFrom + n];
            integral_dCosDx += gaussHermiteWeights[n] * m_dAccumCosPotentialDx[startFrom + n];
            integral_dCosDy += gaussHermiteWeights[n] * m_dAccumCosPotentialDy[startFrom + n];
            integral_dCosDz += gaussHermiteWeights[n] * m_dAccumCosPotentialDz[startFrom + n];
        }
        integralCos *= (4.0 * currEwaldCoeff / M_PI);
        integral_dCosDx *= (4.0 * currEwaldCoeff / M_PI);
        integral_dCosDy *= (4.0 * currEwaldCoeff / M_PI);
        integral_dCosDz *= (4.0 * currEwaldCoeff / M_PI);
        m_integratedCosPotential[k] = integralCos;
        m_integrated_dCosPotentialDx[k] = integral_dCosDx;
        m_integrated_dCosPotentialDy[k] = integral_dCosDy;
        m_integrated_dCosPotentialDz[k] = integral_dCosDz;
        currEwaldCoeff *= 2.0;
    }
}

std::string ErfPointCharge::PrintApproxErfPotAndGrad() const {
    double approxErfPot = 0.0;
    std::vector<double> gradOfApproxErfPot;
    GetApproxErfPotential_and_Grad(approxErfPot, gradOfApproxErfPot);

    std::stringstream ss;
    ss << "approx     Erf pot: " << approxErfPot << std::endl;
    ss << "gradient of Approx: ";
    for (std::size_t i = 0; i < 3; ++i) ss << gradOfApproxErfPot[i] << ", ";
    ss << std::endl;

    return ss.str();
}

void ErfPointCharge::Clear() {
    m_charge = 0.0;
    m_xCoord = std::numeric_limits<double>::max();
    m_yCoord = std::numeric_limits<double>::max();
    m_zCoord = std::numeric_limits<double>::max();
    m_indexInListOfAtoms = std::numeric_limits<std::size_t>::max();

    m_approxErfPotential = 0.0;
    m_gradOfApproxErfPotential.clear();

    m_integratedCosPotential.clear();
    m_integrated_dCosPotentialDx.clear();
    m_integrated_dCosPotentialDy.clear();
    m_integrated_dCosPotentialDz.clear();

    m_accumCosPotential.clear();
    m_dAccumCosPotentialDx.clear();
    m_dAccumCosPotentialDy.clear();
    m_dAccumCosPotentialDz.clear();

    m_cosBasisFunctions.clear();
    m_sinBasisFunctions.clear();
    m_dArg_dDot.clear();
}
