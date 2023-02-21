#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include "atom.h"

// procs for Atom

Atom::Atom() {}

Atom::Atom(const std::string &name, double charge, double xCoord, double yCoord, double zCoord,
           std::size_t indexInListOfAtoms, const std::vector<std::size_t> &exclusionList)
    : m_name(name),
      m_charge(charge),
      m_xCoord(xCoord),
      m_yCoord(yCoord),
      m_zCoord(zCoord),
      m_indexInListOfAtoms(indexInListOfAtoms),
      m_exclusionList(exclusionList) {}

const std::string &Atom::GetName() const { return m_name; }

double Atom::GetCharge() const { return m_charge; }

double Atom::GetXCoord() const { return m_xCoord; }

double Atom::GetYCoord() const { return m_yCoord; }

double Atom::GetZCoord() const { return m_zCoord; }

std::size_t Atom::GetIndexInListOfAtoms() const { return m_indexInListOfAtoms; }

void Atom::SetIndexInListOfCubes(std::size_t indexInListOfCubes) { m_indexInListOfCubes = indexInListOfCubes; }

std::size_t Atom::GetIndexInListOfCubes() const { return m_indexInListOfCubes; }

const std::vector<std::size_t> &Atom::GetExclusionList() const { return m_exclusionList; }

void Atom::LoadExactCoulombPotential_and_Grad(std::vector<Atom> &atoms) {
    m_exactCoulombPot = 0.0;
    m_gradOfExactCoulombPot.resize(3, 0.0);
    for (auto atom : atoms) IncrementCoulombPotential_and_Grad(atom);
}

void Atom::GetExactCoulombPotential_and_Grad(double &exactCoulombPotential,
                                             std::vector<double> &gradOfExactCoulombPotential) const {
    exactCoulombPotential = m_exactCoulombPot;
    gradOfExactCoulombPotential = m_gradOfExactCoulombPot;
}

void Atom::LoadCutoffErfc_and_ExclusionListErf_Potential_and_Grad(std::vector<Atom> &atoms, double cutoff,
                                                                  double ewaldCoeff) {
    m_cutoffErfcPot = 0.0;
    m_exclusionListErfPot = 0.0;

    m_gradOfCutoffErfcPot.resize(3, 0.0);
    m_gradOfExclusionListErfPot.resize(3, 0.0);
    for (auto atom : atoms) IncrementCutoffErfc_and_ExclusionListErf_Potential_and_Grad(atom, cutoff, ewaldCoeff);
}

void Atom::SetApproxErfPotential_and_Grad(double approxErfPotential,
                                          const std::vector<double> &gradOfApproxErfPotential) {
    m_approxErfPotential = approxErfPotential;
    m_gradOfApproxErfPotential = gradOfApproxErfPotential;
}

void Atom::SetApproxCoulombPotential_and_Grad() {
    m_approxCoulombPot = m_cutoffErfcPot + m_approxErfPotential - m_exclusionListErfPot;

    m_gradOfApproxCoulombPot.resize(3, 0.0);
    for (std::size_t j = 0; j < 3; ++j) {
        m_gradOfApproxCoulombPot[j] =
            m_gradOfCutoffErfcPot[j] + m_gradOfApproxErfPotential[j] - m_gradOfExclusionListErfPot[j];
    }
}

void Atom::GetApproxCoulombPotential_and_Grad(double &approxCoulombPotential,
                                              std::vector<double> &gradOfApproxCoulombPotential) const {
    approxCoulombPotential = m_approxCoulombPot;
    gradOfApproxCoulombPotential = m_gradOfApproxCoulombPot;
}

void Atom::GetErrorInApproxCoulombPotential_and_Grad(double &absErrorInApproxCoulombPotential,
                                                     double &sumSqErrorInGradOfApproxCoulombPotential) {
    double exactCoulombPotential;
    std::vector<double> gradOfExactCoulombPotential;
    GetExactCoulombPotential_and_Grad(exactCoulombPotential, gradOfExactCoulombPotential);

    double approxCoulombPotential;
    std::vector<double> gradOfApproxCoulombPotential;
    GetApproxCoulombPotential_and_Grad(approxCoulombPotential, gradOfApproxCoulombPotential);

    absErrorInApproxCoulombPotential = abs(exactCoulombPotential - approxCoulombPotential);
    sumSqErrorInGradOfApproxCoulombPotential = 0.0;
    for (std::size_t j = 0; j < 3; ++j) {
        double error = gradOfExactCoulombPotential[j] - gradOfApproxCoulombPotential[j];
        sumSqErrorInGradOfApproxCoulombPotential += error * error;
    }
}

std::string Atom::PrintPots() const {
    double exactCoulombPotential = 0.0;
    std::vector<double> gradOfExactCoulomb;
    GetExactCoulombPotential_and_Grad(exactCoulombPotential, gradOfExactCoulomb);
    double approxCoulombPotential = 0.0;
    std::vector<double> gradOfApproxCoulombPotential;
    GetApproxCoulombPotential_and_Grad(approxCoulombPotential, gradOfApproxCoulombPotential);

    std::stringstream ss;
    ss << "exact  Coulomb pot: " << exactCoulombPotential << std::endl;
    ss << "gradient of Exact : ";
    for (std::size_t i = 0; i < 3; ++i) ss << gradOfExactCoulomb[i] << ", ";
    ss << std::endl;

    ss << "approx Coulomb pot: " << approxCoulombPotential << std::endl;
    ss << "gradient of Approx: ";
    for (std::size_t i = 0; i < 3; ++i) ss << gradOfApproxCoulombPotential[i] << ", ";
    ss << std::endl;

    return ss.str();
}

void Atom::Clear() {
    m_name = "";
    m_charge = 0.0;
    m_xCoord = std::numeric_limits<double>::max();
    m_yCoord = std::numeric_limits<double>::max();
    m_zCoord = std::numeric_limits<double>::max();
    m_indexInListOfAtoms = std::numeric_limits<std::size_t>::max();
    m_exclusionList.clear();

    m_exactCoulombPot = 0.0;
    m_approxCoulombPot = 0.0;

    m_gradOfExactCoulombPot.clear();
    m_gradOfApproxCoulombPot.clear();

    m_cutoffErfcPot = 0.0;
    m_exclusionListErfPot = 0.0;
    m_approxErfPotential = 0.0;

    m_gradOfCutoffErfcPot.clear();
    m_gradOfExclusionListErfPot.clear();
    m_gradOfApproxErfPotential.clear();
}

// private procs for Atom

void Atom::IncrementCoulombPotential_and_Grad(const Atom &atom) {
    std::size_t index = atom.GetIndexInListOfAtoms();
    bool excluded = false;
    // filter all but exclusion ListErf
    for (auto ind : m_exclusionList) {
        if (ind == index) {
            excluded = true;
            break;
        }
    }
    if (!excluded) {
        double charge = atom.GetCharge();
        // get the distance
        double xCoord = atom.GetXCoord();
        double yCoord = atom.GetYCoord();
        double zCoord = atom.GetZCoord();
        double r12 = sqrt((m_xCoord - xCoord) * (m_xCoord - xCoord) + (m_yCoord - yCoord) * (m_yCoord - yCoord) +
                          (m_zCoord - zCoord) * (m_zCoord - zCoord));
        m_exactCoulombPot += charge / r12;
        double dPot_dr = -charge / (r12 * r12);
        double dr_dMyX = (m_xCoord - xCoord) / r12;
        double dr_dMyY = (m_yCoord - yCoord) / r12;
        double dr_dMyZ = (m_zCoord - zCoord) / r12;
        m_gradOfExactCoulombPot[0] += dPot_dr * dr_dMyX;
        m_gradOfExactCoulombPot[1] += dPot_dr * dr_dMyY;
        m_gradOfExactCoulombPot[2] += dPot_dr * dr_dMyZ;
    }
}

void Atom::IncrementCutoffErfc_and_ExclusionListErf_Potential_and_Grad(const Atom &atom, double cutoff,
                                                                       double ewaldCoeff) {
    std::size_t index = atom.GetIndexInListOfAtoms();
    double charge = atom.GetCharge();
    // get the distance
    double xCoord = atom.GetXCoord();
    double yCoord = atom.GetYCoord();
    double zCoord = atom.GetZCoord();
    double r12 = sqrt((xCoord - m_xCoord) * (xCoord - m_xCoord) + (yCoord - m_yCoord) * (yCoord - m_yCoord) +
                      (zCoord - m_zCoord) * (zCoord - m_zCoord));
    bool excluded = false;
    // filter all but exclusion ListErf
    for (auto ind : m_exclusionList) {
        if (ind == index) {
            excluded = true;
            break;
        }
    }
    if (excluded) {
        if (index == m_indexInListOfAtoms) {
            // r12 = 0.0 -- take limit
            m_exclusionListErfPot += 2.0 * charge * ewaldCoeff / sqrt(M_PI);
            // force increment is zero
        } else {
            m_exclusionListErfPot += charge * erf(ewaldCoeff * r12) / r12;
            double dPot_dr =
                charge * ((2.0 * ewaldCoeff / sqrt(M_PI)) * exp(-ewaldCoeff * ewaldCoeff * r12 * r12) / r12 -
                          erf(ewaldCoeff * r12) / (r12 * r12));
            double dr_dMyX = (m_xCoord - xCoord) / r12;
            double dr_dMyY = (m_yCoord - yCoord) / r12;
            double dr_dMyZ = (m_zCoord - zCoord) / r12;
            m_gradOfExclusionListErfPot[0] += dPot_dr * dr_dMyX;
            m_gradOfExclusionListErfPot[1] += dPot_dr * dr_dMyY;
            m_gradOfExclusionListErfPot[2] += dPot_dr * dr_dMyZ;
        }
    } else {
        if (r12 <= cutoff) m_cutoffErfcPot += charge * erfc(ewaldCoeff * r12) / r12;
        double dCoulPot_dr = -charge / (r12 * r12);
        double dErfPot_dr =
            charge * ((2.0 * ewaldCoeff / sqrt(M_PI)) * exp(-ewaldCoeff * ewaldCoeff * r12 * r12) / r12 -
                      erf(ewaldCoeff * r12) / (r12 * r12));
        double dPot_dr = dCoulPot_dr - dErfPot_dr;
        double dr_dMyX = (m_xCoord - xCoord) / r12;
        double dr_dMyY = (m_yCoord - yCoord) / r12;
        double dr_dMyZ = (m_zCoord - zCoord) / r12;
        m_gradOfCutoffErfcPot[0] += dPot_dr * dr_dMyX;
        m_gradOfCutoffErfcPot[1] += dPot_dr * dr_dMyY;
        m_gradOfCutoffErfcPot[2] += dPot_dr * dr_dMyZ;
    }
}
