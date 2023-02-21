#ifndef ATOM
#define ATOM

#include <vector>

class Atom {
   public:
    Atom();
    Atom(const std::string &name, double charge, double xCoord, double yCoord, double zCoord,
         std::size_t indexInListOfAtoms, const std::vector<std::size_t> &exclusionList);

    const std::string &GetName() const;
    double GetCharge() const;
    double GetXCoord() const;
    double GetYCoord() const;
    double GetZCoord() const;
    std::size_t GetIndexInListOfAtoms() const;

    void SetIndexInListOfCubes(std::size_t indexInListOfCubes);
    std::size_t GetIndexInListOfCubes() const;

    const std::vector<std::size_t> &GetExclusionList() const;

    void LoadExactCoulombPotential_and_Grad(std::vector<Atom> &atomsInVec);

    void GetExactCoulombPotential_and_Grad(double &exactCoulombPotential,
                                           std::vector<double> &gradOfExactCoulombPotential) const;

    void LoadCutoffErfc_and_ExclusionListErf_Potential_and_Grad(std::vector<Atom> &atoms, double cutoff,
                                                                double ewaldCoeff);

    void SetApproxErfPotential_and_Grad(double approxErfPotential, const std::vector<double> &gradOfApproxErfPotential);

    void SetApproxCoulombPotential_and_Grad();

    void GetApproxCoulombPotential_and_Grad(double &approxCoulombPotential,
                                            std::vector<double> &gradOfApproxCoulombPotential) const;

    void GetErrorInApproxCoulombPotential_and_Grad(double &absErrorInApproxCoulombPotential,
                                                   double &sumSqErrorInGradOfApproxCoulombPotential);

    std::string PrintPots() const;
    void Clear();

   private:
    void IncrementCoulombPotential_and_Grad(const Atom &atom);
    void IncrementCutoffErfc_and_ExclusionListErf_Potential_and_Grad(const Atom &atom, double cutoff,
                                                                     double ewaldCoeff);

    std::string m_name = "";
    double m_charge = 0.0;
    double m_xCoord = std::numeric_limits<double>::max();
    double m_yCoord = std::numeric_limits<double>::max();
    double m_zCoord = std::numeric_limits<double>::max();
    std::size_t m_indexInListOfAtoms = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> m_exclusionList;  // eg atoms on same water

    std::size_t m_indexInListOfCubes = std::numeric_limits<std::size_t>::max();

    double m_exactCoulombPot = 0.0;
    double m_approxCoulombPot = 0.0;

    std::vector<double> m_gradOfExactCoulombPot;
    std::vector<double> m_gradOfApproxCoulombPot;

    double m_cutoffErfcPot = 0.0;
    double m_exclusionListErfPot = 0.0;
    double m_approxErfPotential = 0.0;

    std::vector<double> m_gradOfCutoffErfcPot;
    std::vector<double> m_gradOfExclusionListErfPot;
    std::vector<double> m_gradOfApproxErfPotential;
};

#endif  // ATOM
