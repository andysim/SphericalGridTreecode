#ifndef CUBE_H
#define CUBE_H

#include <vector>

class Cube {
   public:
    Cube();
    Cube(std::size_t numLevelsInHeirarchyOfCubes, std::size_t numRadialQuadPts);

    void ZeroFillAccumulatedSourcePotentials();

    void SetUpwardsAccumCosPotentials(const std::vector<double>& upwardsAccumCosPotentials);
    const std::vector<double>& GetUpwardsAccumCosPotentials() const;

    void SetUpwardsAccumSinPotentials(const std::vector<double>& upwardsAccumSinPotentials);
    const std::vector<double>& GetUpwardsAccumSinPotentials() const;

    void SetDownwardsAccumCosPotentials(const std::vector<double>& downwardsAccumCosPotentials);
    const std::vector<double>& GetDownwardsAccumCosPotentials() const;

    void SetDownwardsAccumSinPotentials(const std::vector<double>& downwardsAccumSinPotentials);
    const std::vector<double>& GetDownwardsAccumSinPotentials() const;

   private:
    void AllocateAccumulatedSourcePotentials(std::size_t numLevelsInHeirarchyOfCubes, std::size_t numRadialQuadPts);

    std::vector<double> m_upwardsAccumCosPotentials;
    std::vector<double> m_upwardsAccumSinPotentials;
    std::vector<double> m_downwardsAccumCosPotentials;
    std::vector<double> m_downwardsAccumSinPotentials;
};

#endif  // CUBE_H
