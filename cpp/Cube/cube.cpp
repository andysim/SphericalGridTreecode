#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include "cube.h"

// procs for class Cube
//
Cube::Cube() {}

Cube::Cube(std::size_t numLevelsInHeirarchyOfCubes, std::size_t numRadialQuadPts) {
    AllocateAccumulatedSourcePotentials(numLevelsInHeirarchyOfCubes, numRadialQuadPts);
}

void Cube::ZeroFillAccumulatedSourcePotentials() {
    // re-init the accumulators
    std::fill(m_upwardsAccumCosPotentials.begin(), m_upwardsAccumCosPotentials.end(), 0.0);
    std::fill(m_upwardsAccumSinPotentials.begin(), m_upwardsAccumSinPotentials.end(), 0.0);
    std::fill(m_downwardsAccumCosPotentials.begin(), m_downwardsAccumCosPotentials.end(), 0.0);
    std::fill(m_downwardsAccumSinPotentials.begin(), m_downwardsAccumSinPotentials.end(), 0.0);
}

void Cube::SetUpwardsAccumCosPotentials(const std::vector<double>& upwardsAccumCosPotentials) {
    m_upwardsAccumCosPotentials = upwardsAccumCosPotentials;
}

const std::vector<double>& Cube::GetUpwardsAccumCosPotentials() const { return m_upwardsAccumCosPotentials; }

void Cube::SetUpwardsAccumSinPotentials(const std::vector<double>& upwardsAccumSinPotentials) {
    m_upwardsAccumSinPotentials = upwardsAccumSinPotentials;
}

const std::vector<double>& Cube::GetUpwardsAccumSinPotentials() const { return m_upwardsAccumSinPotentials; }

void Cube::SetDownwardsAccumCosPotentials(const std::vector<double>& downwardsAccumCosPotentials) {
    m_downwardsAccumCosPotentials = downwardsAccumCosPotentials;
}

const std::vector<double>& Cube::GetDownwardsAccumCosPotentials() const { return m_downwardsAccumCosPotentials; }

void Cube::SetDownwardsAccumSinPotentials(const std::vector<double>& downwardsAccumSinPotentials) {
    m_downwardsAccumSinPotentials = downwardsAccumSinPotentials;
}

const std::vector<double>& Cube::GetDownwardsAccumSinPotentials() const { return m_downwardsAccumSinPotentials; }

// private procs

void Cube::AllocateAccumulatedSourcePotentials(std::size_t numLevelsInHeirarchyOfCubes, std::size_t numRadialQuadPts) {
    std::size_t numElements = numLevelsInHeirarchyOfCubes * numRadialQuadPts;

    m_upwardsAccumCosPotentials.resize(numElements, 0.0);
    m_upwardsAccumSinPotentials.resize(numElements, 0.0);
    m_downwardsAccumCosPotentials.resize(numElements, 0.0);
    m_downwardsAccumSinPotentials.resize(numElements, 0.0);
}
