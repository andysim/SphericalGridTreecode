#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../Atom/atom.h"
#include "../Atom/utility.h"

#include "../IO/readInputFiles.h"
#include "../Quadrature/hermiteRootsWeights.h"

#include "../Hierarchy/hierarchy.h"

bool InputQuadraturePoints(std::vector<double> &gaussHermiteRoots, std::vector<double> &gaussHermiteWeights,
                           std::vector<double> &sphPtsX, std::vector<double> &sphPtsY, std::vector<double> &sphPtsZ,
                           double &ewaldCoeff, double THRESH, double directSumCutoff,
                           const std::string &pathFromHereToSpherePts) {
    gaussHermiteRoots.clear();
    gaussHermiteWeights.clear();
    sphPtsX.clear();
    sphPtsY.clear();
    sphPtsZ.clear();

    if (!GetTrimmedGaussHermiteRoots_and_Weights(gaussHermiteRoots, gaussHermiteWeights, ewaldCoeff, THRESH,
                                                 directSumCutoff)) {
        std::cout << "GetTrimmedGaussHermiteRoots_and_Weights failed" << std::endl;
        return false;
    }

    unsigned degree = GetDegreeOfSphericalDesignFromHermiteRootsWeights(THRESH, directSumCutoff, ewaldCoeff,
                                                                        gaussHermiteRoots, gaussHermiteWeights);

    std::cout << "degree = " << degree << std::endl;

    if (!InputSpherePoints(sphPtsX, sphPtsY, sphPtsZ, degree, pathFromHereToSpherePts)) {
        std::cout << "failed to read in sphere pts" << std::endl;
        return false;
    }
    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        std::cout << "arguments: " << std::endl
                  << std::endl
                  << "accuracy_threshold  cutoff  numprocessors "
                  << "relative_path_to_the_directory_containing_sphere_pt_files  "
                  << "waterbox_xyz" << std::endl;
        exit(1);
    }

    // we aim for accuracy of 10^(-THRESH)
    // directSumCutoff is the cutoff for the erfc calculation that complements the erf calculation
    // the input for calculation is a waterbox containing N point charges
    double THRESH = std::stod(argv[1]);
    double directSumCutoff = std::stod(argv[2]);
    int nprocs = std::stoi(argv[3]);
    std::string pathFromHereToSpherePts = argv[4];
    std::string waterBoxName = argv[5];

    //--------------------------------------------------------------//
    // Input part //

    std::vector<double> gaussHermiteRoots;
    std::vector<double> gaussHermiteWeights;

    std::vector<double> sphPtsX;
    std::vector<double> sphPtsY;
    std::vector<double> sphPtsZ;

    double ewaldCoeff = 0.0;

    if (!InputQuadraturePoints(gaussHermiteRoots, gaussHermiteWeights, sphPtsX, sphPtsY, sphPtsZ, ewaldCoeff, THRESH,
                               directSumCutoff, pathFromHereToSpherePts)) {
        std::cout << "InputQuadraturePoints	failed!" << std::endl;
        exit(1);
    }
    std::cout << "size of sphPtsX = " << sphPtsX.size() << std::endl;

    std::vector<Atom> vecOfAtoms;
    if (!InputVecOfAtomsFromWaterBox(vecOfAtoms, waterBoxName)) {
        std::cout << "failed to read in atoms from waterBoxName" << std::endl;
        exit(1);
    }
    std::cout << "num atoms = " << vecOfAtoms.size() << std::endl;

    double lowestXCoordOfMol = 0.0;
    double lowestYCoordOfMol = 0.0;
    double lowestZCoordOfMol = 0.0;

    double highestXCoordOfMol = 0.0;
    double highestYCoordOfMol = 0.0;
    double highestZCoordOfMol = 0.0;

    double cubeWidth = directSumCutoff;

    GetLowerUpperBoundsOnAtomCoords(lowestXCoordOfMol, highestXCoordOfMol, lowestYCoordOfMol, highestYCoordOfMol,
                                    lowestZCoordOfMol, highestZCoordOfMol, vecOfAtoms);

    AdjustMaxMinOfMolForCubeWidth(lowestXCoordOfMol, highestXCoordOfMol, lowestYCoordOfMol, highestYCoordOfMol,
                                  lowestZCoordOfMol, highestZCoordOfMol, cubeWidth);

    // end input part //
    //--------------------------------------------------------------//
    // start of calc over vec of Atoms

    auto start1 = std::chrono::steady_clock::now();
    std::vector<Atom>::iterator iterA = vecOfAtoms.begin();
    std::vector<Atom>::iterator endA = vecOfAtoms.end();
    for (; iterA != endA; ++iterA) iterA->LoadExactCoulombPotential_and_Grad(vecOfAtoms);
    auto end1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds1 = end1 - start1;
    std::cout << "elapsed time to load coulomb: " << elapsed_seconds1.count() << "s\n";

    auto start2 = std::chrono::steady_clock::now();
    iterA = vecOfAtoms.begin();
    for (; iterA != endA; ++iterA)
        iterA->LoadCutoffErfc_and_ExclusionListErf_Potential_and_Grad(vecOfAtoms, directSumCutoff, ewaldCoeff);
    auto end2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds2 = end2 - start2;
    std::cout << "elapsed time to load erfc, excluded erf pots: " << elapsed_seconds2.count() << "s\n";

    // end of calc over atoms
    //--------------------------------------------------------------//
    // note below we are creating a vector of hierarchies --one per "processor" --this is
    // to emulate mpi --at least to show the embarrasingly parallel nature of the
    // sum over quadrature points
    //--------------------------------------------------------------//

    double r12Max = 4.0 * sqrt(3.0) * directSumCutoff;
    double maxDistance = GetLargestDistanceInMolecule(vecOfAtoms);
    std::cout << "r12Max = " << r12Max << std::endl;
    std::cout << "largest atom-atom distance in molecule = " << maxDistance << std::endl;

    size_t numRadialQuadPts = gaussHermiteRoots.size();
    size_t numSpherePoints = sphPtsX.size();
    std::cout << "numRadialQuadPts = " << numRadialQuadPts << std::endl;
    std::cout << "numSpherePts = " << numSpherePoints << std::endl;

    std::vector<Hierarchy> vecOfHierarchies;
    vecOfHierarchies.reserve(nprocs);

    for (size_t n = 0; n < nprocs; ++n) {
        Hierarchy hierarchy(cubeWidth, lowestXCoordOfMol, lowestYCoordOfMol, lowestZCoordOfMol, ewaldCoeff,
                            numSpherePoints, gaussHermiteRoots, gaussHermiteWeights, highestXCoordOfMol,
                            highestYCoordOfMol, highestZCoordOfMol, vecOfAtoms);

        if (!hierarchy.IsValid()) {
            std::cout << "hierarchy not valid!!" << std::endl;
            exit(1);
        }
        vecOfHierarchies.push_back(hierarchy);
    }
    std::cout << "done with hierarchy construction" << std::endl;

    auto start3 = std::chrono::steady_clock::now();
    //--------------------------------------------------------------//
    // start of sum over sphere pts of approx erf
    #pragma omp parallel num_threads(nprocs)
    {
        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif        
        for (size_t m = rank; m < numSpherePoints; m += nprocs) {
            double sphPtx = sphPtsX[m];
            double sphPty = sphPtsY[m];
            double sphPtz = sphPtsZ[m];

            if (!vecOfHierarchies[rank].AccumulatePotentialsForNewSphericalQuadraturePoint(sphPtx, sphPty, sphPtz)) {
                std::cout << "AccumulatePotentialsForNewSphericalQuadraturePoint failed!" << std::endl;
                exit(1);
            }
        }
    }
    //for (size_t m = 0; m < numSpherePoints; ++m) {
    //    double sphPtx = sphPtsX[m];
    //    double sphPty = sphPtsY[m];
    //    double sphPtz = sphPtsZ[m];

    //    size_t jproc = m % nprocs;
    //    if (!vecOfHierarchies[jproc].AccumulatePotentialsForNewSphericalQuadraturePoint(sphPtx, sphPty, sphPtz)) {
    //        std::cout << "AccumulatePotentialsForNewSphericalQuadraturePoint failed!" << std::endl;
    //        exit(1);
    //    }
    //}
    // end of sum over sphere pts of approx erf
    //--------------------------------------------------------------//
    // now get approx erf potentials

    for (size_t n = 0; n < nprocs; ++n) vecOfHierarchies[n].LoadApproxErfPotential_and_Grad();

    size_t sizeOfAtomList = vecOfAtoms.size();

    std::vector<double> approxErfPots(sizeOfAtomList, 0.0);
    std::vector<double> approxGradErfPot_wrt_X(sizeOfAtomList, 0.0);
    std::vector<double> approxGradErfPot_wrt_Y(sizeOfAtomList, 0.0);
    std::vector<double> approxGradErfPot_wrt_Z(sizeOfAtomList, 0.0);

    for (size_t n = 0; n < nprocs; ++n) {
        std::vector<double> partialApproxErfPots;
        std::vector<double> partialApproxGradErfPot_wrt_X;
        std::vector<double> partialApproxGradErfPot_wrt_Y;
        std::vector<double> partialApproxGradErfPot_wrt_Z;

        if (!vecOfHierarchies[n].GetVecsOfApproxErfPot_and_Grads(partialApproxErfPots, partialApproxGradErfPot_wrt_X,
                                                                 partialApproxGradErfPot_wrt_Y,
                                                                 partialApproxGradErfPot_wrt_Z, sizeOfAtomList)) {
            std::cout << "GetVecsOfApproxErfPot_and_Grads failed!" << std::endl;
            exit(1);
        }
        for (size_t i = 0; i < sizeOfAtomList; ++i) {
            approxErfPots[i] += partialApproxErfPots[i];
            approxGradErfPot_wrt_X[i] += partialApproxGradErfPot_wrt_X[i];
            approxGradErfPot_wrt_Y[i] += partialApproxGradErfPot_wrt_Y[i];
            approxGradErfPot_wrt_Z[i] += partialApproxGradErfPot_wrt_Z[i];
        }
    }

    auto end3 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds3 = end3 - start3;
    std::cout << "elapsed time to calculate approx erf: " << elapsed_seconds3.count() << "s\n";

    for (size_t i = 0; i < sizeOfAtomList; ++i) {
        double pot = approxErfPots[i];
        std::vector<double> grad(3, 0.0);
        grad[0] = approxGradErfPot_wrt_X[i];
        grad[1] = approxGradErfPot_wrt_Y[i];
        grad[2] = approxGradErfPot_wrt_Z[i];
        vecOfAtoms[i].SetApproxErfPotential_and_Grad(pot, grad);
        vecOfAtoms[i].SetApproxCoulombPotential_and_Grad();
    }

    if (vecOfAtoms.size() < 10) {
        std::cout << "potentials: " << std::endl;
        for (size_t i = 0; i < vecOfAtoms.size(); ++i) std::cout << vecOfAtoms[i].PrintPots() << std::endl;
    }

    std::cout << "max error over atoms of coulomb pot = " << GetMaxErrorInCoulombPotOverVecOfAtoms(vecOfAtoms)
              << std::endl;
    std::cout << "RMS error over atoms of coulomb pot = " << GetRMSErrorInCoulombPotOverVecOfAtoms(vecOfAtoms)
              << std::endl;
    std::cout << "RMS error over atoms of gradient of coulomb pot = "
              << GetRMSErrorInGradOfCoulombPotOverVecOfAtoms(vecOfAtoms) << std::endl;
    std::cout << "relative RMS error over atoms of gradient of coulomb pot = "
              << GetRelativeRMSErrorOfGradOfCoulombPotOverVecOfAtoms(vecOfAtoms) << std::endl;

    return 0;
}
