#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include "../Atom/atom.h"
#include "../Atom/utility.h"
#include "../IO/readInputFiles.h"

int main(int argc, char *argv[]) {
    if (argc != 6) {
        std::cout << "input filename for waters, output file and x,y,z translations for extra atoms" << std::endl;
        exit(1);
    }

    std::string waterBoxName = argv[1];
    std::string outputName = argv[2];
    double transX = std::stod(argv[3]);
    double transY = std::stod(argv[4]);
    double transZ = std::stod(argv[5]);

    std::vector<Atom> vecOfAtoms;
    if (!InputVecOfAtomsFromWaterBox(vecOfAtoms, waterBoxName)) {
        std::cout << "failed to read in atoms from waterBoxName" << std::endl;
        exit(1);
    }
    std::cout << "num atoms before = " << vecOfAtoms.size() << std::endl;

    double minDistance = std::numeric_limits<double>::max();

    if (!GetSmallestDistanceBetweenWatersInMolecule(minDistance, vecOfAtoms)) {
        std::cout << "GetSmallestDistanceBetweenWatersInMolecule fails!" << std::endl;
        exit(1);
    }
    std::cout << "smallest distance between waters in mol before = " << minDistance << std::endl;

    std::vector<Atom> translatedAtoms = vecOfAtoms;

    std::size_t index = vecOfAtoms.size();

    for (auto atom : vecOfAtoms) {
        std::string name = atom.GetName();
        double q = atom.GetCharge();
        double x = atom.GetXCoord();
        double y = atom.GetYCoord();
        double z = atom.GetZCoord();

        double newX = x + transX;
        double newY = y + transY;
        double newZ = z + transZ;

        std::vector<std::size_t> exclList = atom.GetExclusionList();
        std::vector<std::size_t> newExclList;
        for (auto l : exclList) newExclList.push_back(l + vecOfAtoms.size());

        Atom newAtom(name, q, newX, newY, newZ, index, newExclList);
        translatedAtoms.push_back(newAtom);
        ++index;
    }

    std::cout << "num atoms after = " << translatedAtoms.size() << std::endl;

    minDistance = std::numeric_limits<double>::max();

    if (!GetSmallestDistanceBetweenWatersInMolecule(minDistance, translatedAtoms)) {
        std::cout << "GetSmallestDistanceBetweenWatersInMolecule fails!" << std::endl;
        exit(1);
    }
    std::cout << "smallest distance between waters in mol after = " << minDistance << std::endl;

    if (!OutputVecOfAtoms(translatedAtoms, outputName)) {
        std::cout << "failed to output atoms" << std::endl;
        exit(1);
    }

    return 0;
}
