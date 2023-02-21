#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <filesystem>

#include "readInputFiles.h"

void Tokenize(std::vector<std::string> &fields, const std::string &line) {
    fields.clear();
    std::string delim = " ";
    std::string lineToSplit = line;
    std::string field;

    std::size_t start = line.find_first_not_of(delim, 0);
    std::size_t end = line.find_first_of(delim, start);
    while ((start <= line.size()) && (end <= line.size())) {
        fields.push_back(line.substr(start, end - start));
        start = line.find_first_not_of(delim, end);
        if (start == std::string::npos) break;
        end = line.find_first_of(delim, start);
        if (end == std::string::npos) {
            end = line.size();
            fields.push_back(line.substr(start, end - start));
            break;
        }
    }
}

bool InputAtomInfoFromWaterBox(std::vector<std::string> &names, std::vector<double> &charges,
                               std::vector<double> &xCoords, std::vector<double> &yCoords, std::vector<double> &zCoords,
                               std::vector<std::size_t> &indicesInVecOfAtoms,
                               std::vector<std::vector<std::size_t>> &exclusionLists, const std::string &fileName) {
    names.clear();
    charges.clear();
    xCoords.clear();
    yCoords.clear();
    zCoords.clear();
    indicesInVecOfAtoms.clear();
    std::vector<std::vector<std::size_t>>::iterator iterE = exclusionLists.begin();
    std::vector<std::vector<std::size_t>>::iterator endE = exclusionLists.end();
    for (; iterE != endE; ++iterE) iterE->clear();
    exclusionLists.clear();

    std::ifstream ifs;
    ifs.open(fileName, std::ifstream::in);
    if (!ifs.is_open()) return false;

    std::string line;
    if (!getline(ifs, line)) return false;
    std::cout << line << std::endl;
    std::vector<std::string> fields;
    Tokenize(fields, line);
    std::size_t num = std::stoi(fields[0]);

    names.reserve(num);
    charges.reserve(num);
    xCoords.reserve(num);
    yCoords.reserve(num);
    zCoords.reserve(num);
    indicesInVecOfAtoms.reserve(num);
    exclusionLists.reserve(num);

    std::size_t numWaters = num / 3;
    if (numWaters == 0) numWaters = 1;

    std::cout << "numWaters = " << numWaters << std::endl;

    for (std::size_t n = 0; n < numWaters; ++n) {
        std::vector<std::size_t> exclusionList;
        exclusionList.push_back(3 * n);
        exclusionList.push_back(3 * n + 1);
        exclusionList.push_back(3 * n + 2);

        if (!getline(ifs, line)) return false;
        if (line == "continue") continue;
        Tokenize(fields, line);
        if (fields.size() < 5) return false;
        std::string name = fields[1];
        double charge = 0.0;
        if (name == "O") {
            charge = -0.82;  // SPC oxygen
        } else if (name == "H") {
            charge = 0.41;  // SPC hydrogen
        } else if (name == "N") {
            charge = 1.00;  // SPC hydrogen
        } else if (name == "C") {
            charge = 0.00;  // SPC hydrogen
        }
        double xCoord = std::stod(fields[2]);
        double yCoord = std::stod(fields[3]);
        double zCoord = std::stod(fields[4]);

        names.push_back(name);
        charges.push_back(charge);
        xCoords.push_back(xCoord);
        yCoords.push_back(yCoord);
        zCoords.push_back(zCoord);
        std::size_t indexInListOfAtoms = 3 * n;
        indicesInVecOfAtoms.push_back(indexInListOfAtoms);
        exclusionLists.push_back(exclusionList);

        if (!getline(ifs, line)) return false;
        if (line == "continue") continue;
        Tokenize(fields, line);
        if (fields.size() < 5) return false;
        name = fields[1];
        if (name == "O") {
            charge = -0.82;  // SPC oxygen
        } else if (name == "H") {
            charge = 0.41;  // SPC hydrogen
        } else if (name == "N") {
            charge = 1.00;  // SPC hydrogen
        } else if (name == "C") {
            charge = 0.00;  // SPC hydrogen
        }
        indexInListOfAtoms = 3 * n + 1;
        xCoord = std::stod(fields[2]);
        yCoord = std::stod(fields[3]);
        zCoord = std::stod(fields[4]);

        names.push_back(name);
        charges.push_back(charge);
        xCoords.push_back(xCoord);
        yCoords.push_back(yCoord);
        zCoords.push_back(zCoord);
        indicesInVecOfAtoms.push_back(indexInListOfAtoms);
        exclusionLists.push_back(exclusionList);

        if (!getline(ifs, line)) return false;
        if (line == "continue") continue;
        Tokenize(fields, line);
        if (fields.size() < 5) return false;
        name = fields[1];
        if (name == "O") {
            charge = -0.82;  // SPC oxygen
        } else if (name == "H") {
            charge = 0.41;  // SPC hydrogen
        } else if (name == "N") {
            charge = 1.00;  // SPC hydrogen
        } else if (name == "C") {
            charge = 0.00;  // SPC hydrogen
        }
        indexInListOfAtoms = 3 * n + 2;
        xCoord = std::stod(fields[2]);
        yCoord = std::stod(fields[3]);
        zCoord = std::stod(fields[4]);

        names.push_back(name);
        charges.push_back(charge);
        xCoords.push_back(xCoord);
        yCoords.push_back(yCoord);
        zCoords.push_back(zCoord);
        indicesInVecOfAtoms.push_back(indexInListOfAtoms);
        exclusionLists.push_back(exclusionList);
    }
    ifs.close();

    return true;
}

void GetFileNameInDirectory(std::string &fileName, std::size_t &numPts, unsigned degree,
                            const std::string &pathToDirectory) {
    fileName = "";
    numPts = 0;
    std::stringstream ss;
    ss << degree;
    std::string pattern = "ss" + ss.str();
    if (degree < 100) pattern = "ss0" + ss.str();
    std::cout << "pattern = " << pattern << std::endl;
    for (const auto &file : std::filesystem::directory_iterator{pathToDirectory}) {
        std::string name = file.path();
        std::string subName = name.substr(pathToDirectory.size() + 1, 5);
        std::string endName = name.substr(pathToDirectory.size() + 7, std::string::npos);

        if (pattern == subName) {
            numPts = static_cast<std::size_t>(std::stoi(endName));
            numPts /= 2;  // only take upper hemisphere
            fileName = name;
            return;
        }
    }
}

bool LoadSpherePoints(std::vector<double> &sx, std::vector<double> &sy, std::vector<double> &sz, std::size_t numPts,
                      const std::string &fileName) {
    sx.clear();
    sy.clear();
    sz.clear();
    std::ifstream ifs;
    ifs.open(fileName, std::ifstream::in);
    if (!ifs.is_open()) return false;

    std::string line;
    std::vector<std::string> fields;

    sx.reserve(numPts);
    sy.reserve(numPts);
    sz.reserve(numPts);
    for (int n = 0; n < numPts; ++n) {
        if (!getline(ifs, line)) return false;
        Tokenize(fields, line);
        if (fields.size() != 3) return false;
        sx.push_back(std::stod(fields[0]));
        sy.push_back(std::stod(fields[1]));
        sz.push_back(std::stod(fields[2]));
    }

    ifs.close();
    return true;
}

bool InputSpherePoints(std::vector<double> &sx, std::vector<double> &sy, std::vector<double> &sz, unsigned degree,
                       const std::string &pathToDirectory) {
    std::string fileName;
    std::size_t numPts = 0;
    GetFileNameInDirectory(fileName, numPts, degree, pathToDirectory);
    if (numPts == 0) return false;
    if (!LoadSpherePoints(sx, sy, sz, numPts, fileName)) return false;
    return true;
}
