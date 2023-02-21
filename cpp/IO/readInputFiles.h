#ifndef READ_INPUT_FILES_H
#define READ_INPUT_FILES_H

#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>

void Tokenize(std::vector<std::string> &fields, const std::string &line);

void GetFileNameInDirectory(std::string &fileName, unsigned degree, const std::string &pathToDirectory);

bool InputAtomInfoFromWaterBox(std::vector<std::string> &names, std::vector<double> &charges,
                               std::vector<double> &xCoords, std::vector<double> &yCoords, std::vector<double> &zCoords,
                               std::vector<std::size_t> &indicesInVecOfAtoms,
                               std::vector<std::vector<std::size_t>> &exclusionLists, const std::string &fileName);

void GetFileNameInDirectory(std::string &fileName, std::size_t &numPts, unsigned degree,
                            const std::string &pathToDirectory);

bool LoadSpherePoints(std::vector<double> &sx, std::vector<double> &sy, std::vector<double> &sz, std::size_t numPts,
                      const std::string &fileName);

bool InputSpherePoints(std::vector<double> &sx, std::vector<double> &sy, std::vector<double> &sz, unsigned degree,
                       const std::string &pathToDirectory);

#endif  // READ_INPUT_FILES_H
