#include "AneuMeshLoader.h"
#include <iostream>

enum ReadMode {
  None = 0,
  ReadNode = 1,
  ReadFiniteElem = 2,
  ReadPlaneFiniteElem = 3
};


void AneuMeshLoader::LoadMesh(const std::string &path_to_file) {
  m_File.open(path_to_file, std::ios_base::in);
  if (!m_File.is_open()) {
    std::cout << "DUMB)))";
  }
  std::string line{};
  short mode = ReadMode::None;

  while (std::getline(m_File, line)) {
    if (not line.starts_with(' ')) {
      mode++;
      size_t spase_pos = line.find(' ');
      int FirstNumber = std::stoi(line.substr(0, spase_pos));
      int SecondNumber = std::stoi(line.substr(spase_pos + 1, line.length()));
      switch (mode) {
        case ReadMode::ReadNode: {
          NodesNumber = FirstNumber;
          Dimension = SecondNumber;
          LoadNodes();
          break;
        }
        case ReadMode::ReadFiniteElem: {
          FiniteElemNumber = FirstNumber;
          FiniteElemSize = SecondNumber;
          LoadFiniteElements();
          break;
        }
        case ReadMode::ReadPlaneFiniteElem: {
          PlaneFiniteElemNumber = FirstNumber;
          PlaneFiniteElemSize = SecondNumber;
          LoadPlateFiniteElements();
          break;
        }
        default:
          break;
      }
    }
  }
}


void AneuMeshLoader::LoadNodes() {
  std::string line{};
  m_Nodes.reserve(NodesNumber);
  for (int i = 0; i < NodesNumber; ++i) {
    std::getline(m_File, line);
    m_Nodes.emplace_back(line, Dimension, i + 1);
  }
}


void AneuMeshLoader::LoadFiniteElements() {
  std::string line{};
  m_FiniteElements.reserve(FiniteElemNumber);
  for (int i = 0; i < FiniteElemNumber; ++i) {
    std::getline(m_File, line);
    m_FiniteElements.emplace_back(line, FiniteElemSize, i + 1);
  }
}


void AneuMeshLoader::LoadPlateFiniteElements() {
  std::string line{};
  m_PlaneFiniteElements.reserve(PlaneFiniteElemNumber);
  for (int i = 0; i < PlaneFiniteElemNumber; ++i) {
    std::getline(m_File, line);
    m_PlaneFiniteElements.emplace_back(line, PlaneFiniteElemSize, i + 1);
    m_PlaneFiniteElements.back().ID = FindFiniteElementsByEdgeBySurfaceNodes(m_PlaneFiniteElements.back().NodeIDs);
  }
}


