#ifndef FINITE_ELEMENT_GRID_LOADER__ANEUMESHLOADER_H_
#define FINITE_ELEMENT_GRID_LOADER__ANEUMESHLOADER_H_
#include <fstream>
#include "MeshLoader.h"
#include "Node.h"
#include "FiniteElement.h"

class AneuMeshLoader : public MeshLoader {
 public:
  AneuMeshLoader() = default;
  void LoadMesh(const std::string &path_to_file) override;

 private:
  std::ifstream m_File;

 private:

  void LoadNodes();
  void LoadFiniteElements();
  void LoadPlateFiniteElements();
};

#endif //FINITE_ELEMENT_GRID_LOADER__ANEUMESHLOADER_H_
