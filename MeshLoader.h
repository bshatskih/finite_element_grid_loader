#pragma once
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <ranges>
#include "Node.h"
#include "FiniteElement.h"

class MeshLoader {
 protected:
  uint32_t NodesNumber = 0;
  uint32_t Dimension = 0;

  uint32_t FiniteElemNumber = 0;
  uint32_t FiniteElemSize = 0;

  uint32_t PlaneFiniteElemNumber = 0;
  uint32_t PlaneFiniteElemSize = 0;

  std::vector<Node> m_Nodes;
  std::vector<FiniteElement> m_FiniteElements;
  std::vector<FiniteElement> m_PlaneFiniteElements;

 public:

  MeshLoader() = default;

  virtual ~MeshLoader() = default;

  virtual void LoadMesh(const std::string &path_to_file) = 0;

  [[nodiscard]] std::vector<uint32_t> FindFiniteElementsByVertices(uint32_t id1, uint32_t id2, uint32_t id3) const;

  [[nodiscard]] std::vector<uint32_t> FindFiniteElementsByEdge(uint32_t id1, uint32_t id2) const;


  [[nodiscard]] uint32_t FindFiniteElementsByEdgeBySurfaceNodes(const std::vector<uint32_t>& node_id) const;

  [[nodiscard]] const std::vector<Node>& GetNodes() const;

  [[nodiscard]] const std::vector<FiniteElement>& GetFiniteElements() const;

  [[nodiscard]] const std::vector<FiniteElement>& GetPlaneFiniteElements() const;

  [[nodiscard]] std::unordered_set<uint32_t> GetNodesFromSurface(uint32_t surface_id) const;

  [[nodiscard]] std::vector<uint32_t> GetPlaneFiniteElementsFromSurface(uint32_t surface_id) const;

  [[nodiscard]] std::vector<uint32_t> GetFiniteElementsFromMaterial(uint32_t material_id) const;

  [[nodiscard]] std::unordered_map<uint32_t, std::unordered_set<uint32_t>> GetAllNeighbors() const;

  void BisectAllEdges();

  [[nodiscard]] double MinAngel(const FiniteElement& finite_element) const;

  [[nodiscard]] double VolumeFiniteElement(const FiniteElement& finite_element) const;

  [[nodiscard]] double SquareSurfaceFiniteElement(const FiniteElement& finite_element) const;

  [[nodiscard]] std::vector<uint32_t> GetFiniteVertexFiniteElement(const FiniteElement& finite_element) const;

  friend std::ostream &operator<<(std::ostream &os, const MeshLoader &mesh_loader);
};


