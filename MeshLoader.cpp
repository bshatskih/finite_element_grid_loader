#include <set>
#include "MeshLoader.h"
#include "PlainClass.h"

std::ostream &operator<<(std::ostream &os, const MeshLoader &mesh_loader) {
  os << mesh_loader.NodesNumber << ' ' << mesh_loader.Dimension << '\n';
  for (const Node &m_Node : mesh_loader.m_Nodes) {
    os << std::setprecision(6) << std::fixed << m_Node << '\n';
  }
  os << mesh_loader.FiniteElemNumber << ' ' << mesh_loader.FiniteElemSize << '\n';
  for (const FiniteElement &finite_element : mesh_loader.m_FiniteElements) {
    os << finite_element << '\n';
  }
  os << mesh_loader.PlaneFiniteElemNumber << ' ' << mesh_loader.PlaneFiniteElemSize << '\n';
  for (const FiniteElement &plate_finite_element : mesh_loader.m_PlaneFiniteElements) {
    os << plate_finite_element << '\n';
  }
  return os;
}

const std::vector<Node> &MeshLoader::GetNodes() const {
  return m_Nodes;
}

const std::vector<FiniteElement> &MeshLoader::GetFiniteElements() const {
  return m_FiniteElements;
}

const std::vector<FiniteElement> &MeshLoader::GetPlaneFiniteElements() const {
  return m_PlaneFiniteElements;
}


///first initialization without functors and predicates
//std::vector<uint32_t> MeshLoader::FindFiniteElementsByVertices(uint32_t id1, uint32_t id2, uint32_t id3) const {
//  std::vector<uint32_t> result;
//  for (std::vector<FiniteElement>::const_iterator it = m_FiniteElements.cbegin(); it != m_FiniteElements.cend(); it++) {
//    std::vector<uint32_t>::const_iterator beg = it->NodeIDs.cbegin();
//    std::vector<uint32_t>::const_iterator end = it->NodeIDs.cend();
//    bool flag = std::find(beg, end, id1) != end && std::find(beg, end, id2) != end && std::find(beg, end, id3) != end;
//    if (flag) {
//      result.push_back(it->ID);
//    }
//  }
//  return result;
//}

///initialization option with functor
//struct FindFiniteElementsByVertices_functor {
//  const uint32_t ID1;
//  const uint32_t ID2;
//  const uint32_t ID3;
//
//  FindFiniteElementsByVertices_functor(uint32_t id1, uint32_t id2, uint32_t id3) : ID1(id1), ID2(id2), ID3(id3) {}
//
//  bool operator()(const FiniteElement &finite_element) {
//    auto beg = finite_element.NodeIDs.cbegin();
//    auto end = finite_element.NodeIDs.cend();
//    return std::find(beg, end, ID1) != end && std::find(beg, end, ID2) != end && std::find(beg, end, ID3) != end;
//  }
//};
//
//std::vector<uint32_t> MeshLoader::FindFiniteElementsByVertices(uint32_t id1, uint32_t id2, uint32_t id3) const {
//  std::vector<uint32_t> result;
//  auto beg = m_FiniteElements.begin();
//  auto end = m_FiniteElements.end();
//
//  FindFiniteElementsByVertices_functor functor(id1, id2, id3);
//  bool flag;
//  while (beg != end) {
//    beg = std::find_if(beg, end, functor);
//    if (beg != end) {
//      result.push_back(beg->ID);
//      ++beg;
//    }
//  }
//  return result;
//}

///initialization option with predicate
std::vector<uint32_t> MeshLoader::FindFiniteElementsByVertices(uint32_t id1, uint32_t id2, uint32_t id3) const {
  std::vector<uint32_t> result;
  auto beg = m_FiniteElements.begin();
  auto end = m_FiniteElements.end();

  while (beg != end) {
    beg = std::find_if(beg, end, [=](const FiniteElement &finite_element) {
      auto b = finite_element.NodeIDs.cbegin();
      auto e = finite_element.NodeIDs.cend();
      return std::find(b, e, id1) != e && std::find(b, e, id2) != e && std::find(b, e, id3) != e;
    });
    if (beg != end) {
      result.push_back(beg->ID);
      ++beg;
    }
  }
  return result;
}

uint32_t MeshLoader::FindFiniteElementsByEdgeBySurfaceNodes(const std::vector<uint32_t> &node_id) const {
  auto pos = std::find_if(m_FiniteElements.cbegin(), m_FiniteElements.cend(), [&](const FiniteElement &finite_element) {
    bool res = true;
    for (unsigned int id : node_id) {
      res &= std::find(finite_element.NodeIDs.cbegin(), finite_element.NodeIDs.cend(), id)
          != finite_element.NodeIDs.cend();
    }
    return res;
  });
  return pos->ID;
}

std::vector<uint32_t> MeshLoader::FindFiniteElementsByEdge(uint32_t id1, uint32_t id2) const {
  std::vector<uint32_t> result;
  auto beg = m_FiniteElements.begin();
  auto end = m_FiniteElements.end();

  while (beg != end) {
    beg = std::find_if(beg, end, [=](const FiniteElement &finite_element) {
      auto b = finite_element.NodeIDs.cbegin();
      auto e = finite_element.NodeIDs.cend();
      return std::find(b, e, id1) != e && std::find(b, e, id2) != e;
    });
    if (beg != end) {
      result.push_back(beg->ID);
      ++beg;
    }
  }
  return result;
}



///initialization option number one
//std::unordered_set<uint32_t> MeshLoader::GetNodesFromSurface(uint32_t surface_id) const {
//  std::unordered_set<uint32_t> result_set;
//  std::vector<Node> result_vec;
//  auto beg = m_PlaneFiniteElements.cbegin();
//  auto end = m_PlaneFiniteElements.cend();
//  while (beg != end) {
//    beg = std::find_if(beg, end, [&](const FiniteElement& finite_element) {
//      return finite_element.ID == surface_id;
//    });
//    if (beg != end) {
//      result_set.insert(beg->NodeIDs.cbegin(), beg->NodeIDs.cend());
//      std::copy_if(m_Nodes.cbegin(), m_Nodes.cend(), std::back_inserter(result_vec), [&](const Node& node) {
//        return std::find(beg->NodeIDs.cbegin(), beg->NodeIDs.cend(), node.ID) != beg->NodeIDs.cend();
//      });
//      ++beg;
//    }
//  }
//}

///initialization option number two
std::unordered_set<uint32_t> MeshLoader::GetNodesFromSurface(uint32_t surface_id) const {
  std::unordered_set<uint32_t> result_set;
  std::for_each(m_PlaneFiniteElements.cbegin(), m_PlaneFiniteElements.cend(), [&](const FiniteElement &finite_element) {
    if (finite_element.PositionID == surface_id) {
      result_set.insert(finite_element.NodeIDs.cbegin(), finite_element.NodeIDs.cend());
    }
  });
  return result_set;
}

std::vector<uint32_t> MeshLoader::GetPlaneFiniteElementsFromSurface(uint32_t surface_id) const {
  std::vector<uint32_t> result_vec;
  std::for_each(m_PlaneFiniteElements.cbegin(), m_PlaneFiniteElements.cend(), [&](const FiniteElement &finite_element) {
    if (finite_element.PositionID == surface_id) {
      result_vec.push_back(finite_element.ID);
    }
  });
  return result_vec;
}

std::vector<uint32_t> MeshLoader::GetFiniteElementsFromMaterial(uint32_t material_id) const {
  std::vector<uint32_t> result_vec;
  std::for_each(m_FiniteElements.cbegin(), m_FiniteElements.cend(), [&](const FiniteElement &finite_element) {
    if (finite_element.PositionID == material_id) {
      result_vec.push_back(finite_element.ID);
    }
  });
  return result_vec;
}

///not ideal way
//std::unordered_map<uint32_t, std::unordered_set<uint32_t>> MeshLoader::GetAllNeighbors() const {
//  std::unordered_map<uint32_t, std::unordered_set<uint32_t>> result_map;
//  result_map.reserve(m_Nodes.size());
//  std::for_each(m_Nodes.cbegin(), m_Nodes.cend(), [&](const Node &node) {
//    std::for_each(m_FiniteElements.cbegin(), m_FiniteElements.cend(), [&](const FiniteElement &finite_element) {
//      if (std::find(finite_element.NodeIDs.cbegin(), finite_element.NodeIDs.cend(), node.ID) != finite_element.NodeIDs.cend()) {
//        result_map[node.ID].insert(finite_element.NodeIDs.cbegin(), finite_element.NodeIDs.cend());
//      }
//    });
//    result_map[node.ID].erase(node.ID);
//  });
//  return result_map;
//}

///more better, but i don't think, that is perfect realization of this method
std::unordered_map<uint32_t, std::unordered_set<uint32_t>> MeshLoader::GetAllNeighbors() const {
  std::unordered_map<uint32_t, std::unordered_set<uint32_t>> result_map;
  result_map.reserve(m_Nodes.size());
  std::for_each(m_FiniteElements.cbegin(), m_FiniteElements.cend(), [&](const FiniteElement &finite_element) {
    std::for_each(finite_element.NodeIDs.cbegin(), finite_element.NodeIDs.cend(), [&](const uint32_t &node) {
      result_map[node].insert(finite_element.NodeIDs.cbegin(), finite_element.NodeIDs.cend());
    });
  });
  std::for_each(result_map.begin(),
                result_map.end(),
                [](std::pair<const uint32_t, std::unordered_set<uint32_t>> &nodes) {
                  nodes.second.erase(nodes.first);
                });
  return result_map;
}

///The problem that arises when implementing the method of inserting new nodes
/// into the new_vertex of an existing edge is that the same edge can belong to several
/// finite elements, which prompts us to somehow save information about which edge
/// we have already added a new node to.


template<typename T>
void HashCombine(size_t &seed, const T &vertex) {
  seed ^= std::hash<T>()(vertex) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template<typename T>
void HashVal(size_t &seed, const T &vertex) {
  HashCombine(seed, vertex);
}

template<typename T, typename ... Types>
void HashVal(size_t &seed, const T &vertex, const Types &... vertices) {
  HashCombine(seed, vertex);
  HashVal(seed, vertices ...);
};

template<typename ... Types>
size_t HashVal(const Types &... vertices) {
  size_t seed = 0;
  HashVal(seed, vertices ...);
  return seed;
};

struct new_node {
  uint32_t vertex_1;
  uint32_t vertex_2;
  uint32_t new_vertex;

  new_node(uint32_t id1, uint32_t id2, uint32_t mid) : vertex_1(id1), vertex_2(id2), new_vertex(mid) {}

  bool operator==(const new_node &other) const {
    return (vertex_1 == other.vertex_1 && vertex_2 == other.vertex_2)
        || (vertex_2 == other.vertex_1 && vertex_1 == other.vertex_2);
  }

  friend std::ostream &operator<<(std::ostream &os, const new_node &other) {
    os << "vertex_1: " << other.vertex_1 << " vertex_2: " << other.vertex_2 << " new_vertex: " << other.new_vertex;
    return os;
  }
};

class NewNodeEqual {
 public:
  bool operator()(const new_node &node1, const new_node &node2) const {
    return node1 == node2;
  }
};

class NewNodeHash {
 public:
  size_t operator()(const new_node &node) const {
    uint32_t vertex_1 = node.vertex_1;
    uint32_t vertex_2 = node.vertex_2;

    if (vertex_1 > vertex_2) {
      std::swap(vertex_1, vertex_2);
    }
    return HashVal(vertex_1, vertex_2);
  }
};

void MeshLoader::BisectAllEdges() {
  FiniteElemSize += (FiniteElemSize - 1) * FiniteElemSize / 2;
  PlaneFiniteElemSize += (PlaneFiniteElemSize - 1) * PlaneFiniteElemSize / 2;
  std::unordered_set<new_node, NewNodeHash, NewNodeEqual> new_nodes;
  std::for_each(m_FiniteElements.begin(), m_FiniteElements.end(), [&](FiniteElement &finite_element) {
    finite_element.NodeIDs.reserve(FiniteElemSize);
    size_t size = finite_element.NodeIDs.size();
    for (int i = 0; i < size - 1; ++i) {
      for (int j = i + 1; j < size; ++j) {
        uint32_t id1 = finite_element.NodeIDs[i];
        uint32_t id2 = finite_element.NodeIDs[j];
        new_node node_to_insert(id1, id2, m_Nodes.size());
        auto [inserted_iter, result] = new_nodes.insert(node_to_insert);
        if (result) {
          auto new_vertex = m_Nodes.emplace_back(m_Nodes[id1 - 1], m_Nodes[id2 - 1], Dimension, m_Nodes.size() + 1);
          NodesNumber++;
          finite_element.NodeIDs.push_back(new_vertex.ID);
        } else {
          finite_element.NodeIDs.push_back(inserted_iter->new_vertex);
        }
      }
    }
  });
  std::for_each(m_PlaneFiniteElements.begin(), m_PlaneFiniteElements.end(), [&](FiniteElement &finite_element) {
    finite_element.NodeIDs.reserve(PlaneFiniteElemSize);
    size_t size = finite_element.NodeIDs.size();
    for (int i = 0; i < size - 1; ++i) {
      for (int j = i + 1; j < size; ++j) {
        uint32_t id1 = finite_element.NodeIDs[i];
        uint32_t id2 = finite_element.NodeIDs[j];
        new_node node_to_find(id1, id2, 0);
        std::unordered_set<new_node, NewNodeHash, NewNodeEqual>::iterator it = new_nodes.find(node_to_find);
        finite_element.NodeIDs.push_back(it->new_vertex);
      }
    }
  });
}

class PlainHesh {
 public:
  size_t operator()(const Plain &plain) const {
    return HashVal(plain.node_1.ID, plain.node_2.ID, plain.node_3.ID);
  }
};

class PlainEqual {
 public:
  bool operator()(const Plain &plain1, const Plain &plain2) const {
    return plain1 == plain2;
  }
};

double MeshLoader::MinAngel(const FiniteElement &finite_element) const {
  std::unordered_set<Plain, PlainHesh, PlainEqual> planes;
  for (const auto &id1 : finite_element.NodeIDs) {
    auto other_1 = finite_element.NodeIDs | std::ranges::views::filter([&](const auto &id) { return id != id1; });
    for (const auto &id2 : other_1) {
      auto other_2 = other_1 | std::ranges::views::filter([&](const auto &id) { return id != id2; });
      for (const auto &id3 : other_2) {
        uint32_t ids[3]{id1, id2, id3};
//        std::sort(ids, ids + 3);
        std::ranges::sort(ids);
        Node node_1 = *std::ranges::find_if(m_Nodes, [&](const Node &node) { return node.ID == id1; });
//        Node node_1 = *std::find_if(std::begin(m_Nodes), std::end(m_Nodes),
//                                    [&](const Node& node) { return node.ID == id1; });
        Node node_2 = *std::find_if(std::begin(m_Nodes), std::end(m_Nodes),
                                    [&](const Node &node) { return node.ID == id2; });
        Node node_3 = *std::find_if(std::begin(m_Nodes), std::end(m_Nodes),
                                    [&](const Node &node) { return node.ID == id3; });
        planes.insert(Plain(node_1, node_2, node_3));
      }
    }
  }

  double min_angel = std::numeric_limits<double>::max();
  for (const auto &plane_1 : planes) {
    auto other = planes | std::ranges::views::filter([&](const Plain &plane) { return plane == plane_1; });
    for (const auto &plane_2 : other) {
      double angel = AngelBetweenPlains(plane_1, plane_2);
      if (angel < min_angel && angel != 0) {
        min_angel = angel;
      }
    }
  }
  constexpr double to_degrees_multiplier = 180.0 / 3.141592653589793238463;
  return min_angel * to_degrees_multiplier;
}

double eps = 0.00000001;

bool find_node_in_line(const std::vector<double> &line_vector,
                       const std::vector<double> &point,
                       const std::vector<double> &point_to_check) {
  double x = (point_to_check[0] - point[0]) / line_vector[0];
  double y = (point_to_check[1] - point[1]) / line_vector[1];
  double z = (point_to_check[2] - point[2]) / line_vector[2];
  bool flag = fabs(x - y) < eps && fabs(y - z) < eps;
  return fabs(x - y) < eps && fabs(y - z) < eps;
}

double Distance(const std::vector<double> &v1, const std::vector<double> &v2) {
  return sqrt(
      (v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1]) + (v1[2] - v2[2]) * (v1[2] - v2[2]));
}

void make_line(auto other,
               std::vector<uint32_t> &index_node_in_line,
               const std::vector<Node> &m_Nodes,
               uint32_t node_id1,
               uint32_t node_id2) {
  std::vector<double> coord_1 = m_Nodes[node_id1 - 1].Coordinates;
  std::vector<double> coord_2 = m_Nodes[node_id2 - 1].Coordinates;
//      std::vector<double> coord_2 = std::ranges::find_if(m_Nodes, [&](const Node &node) { return node.ID == node_id2; })->Coordinates;
  std::vector<double> line_vector{coord_1[0] - coord_2[0],
                                  coord_1[1] - coord_2[1],
                                  coord_1[2] - coord_2[2]};
  std::ranges::for_each(other, [&](const uint32_t &node_idn) {
    std::vector<double> coord_n = m_Nodes[node_idn - 1].Coordinates;
    if (find_node_in_line(line_vector, coord_1, coord_n)) {
      index_node_in_line.push_back(node_idn);
    }
  });
}

//ID: 121 Coordinates: (-0.586747; 0.483391; -0.127803)
//ID: 84 Coordinates: (-0.379073; 0.429080; -0.256691)

//ID: 224 Coordinates: (-0.447556; 0.216329; 0.015566)
//ID: 32 Coordinates: (-0.399999; 0.500000; -0.000001)
//ID: 237 Coordinates: (-0.423778; 0.358164; 0.007783)

//ID: 232 Coordinates: (-0.413315; 0.322705; -0.120563)



std::pair<uint32_t, uint32_t> find_vertex(const std::vector<uint32_t> &index_node_in_line,
                                          const std::vector<Node> &m_Nodes) {
  std::pair<uint32_t, uint32_t> res;
  double max_distance = std::numeric_limits<double>::min();
  for (int i = 0; i < index_node_in_line.size() - 1; ++i) {
    for (int j = i + 1; j < index_node_in_line.size(); ++j) {
      double distance =
          Distance(m_Nodes[index_node_in_line[i] - 1].Coordinates, m_Nodes[index_node_in_line[j] - 1].Coordinates);
      if (distance > max_distance) {
        res = std::make_pair(index_node_in_line[i], index_node_in_line[j]);
        max_distance = distance;
      }
    }
  }
  return res;
}

uint32_t find_D(auto other_3, const std::vector<Node> &m_Nodes, const Plain &plain) {
  uint32_t node_id;
  double max_distance = std::numeric_limits<double>::min();
  std::ranges::for_each(other_3, [&](const uint32_t &ids) {
    double distance = DistanceBetweenPlainAndPoint(plain, m_Nodes[ids - 1].Coordinates);
    if (distance > max_distance) {
      max_distance = distance;
      node_id = ids;
    }
  });
  return node_id;
}

std::vector<uint32_t> MeshLoader::GetFiniteVertexFiniteElement(const FiniteElement &finite_element) const {
  uint32_t arbitrary_point_id = finite_element.NodeIDs.back();
  auto other = finite_element.NodeIDs
      | std::ranges::views::filter([&](const uint32_t &node_id) { return node_id != arbitrary_point_id; });
  std::vector<uint32_t> result;

  for (const uint32_t &node_id1 : other) {
    if (result.empty()) {
      ///first line
      std::vector<uint32_t> index_node_in_line_AB{arbitrary_point_id};
      make_line(other, index_node_in_line_AB, m_Nodes, arbitrary_point_id, node_id1);
      std::pair<uint32_t, uint32_t> points_AB = find_vertex(index_node_in_line_AB, m_Nodes);
      auto other_1 = other | std::ranges::views::filter([&](const uint32_t &node_id) {
        return std::ranges::find(index_node_in_line_AB, node_id) == std::end(index_node_in_line_AB);
      });

      for (const uint32_t &node_id2 : other_1) {
        if (result.empty()) {
          ///second line
          std::vector<uint32_t> index_node_in_line_AC{points_AB.first};
          make_line(other_1, index_node_in_line_AC, m_Nodes, points_AB.first, node_id2);
          std::pair<uint32_t, uint32_t> points_AC = find_vertex(index_node_in_line_AC, m_Nodes);
          auto other_2 = other_1 | std::ranges::views::filter([&](const uint32_t &node_id) {
            return std::ranges::find(index_node_in_line_AC, node_id) == std::end(index_node_in_line_AC);
          });

          if (points_AC.first != points_AB.first)
            break;

          std::vector<uint32_t> index_node_in_line_BC{points_AB.second, points_AC.second};
          make_line(other_2, index_node_in_line_BC, m_Nodes, points_AB.second, points_AC.second);
          std::pair<uint32_t, uint32_t> points_BC = find_vertex(index_node_in_line_BC, m_Nodes);
          auto other_3 = other_2 | std::ranges::views::filter([&](const uint32_t &node_id) {
            return std::ranges::find(index_node_in_line_BC, node_id) == std::end(index_node_in_line_BC);
          });

          if (points_BC.first != points_AB.second)
            continue;

          if (points_BC.second != points_AC.second)
            continue;

          Node A = m_Nodes[points_AB.first - 1];
          Node B = m_Nodes[points_AB.second - 1];
          Node C = m_Nodes[points_AC.second - 1];

          Plain ABC(A, B, C);

          uint32_t D_id = find_D(other_3, m_Nodes, ABC);

          std::vector<uint32_t> index_node_in_line_CD{points_AC.second};
          make_line(other_3, index_node_in_line_CD, m_Nodes, points_AC.second, D_id);

          std::vector<uint32_t> index_node_in_line_BD{points_AB.second};
          make_line(other_2, index_node_in_line_BD, m_Nodes, points_AB.second, D_id);

          std::vector<uint32_t> index_node_in_line_AD{points_AB.second};
          make_line(other_2, index_node_in_line_AD, m_Nodes, points_AB.first, D_id);

          std::vector<uint32_t> test{};
          std::ranges::copy(index_node_in_line_AB, std::back_inserter(test));
          std::ranges::copy(index_node_in_line_AC, std::back_inserter(test));
          std::ranges::copy(index_node_in_line_BC, std::back_inserter(test));
          std::ranges::copy(index_node_in_line_AD, std::back_inserter(test));
          std::ranges::copy(index_node_in_line_BD, std::back_inserter(test));
          std::ranges::copy(index_node_in_line_CD, std::back_inserter(test));

          if (std::set(test.cbegin(), test.cend())
              == std::set(finite_element.NodeIDs.cbegin(), finite_element.NodeIDs.cend())) {
            result.reserve(4);
            result.push_back(points_AB.first);
            result.push_back(points_AB.second);
            result.push_back(points_AC.second);
            result.push_back(D_id);
          }
        }
      }
    }
  }
  return result;
}

double MeshLoader::VolumeFiniteElement(const FiniteElement &finite_element) const {
  std::vector<uint32_t> vertex = GetFiniteVertexFiniteElement(finite_element);
  std::vector<double> A = m_Nodes[vertex[0] - 1].Coordinates;
  std::vector<double> B = m_Nodes[vertex[1] - 1].Coordinates;
  std::vector<double> C = m_Nodes[vertex[2] - 1].Coordinates;
  std::vector<double> D = m_Nodes[vertex[3] - 1].Coordinates;

  std::vector<double> AB{B[0] - A[0], B[1] - A[1], B[2] - A[2]};
  std::vector<double> AC{C[0] - A[0], C[1] - A[1], C[2] - A[2]};
  std::vector<double> AD{D[0] - A[0], D[1] - A[1], D[2] - A[2]};

  return fabs(AB[0] * (AC[1] * AD[2] - AC[2] * AD[1]) - AB[1] * (AC[0] * AD[2] - AC[2] * AD[0])
                  + AB[2] * (AC[0] * AD[1] - AC[1] * AD[0])) / 6;
}

double MeshLoader::SquareSurfaceFiniteElement(const FiniteElement &finite_element) const {
  std::vector<uint32_t> vertex = GetFiniteVertexFiniteElement(finite_element);
  Node A = m_Nodes[vertex[0] - 1];
  Node B = m_Nodes[vertex[1] - 1];
  Node C = m_Nodes[vertex[2] - 1];
  Node D = m_Nodes[vertex[3] - 1];

  return SquarePlain(Plain(A, B, C)) + SquarePlain(Plain(A, B, D)) + SquarePlain(Plain(B, C, D)) + SquarePlain(Plain(A, C, D));
}
