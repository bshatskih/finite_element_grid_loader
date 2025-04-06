#ifndef FINITE_ELEMENT_GRID_LOADER_NODE_H_
#define FINITE_ELEMENT_GRID_LOADER_NODE_H_
#include <cstdint>
#include <sstream>
#include <ostream>
#include <vector>
#include <string>


struct Node {
  uint32_t ID = 0;
  std::vector<double> Coordinates;

  Node(const std::string& str, uint32_t dimension, uint32_t id) : ID(id) {
    std::stringstream s{str};
    std::string number_str;
    Coordinates.reserve(dimension);

    for (int i = 0; i < dimension; ++i, number_str.clear()) {
      s >> number_str;
      Coordinates.emplace_back(std::stod(number_str));
    }
  }

  Node(const Node& node_1, const Node& node_2, uint32_t dimension, uint32_t id) : ID(id) {
    Coordinates.reserve(dimension);
    for (int i = 0; i < dimension; i++) {
      Coordinates.emplace_back((node_1.Coordinates[i] + node_2.Coordinates[i]) / 2);
    }
  }

  friend std::ostream& operator<<(std::ostream &os, const Node& node) {
    os << "ID: " << node.ID << " Coordinates: (";
    for (int i = 0; i < node.Coordinates.capacity() - 1; i++) {
      os << node.Coordinates[i] << "; ";
    }
    os << node.Coordinates[node.Coordinates.capacity() - 1];
    os << ')';
    return os;
  }

//  bool operator==(const Node& other) const {
//    return ID == other.ID;
//  }
};


#endif //FINITE_ELEMENT_GRID_LOADER_NODE_H_
