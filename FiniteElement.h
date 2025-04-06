#ifndef FINITE_ELEMENT_GRID_LOADER__FINITEELEMENT_H_
#define FINITE_ELEMENT_GRID_LOADER__FINITEELEMENT_H_
#include <cstdint>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>


struct FiniteElement {
  uint32_t ID = 0;
  uint32_t  PositionID = 0;
  std::vector<uint32_t> NodeIDs;

  FiniteElement(const std::string& str, uint32_t size, uint32_t id) : ID(id) {
    std::stringstream s{str};
    std::string number_str;
    NodeIDs.reserve(size);

    s >> number_str;
    PositionID = std::stoi(number_str);
    number_str.clear();

    for (int i = 0; i < size; ++i, number_str.clear()) {
      s >> number_str;
      NodeIDs.emplace_back(std::stoi(number_str));
    }
  }

  friend std::ostream& operator<<(std::ostream &os, const FiniteElement& finite_element) {
    os << "ID: " << finite_element.ID << " PositionID: " << finite_element.PositionID << " NodeIDs: [";
    for (int i = 0; i < finite_element.NodeIDs.capacity() - 1; i++) {
      os << finite_element.NodeIDs[i] << "; ";
    }
    os << finite_element.NodeIDs[finite_element.NodeIDs.capacity() - 1];
    os << ']';
    return os;
  }

};

#endif //FINITE_ELEMENT_GRID_LOADER__FINITEELEMENT_H_
