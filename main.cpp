#include <iostream>
#include "AneuMeshLoader.h"

int main() {
  AneuMeshLoader test;
  test.LoadMesh("../two_cylinders.anue");
  std::cout << test << '\n';
//  std::vector<uint32_t> res(test.FindFiniteElementsByEdge(232, 232));
//  for (auto i: res) {
//    std::cout << i << ' ';
//  }
//  auto res = test.GetAllNeighbors();
//  for (const auto &i: res) {
//    std::cout << i.first << ": ";
//    for (const uint32_t& j: i.second) {
//      std::cout << j <<  ' ';
//    }
//    std::cout << '\n';
//  }
  //18: 17 19 167 139 191
  //299 422 540

//    std::vector<uint32_t> res_1(test.FindFiniteElementsByEdge(18, 18));
//  for (auto i: res_1) {
//    std::cout << i << ' ';
//  }
  test.BisectAllEdges();
  std::cout << "\n\n" << test << "\n\n";

  FiniteElement finite_element = test.GetFiniteElements()[10];
//  std::cout << finite_element;
  std::vector<uint32_t> res = test.GetFiniteVertexFiniteElement(finite_element);

  std::ranges::copy(res, std::ostream_iterator<int>(std::cout, " "));
//ID: 84 Coordinates: (-0.379073; 0.429080; -0.256691)
//ID: 224 Coordinates: (-0.447556; 0.216329; 0.015566)
  return 0;
}
