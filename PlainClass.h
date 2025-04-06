#ifndef FINITE_ELEMENT_GRID_LOADER__PLAINCLASS_H_
#define FINITE_ELEMENT_GRID_LOADER__PLAINCLASS_H_

struct Plain {
  Node node_1;
  Node node_2;
  Node node_3;

  Plain(const Node &node1, const Node &node2, const Node &node3) : node_1(node1), node_2(node2), node_3(node3) {};

  bool operator==(const Plain &other) const {
    return std::tie(node_1.ID, node_2.ID, node_3.ID) == std::tie(other.node_1.ID, other.node_2.ID, other.node_3.ID);
  }
};

double VectorLength(std::vector<double> &v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] + v[2]);
}

double ScalarProduct(const std::vector<double> &v1, const std::vector<double> &v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

std::vector<double> NormalVector(const Plain &plain) {
  std::vector<double> A = plain.node_1.Coordinates;
  std::vector<double> B = plain.node_2.Coordinates;
  std::vector<double> C = plain.node_3.Coordinates;
  std::vector<double> AB{B[0] - A[0], B[1] - A[1], B[2] - A[2]};
  std::vector<double> AC{C[0] - A[0], C[1] - A[1], C[2] - A[2]};
  return {AB[1] * AC[2] - AB[2] * AC[1], -(AB[0] * AC[2] - AB[2] * AC[0]), AB[0] * AC[1] - AB[1] * AC[0]};
}

std::vector<double> EquationPlain(const Plain &plain) {
  std::vector<double> A = plain.node_1.Coordinates;
  std::vector<double> B = plain.node_2.Coordinates;
  std::vector<double> C = plain.node_3.Coordinates;
  std::vector<double> AB{B[0] - A[0], B[1] - A[1], B[2] - A[2]};
  std::vector<double> AC{C[0] - A[0], C[1] - A[1], C[2] - A[2]};
  std::vector<double>
      normal_vector{AB[1] * AC[2] - AB[2] * AC[1], -(AB[0] * AC[2] - AB[2] * AC[0]), AB[0] * AC[1] - AB[1] * AC[0]};
  double D = -(A[0] * normal_vector[0] + A[1] * normal_vector[1] + A[2] * normal_vector[2]);

  normal_vector.push_back(D);
  return normal_vector;
}

double AngelBetweenPlains(const Plain &plain1, const Plain &plain2) {
  std::vector<double> n1(NormalVector(plain1));
  std::vector<double> n2(NormalVector(plain2));
  return acos(ScalarProduct(n1, n2) / VectorLength(n1) / VectorLength(n2));
}

double DistanceBetweenPlainAndPoint(const Plain &plain, const std::vector<double> &point) {
  const std::vector<double> &equation_plain = EquationPlain(plain);
  return fabs(
      equation_plain[0] * point[0] + equation_plain[1] * point[1] + equation_plain[2] * point[2] + equation_plain[3]) /
      sqrt(equation_plain[0] * equation_plain[0] + equation_plain[1] * equation_plain[1]
               + equation_plain[2] * equation_plain[2]);
}

double SquarePlain(const Plain &plain) {
  std::vector<double> det = NormalVector(plain);
  return sqrt(det[0] * det[0] + det[1] * det[1] + det[2] * det[2]) / 2;
}


#endif //FINITE_ELEMENT_GRID_LOADER__PLAINCLASS_H_
