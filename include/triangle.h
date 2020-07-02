#pragma once

#include <vec3.h>

class triangle {
public:
  /**
   * @brief Construct a new triangle object
   *
   * @param a first vertex
   * @param b second vertex
   * @param c third vertex
   *
   * Vertices must be given in counter-clockwise order with respect to the
   * right-hand rule.
   */
  triangle(const vec3 &a, const vec3 &b, const vec3 &c) : vert_{a, b, c} {}

  inline vec3 normal(bool normalize = true) const {
    vec3 normal = (vert_[2] - vert_[1]) * (vert_[0] - vert_[1]);
    if (normalize) {
      normal.normalize();
    }
    return normal;
  }

  /**
   * @brief returns the vertice at index pos
   *
   * @throw if pos is not in {0, 1, 2}
   *
   * @param pos
   * @return const vec3&
   */
  inline const vec3 &at(int pos) const {
    if (pos < 0 || pos > 2) {
      throw std::out_of_range("pos " + std::to_string(pos) +
                              " is out of range (max 2)");
    }
    return vert_[pos];
  }

private:
  vec3 vert_[3];
};