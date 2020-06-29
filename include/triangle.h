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

private:
  vec3 vert_[3];
};