#pragma once

#include <cmath>

/**
 * @brief A basic 3D vector
 *
 */
class vec3 {
public:
  /**
   * @brief Construct a new vec3 with all coordinates at zero
   *
   */
  vec3() : e_{0, 0, 0} {};
  /**
   * @brief Construct a new vec3 with coordinates a, b and c
   *
   * @param a x-coordinate
   * @param b y-coordinate
   * @param c z-coordinate
   */
  vec3(double a, double b, double c) : e_{a, b, c} {}

  /**
   * @brief assign the vector to other
   */
  inline vec3 &operator=(const vec3 &other) {
    e_[0] = other.e_[0];
    e_[1] = other.e_[1];
    e_[2] = other.e_[2];
    return *this;
  }
  /**
   * @brief returns the cross product of the vector with other
   */
  inline vec3 operator*(const vec3 &other) const {
    return vec3(e_[1] * other.e_[2] - e_[2] * other.e_[1],
                e_[2] * other.e_[0] - e_[0] * other.e_[2],
                e_[0] * other.e_[1] - e_[1] * other.e_[0]);
  }

  /**
   * @brief scales the vector by alpha
   */
  inline vec3 operator*(double alpha) const {
    return vec3(alpha * e_[0], alpha * e_[1], alpha * e_[2]);
  }

  /**
   * @brief adds other to the vector
   */
  inline vec3 operator+(const vec3 &other) const {
    return vec3(e_[0] + other.e_[0], e_[1] + other.e_[1], e_[2] + other.e_[2]);
  }

  /**
   * @brief subtracts other from the vector
   */
  inline vec3 operator-(const vec3 &other) const {
    return vec3(e_[0] - other.e_[0], e_[1] - other.e_[1], e_[2] - other.e_[2]);
  }

  /**
   * @brief returns the squared norm of the vector
   */
  inline double sqnorm() const {
    return e_[0] * e_[0] + e_[1] * e_[1] + e_[2] * e_[2];
  }

  /**
   * @brief returns the norm of the vector
   */
  inline double norm() const { return std::sqrt(sqnorm()); }

  /**
   * @brief normalizes the vector to unit length
   */
  inline vec3 &normalize() {
    double n = norm();
    e_[0] /= n;
    e_[1] /= n;
    e_[2] /= n;
    return *this;
  }

  /**
   * @brief populates the given array with the vector coordinates
   */
  template <typename T> inline void array(T e[3]) const {
    e[0] = e_[0];
    e[1] = e_[1];
    e[2] = e_[2];
  };

  /**
   * @brief returns the x coordinate
   */
  inline double x() const { return e_[0]; }

  /**
   * @brief returns the y coordinate
   */
  inline double y() const { return e_[1]; }

  /**
   * @brief returns the z coordinate
   */
  inline double z() const { return e_[2]; }

private:
  double e_[3];
};