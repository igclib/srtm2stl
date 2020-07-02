#pragma once

#include <cmath>

class vec3 {
public:
  vec3() : e_{0, 0, 0} {};
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

  template <typename T> inline void array(T e[3]) const {
    e[0] = e_[0];
    e[1] = e_[1];
    e[2] = e_[2];
  };

  inline double x() const { return e_[0]; }
  inline double y() const { return e_[1]; }
  inline double z() const { return e_[2]; }

private:
  double e_[3];
};