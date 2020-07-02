#pragma once

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits.h>
#include <numeric>
#include <string>
#include <vector>

#include <triangle.h>
#include <vec3.h>

/**
 * @brief The hgt class encapsulates elevation models of SRTM files.
 */
class hgt {

  static const int SRTM_1_SIZE = 3601;
  static const int SRTM_3_SIZE = 1201;
  static const int SRTM_1_FILE_SIZE = SRTM_1_SIZE * SRTM_1_SIZE * 2;
  static const int SRTM_3_FILE_SIZE = SRTM_3_SIZE * SRTM_3_SIZE * 2;
  static const short SRTM_NO_DATA = -32768;
  static const int TILENAME_SIZE = sizeof("N00E000.hgt");

  using interpolation_t = enum {
    SWNE_DIAGONAL,
    NWSE_DIAGONAL,
    ALTERNATING_DIAGONALS,
    SHORTEST_DIAGONAL,
    DELAUNAY,
    DEFAULT = SWNE_DIAGONAL,
  };

public:
  /**
   * @brief Construct a new hgt object from a complete tile
   *
   * @param filename path to a .hgt tile
   */
  hgt(const std::string &filename) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file) {
      throw std::runtime_error("Could not open " + filename);
    }

    /**
     * @todo handle SRTM-1
     */
    auto file_size = std::filesystem::file_size(filename);
    if (file_size == SRTM_1_FILE_SIZE) {
      throw std::runtime_error("SRTM-1 is not yet supported");
    } else if (file_size != SRTM_3_FILE_SIZE) {
      throw std::runtime_error("File size " + std::to_string(file_size) +
                               " does not match SRTM-3 file size");
    }

    unsigned char buffer[2];
    short alt;
    for (int i = 0; i < SRTM_3_SIZE; ++i) {
      for (int j = 0; j < SRTM_3_SIZE; ++j) {
        if (!file.read(reinterpret_cast<char *>(buffer), sizeof(buffer))) {
          std::string error = "Error reading " + filename + " (byte " +
                              std::to_string(file.tellg()) + ")";
          throw std::runtime_error(error);
        }
        // swap bytes as they are stored in big endian
        alt = (buffer[0] << 8) | buffer[1];
        heightmap_[i][j] = alt;
      }
    }

    // fix missing data
    std::vector<short> nearest;
    int delta;
    for (int i = 0; i < SRTM_3_SIZE; ++i) {
      for (int j = 0; j < SRTM_3_SIZE; ++j) {
        if (heightmap_[i][j] == SRTM_NO_DATA) {
          // nearest north
          delta = 0;
          while ((i - delta) > 0 && heightmap_[i - delta][j] == SRTM_NO_DATA) {
            ++delta;
          }
          if (heightmap_[i - delta][j] != SRTM_NO_DATA) {
            nearest.push_back(heightmap_[i - delta][j]);
          }
          // nearest south
          delta = 0;
          while ((i + delta) < SRTM_3_SIZE - 1 &&
                 heightmap_[i + delta][j] == SRTM_NO_DATA) {
            ++delta;
          }
          if (heightmap_[i + delta][j] != SRTM_NO_DATA) {
            nearest.push_back(heightmap_[i + delta][j]);
          }
          // nearest west
          delta = 0;
          while ((j - delta) > 0 && heightmap_[i][j - delta] == SRTM_NO_DATA) {
            ++delta;
          }
          if (heightmap_[i][j - delta] != SRTM_NO_DATA) {
            nearest.push_back(heightmap_[i][j - delta]);
          }
          // nearest east
          delta = 0;
          while ((j + delta) < SRTM_3_SIZE - 1 &&
                 heightmap_[i][j + delta] == SRTM_NO_DATA) {
            ++delta;
          }
          if (heightmap_[i][j + delta] != SRTM_NO_DATA) {
            nearest.push_back(heightmap_[i][j + delta]);
          }

          if (nearest.empty()) {
            std::cerr << "no data fill" << std::endl;
          } else {
            heightmap_[i][j] =
                std::accumulate(nearest.begin(), nearest.end(), 0.0) /
                nearest.size();
            nearest.clear();
          }
        }
      }
    }

    std::string basename = std::filesystem::path(filename).filename();
    bool north = basename.substr(0, 1) == "N";
    bool east = basename.substr(3, 1) == "E";
    double lat = std::stod(basename.substr(1, 2));
    double lon = std::stod(basename.substr(4, 3));

    swlat_ = north ? lat : -lat;
    swlon_ = east ? lon : -lon;
  }

  void pixel(double lat, double lon, short &xtile, short &ytile) const {
    // TODO check if it is the right tile
    double intlat, intlon;
    double latdec = modf(lat, &intlat);
    double londec = modf(lon, &intlon);
    xtile = std::round(SRTM_3_SIZE - 1 - SRTM_3_SIZE * latdec);
    ytile = std::round(SRTM_3_SIZE * londec);
  }

  inline short at(double lat, double lon) const {
    short x, y;
    pixel(lat, lon, x, y);
    return heightmap_[x][y];
  }

  std::vector<triangle>
  triangulate(interpolation_t interpolation = DEFAULT) const {
    std::vector<triangle> triangles;
    int stepsize = 90; // meters
    double ax, ay, az;
    double bx, by, bz;
    double cx, cy, cz;
    double dx, dy, dz;
    vec3 a, b, c, d;

    /*
     *    A - B
     *    |   |
     *    C - D
     */

    for (int i = 0; i < SRTM_3_SIZE - 1; ++i) {
      for (int j = 0; j < SRTM_3_SIZE - 1; ++j) {
        ax = j * stepsize;
        ay = i * stepsize;
        az = heightmap_[i][j];
        a = vec3(ax, ay, az);

        bx = ax + stepsize;
        // bx = ax;
        by = ay;
        bz = heightmap_[i][j + 1];
        b = vec3(bx, by, bz);

        cx = ax;
        cy = ay + stepsize;
        cz = heightmap_[i + 1][j];
        c = vec3(cx, cy, cz);

        dx = ax + stepsize;
        dy = ay + stepsize;
        dz = heightmap_[i + 1][j + 1];
        d = vec3(dx, dy, dz);

        int max = 4000;
        if (az > max || bz > max || cz > max || dz > max) {
          std::cerr << "outps" << std::endl;
        }

        switch (interpolation) {
        case interpolation_t::SWNE_DIAGONAL: {
          /*
           *    A - B
           *    | / |
           *    C - D
           */
          triangles.emplace_back(a, c, b);
          triangles.emplace_back(b, c, d);
          break;
        }
        case interpolation_t::NWSE_DIAGONAL: {
          /*
           *    A - B
           *    | \ |
           *    C - D
           */
          triangles.emplace_back(a, c, d);
          triangles.emplace_back(d, b, a);
          break;
        }
        default:
          break;
        }
      }
    }

    return triangles;
  }

  void toASCII(const std::string &outfile) const {
    std::ofstream f(outfile);
    vec3 vec;
    f << "solid tile" << std::endl;
    for (const triangle &tri : triangulate()) {
      vec = tri.normal(true);
      f << "  facet normal " << vec.x() << " " << vec.y() << " " << vec.z()
        << std::endl;
      f << "    outer loop" << std::endl;
      for (int i = 0; i < 3; ++i) {
        f << "      vertex " << std::scientific << tri.at(i).x() << " "
          << tri.at(i).y() << " " << tri.at(i).z() << std::endl;
      }
      f << "    endloop" << std::endl;
      f << "  endfacet" << std::endl;
    }
    f << "endsolid tile" << std::endl;
  }

  void toSTL(const std::string &outfile) {
    static_assert(sizeof(float) * CHAR_BIT == 32, "require 32 bits floats");
    std::vector<triangle> tri = triangulate();

    std::ofstream f(outfile, std::ios::binary);

    uint8_t header[80] = {};
    uint32_t n_tri = tri.size();
    // uint32_t n_tri = 100;
    float arr[3] = {};
    char attribute[2] = {};

    f.write(reinterpret_cast<char *>(header), sizeof(header));
    f.write(reinterpret_cast<char *>(&n_tri), sizeof(n_tri));

    for (uint32_t i = 0; i < n_tri; ++i) {
      tri.at(i).normal(true).array(arr);
      f.write(reinterpret_cast<char *>(arr), sizeof(arr)); // 3*4 bytes = 12

      tri.at(i).at(0).array(arr);
      f.write(reinterpret_cast<char *>(arr), sizeof(arr)); // 3*4 bytes = 12

      tri.at(i).at(1).array(arr);
      f.write(reinterpret_cast<char *>(arr), sizeof(arr)); // 3*4 bytes = 12

      tri.at(i).at(2).array(arr);
      f.write(reinterpret_cast<char *>(arr), sizeof(arr)); // 3*4 bytes = 12

      f.write(reinterpret_cast<char *>(attribute),
              sizeof(attribute)); // 2 bytes
    }

    std::cerr << n_tri << std::endl;
    std::cerr << (sizeof(header) + sizeof(n_tri) + n_tri * 50) << std::endl;
    std::cerr << f.tellp() << std::endl;
  }

  static std::string tilename(double lat, double lon) {
    char tilename[TILENAME_SIZE];
    char latchar = lat > 0.0 ? 'N' : 'S';
    char lonchar = lon > 0.0 ? 'E' : 'W';
    double intlat, intlon;
    modf(lat, &intlat);
    modf(lon, &intlon);
    snprintf(tilename, TILENAME_SIZE, "%c%02d%c%03d.hgt", latchar, (int)intlat,
             lonchar, (int)intlon);
    /*
    once C++20 formatting features are implemented by g++
    return std::format("{}{:0>3}{}{:0>2}.hgt", latchar, intlat, lonchar,
                       intlon);
    */
    return std::string(tilename);
  }

private:
  short swlat_;
  short swlon_;
  short heightmap_[SRTM_3_SIZE][SRTM_3_SIZE];
};