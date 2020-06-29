#pragma once

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <triangle.h>

class hgt {

  static const int SRTM_SIZE = 1201;
  static const int TILENAME_SIZE = sizeof("N00E000.hgt");

  using interpolation_t = enum {
    SWNE_DIAGONAL,
    NWSE_DIAGONAL,
    ALTERNATING_DIAGONALS,
    SHORTEST_DIAGONAL,
    DELAUNAY
  };

public:
  hgt(const std::string &filename) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file) {
      throw std::runtime_error("Could not open " + filename);
    }

    unsigned char buffer[2];
    short alt;
    for (int i = 0; i < SRTM_SIZE; ++i) {
      for (int j = 0; j < SRTM_SIZE; ++j) {
        if (!file.read(reinterpret_cast<char *>(buffer), sizeof(buffer))) {
          std::string error = "Error reading " + filename + " (byte" +
                              std::to_string(file.tellg()) + ")";
          throw std::runtime_error(error);
        }
        // swap bytes as they are stored in big endian
        alt = (buffer[0] << 8) | buffer[1];
        heightmap_[i][j] = alt;
      }
    }

    bool north = filename.substr(0, 1) == "N";
    bool east = filename.substr(3, 1) == "E";
    double lat = std::stod(filename.substr(1, 2));
    double lon = std::stod(filename.substr(4, 3));

    swlat_ = north ? lat : -lat;
    swlon_ = east ? lon : -lon;
  }

  void pixel(double lat, double lon, short &xtile, short &ytile) const {
    // TODO check if it is the right tile
    double intlat, intlon;
    double latdec = modf(lat, &intlat);
    double londec = modf(lon, &intlon);
    xtile = std::round(SRTM_SIZE - 1 - SRTM_SIZE * latdec);
    ytile = std::round(SRTM_SIZE * londec);
  }

  inline short at(double lat, double lon) const {
    short x, y;
    pixel(lat, lon, x, y);
    return heightmap_[x][y];
  }

  std::vector<triangle> triangulate(hgt::interpolation_t interpolation) {
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

    for (int i = 0; i < SRTM_SIZE - 1; ++i) {
      for (int j = 0; j < SRTM_SIZE - 1; ++j) {
        ax = i * stepsize;
        ay = j * stepsize;
        az = heightmap_[i][j];
        a = vec3(ax, ay, az);

        bx = ax + stepsize;
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

        switch (interpolation) {
        case interpolation_t::SWNE_DIAGONAL:
          /*
           *    A - B
           *    | / |
           *    C - D
           */
          triangles.emplace_back(a, c, b);
          triangles.emplace_back(b, c, d);
          break;

        default:
          break;
        }
      }
    }

    return triangles;
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
    return std::string(tilename);
  }

private:
  short swlat_;
  short swlon_;
  short heightmap_[SRTM_SIZE][SRTM_SIZE];
};