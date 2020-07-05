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

  static const short TWO_BYTES = 2;
  static const int SRTM_1_ROWS = 3601;
  static const int SRTM_1_COLUMNS = 3601;
  static const int SRTM_3_ROWS = 1201;
  static const int SRTM_3_COLUMNS = 1201;
  static const int SRTM_1_FILE_SIZE = SRTM_1_ROWS * SRTM_1_COLUMNS * TWO_BYTES;
  static const int SRTM_3_FILE_SIZE = SRTM_3_ROWS * SRTM_3_COLUMNS * TWO_BYTES;
  static const short SRTM_NO_DATA = -32768;
  static const int TILENAME_SIZE = sizeof("N00E000.hgt");

  enum class interpolation {
    SWNE_DIAGONAL,
    NWSE_DIAGONAL,
    // ALTERNATING_DIAGONALS,
    // SHORTEST_DIAGONAL,
    // DELAUNAY,
    DEFAULT = SWNE_DIAGONAL,
  };

  enum class fix_missing {
    NONE,
    NEAREST,
    DEFAULT = NEAREST,
  };

public:
  /**
   * @brief Construct a new hgt object from a complete tile
   *
   * @param filename path to a HGT tile
   * @param strategy what to do with missing data points
   */
  hgt(const std::string &filename,
      fix_missing strategy = fix_missing::DEFAULT) {
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

    unsigned char buffer[TWO_BYTES];
    short alt;
    for (int i = 0; i < SRTM_3_ROWS; ++i) {
      for (int j = 0; j < SRTM_3_COLUMNS; ++j) {
        if (!file.read(reinterpret_cast<char *>(buffer), sizeof(buffer))) {
          std::string error = "Error reading " + filename + " (byte " +
                              std::to_string(file.tellg()) + ")";
          throw std::runtime_error(error);
        }
        // bytes are stored in big endian
        alt = (buffer[0] << 8) | buffer[1];
        heightmap_[i][j] = alt;
      }
    }

    // fix missing data
    fix_missing_data(strategy);

    // store tile misc information
    std::string basename = std::filesystem::path(filename).filename();
    bool north = basename.substr(0, 1) == "N";
    bool east = basename.substr(3, 1) == "E";
    double lat = std::stod(basename.substr(1, 2));
    double lon = std::stod(basename.substr(4, 3));

    swlat_ = north ? lat : -lat;
    swlon_ = east ? lon : -lon;
  }

  void fix_missing_data(fix_missing strategy) {
    switch (strategy) {
    case fix_missing::NEAREST: {
      std::vector<short> nearest;
      int delta;
      for (int i = 0; i < SRTM_3_ROWS; ++i) {
        for (int j = 0; j < SRTM_3_ROWS; ++j) {
          if (heightmap_[i][j] == SRTM_NO_DATA) {
            // nearest north
            delta = 0;
            while ((i - delta) > 0 &&
                   heightmap_[i - delta][j] == SRTM_NO_DATA) {
              ++delta;
            }
            if (heightmap_[i - delta][j] != SRTM_NO_DATA) {
              nearest.push_back(heightmap_[i - delta][j]);
            }
            // nearest south
            delta = 0;
            while ((i + delta) < SRTM_3_ROWS - 1 &&
                   heightmap_[i + delta][j] == SRTM_NO_DATA) {
              ++delta;
            }
            if (heightmap_[i + delta][j] != SRTM_NO_DATA) {
              nearest.push_back(heightmap_[i + delta][j]);
            }
            // nearest west
            delta = 0;
            while ((j - delta) > 0 &&
                   heightmap_[i][j - delta] == SRTM_NO_DATA) {
              ++delta;
            }
            if (heightmap_[i][j - delta] != SRTM_NO_DATA) {
              nearest.push_back(heightmap_[i][j - delta]);
            }
            // nearest east
            delta = 0;
            while ((j + delta) < SRTM_3_ROWS - 1 &&
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
    }
    }
  }

  /**
   * @brief computes the tile's xy-coordinates corresponding to the given
   * geographic point
   *
   * @param lat
   * @param lon
   * @param xtile
   * @param ytile
   *
   * @throw std::runtime_error if the point does not belong to the tile
   */
  inline void pixel(double lat, double lon, short &xtile, short &ytile) const {
    double intlat, intlon;
    double latdec = modf(lat, &intlat);
    double londec = modf(lon, &intlon);
    if (intlat != swlat_ || intlon != swlon_) {
      throw std::runtime_error("not the right tile");
    }
    xtile = std::round(SRTM_3_ROWS - 1 - SRTM_3_ROWS * latdec);
    ytile = std::round(SRTM_3_ROWS * londec);
  }

  /**
   * @brief computes the altitude of the tile at the given geographic point
   *
   * @param lat geographic latitude
   * @param lon geographic longitude
   * @return short altitude of the point
   */
  inline short at(double lat, double lon) const {
    short x, y;
    pixel(lat, lon, x, y);
    return heightmap_[x][y];
  }

  /**
   * @brief computes a triangulation of the heightmap
   *
   * @param strategy is the interpolation strategy
   * @return std::vector<triangle>
   */
  std::vector<triangle>
  triangulate(interpolation strategy = interpolation::DEFAULT) const {
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

    for (int i = 0; i < SRTM_3_ROWS - 1; ++i) {
      for (int j = 0; j < SRTM_3_ROWS - 1; ++j) {
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

        switch (strategy) {
        case interpolation::SWNE_DIAGONAL: {
          /*
           *    A - B
           *    | / |
           *    C - D
           */
          triangles.emplace_back(a, c, b);
          triangles.emplace_back(b, c, d);
          break;
        }
        case interpolation::NWSE_DIAGONAL: {
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
          throw ::std::runtime_error("Interpolation not implemented");
          break;
        }
      }
    }

    return triangles;
  }

  /**
   * @brief generates an STL mesh from the heightmap (ASCII version)
   *
   * @param outfile path to the output file
   */
  void toASCII(const std::string &outfile) const {
    std::ofstream f(outfile);
    vec3 vec;
    f << "solid tile" << std::endl;
    for (const triangle &tri : triangulate()) {
      vec = tri.normal(true);
      f << "\tfacet normal " << vec.x() << " " << vec.y() << " " << vec.z()
        << "\n";
      f << "\t\touter loop\n";
      for (int i = 0; i < 3; ++i) {
        f << "\t\t\tvertex " << std::scientific << tri.at(i).x() << " "
          << tri.at(i).y() << " " << tri.at(i).z() << "\n";
      }
      f << "\t\tendloop\n\tendfacet\n";
    }
    f << "endsolid tile\n";
  }

  /**
   * @brief generates an STL mesh from the heightmap (binary version)
   *
   * @param outfile path to the output file
   */
  void toSTL(const std::string &outfile) {
    static_assert(sizeof(float) * CHAR_BIT == 32, "require 32 bits floats");
    std::vector<triangle> tri = triangulate();

    std::ofstream f(outfile, std::ios::binary);

    uint8_t header[80] = {};
    uint32_t n_tri = tri.size();
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
  }

  /**
   * @brief get the tile filename corresponding to a geographic point
   *
   * @param lat latitude of the point
   * @param lon longitude of the point
   * @return std::string filename of the tile
   */
  static std::string tilename(double lat, double lon) {
    char tilename[TILENAME_SIZE];
    char latchar = lat > 0.0 ? 'N' : 'S';
    char lonchar = lon > 0.0 ? 'E' : 'W';
    double intlat, intlon;
    modf(lat, &intlat);
    modf(lon, &intlon);
    snprintf(tilename, TILENAME_SIZE, "%c%02d%c%03d.hgt", latchar, int(intlat),
             lonchar, int(intlon));
    /*
    once C++20 formatting features are implemented by g++
    return std::format("{}{:0>3}{}{:0>2}.hgt", latchar, intlat, lonchar,
                       intlon);
    */
    return std::string(tilename);
  }

  /**
   * @brief computes the elevation of a single geographic point
   *
   * @param lat latitude of the point
   * @param lon longitude of the point
   * @param prefix path to a directory containing HGT files
   * @return double AMSL altitude of the point
   */
  static double elevation(double lat, double lon, const std::string &prefix) {
    std::string tilename = hgt::tilename(lat, lon);
    std::string filename = prefix + tilename;
    hgt heightmap(filename);
    return heightmap.at(lat, lon);
  }

  /**
   * @brief computes the elevation of a collection of geographic points
   *
   * @param latlon vector of lat lon pairs
   * @param prefix path to a directory containing HGT files
   *
   * @todo this is quite inefficient. maybe
   */
  static std::vector<double>
  elevations(std::vector<std::pair<double, double>> latlon,
             const std::string &prefix) {
    std::vector<double> altitudes;
    altitudes.reserve(latlon.size());
    for (const auto &ll : latlon) {
      altitudes.push_back(elevation(ll.first, ll.second, prefix));
    }
    return altitudes;
  }

private:
  short swlat_;
  short swlon_;
  short heightmap_[SRTM_3_ROWS][SRTM_3_ROWS];
};