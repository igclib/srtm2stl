#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

namespace srtm {
// constants
const int SRTM_SIZE = 1201;
const size_t TILENAME_SIZE = sizeof("N00E000.hgt");

// typedefs
using heightmap_t = short[SRTM_SIZE][SRTM_SIZE];
using tilename_t = char[TILENAME_SIZE];
using vec3_t = short[3];
using triangle_t = vec3_t[4]; // 3 vertices + 1 normal
using trianglelist_t = std::vector<triangle_t>;
using interpolation_t = enum {
  SWNE_DIAGONAL,
  NWSE_DIAGONAL,
  ALTERNATING_DIAGONALS,
  SHORTEST_DIAGONAL,
  DELAUNAY
};

void get_tile_name(double lat, double lon, tilename_t &tile_name) {
  char latchar = lat > 0.0 ? 'N' : 'S';
  char lonchar = lon > 0.0 ? 'E' : 'W';
  double intlat, intlon;
  modf(lat, &intlat);
  modf(lon, &intlon);
  snprintf(tile_name, TILENAME_SIZE, "%c%02d%c%03d.hgt", latchar, (int)intlat,
           lonchar, (int)intlon);
}

void get_tile_pixel(double lat, double lon, short &xtile, short &ytile) {
  double intlat, intlon;
  double latdec = modf(lat, &intlat);
  double londec = modf(lon, &intlon);
  xtile = std::round(SRTM_SIZE - 1 - SRTM_SIZE * latdec);
  ytile = std::round(SRTM_SIZE * londec);
}

void read_tile(const std::string &filename, heightmap_t &heightmap) {
  std::ifstream file(filename, std::ios::in | std::ios::binary);
  if (!file) {
    throw std::runtime_error("Could not open " + filename);
  }

  unsigned char buffer[2];
  short alt;
  for (int i = 0; i < SRTM_SIZE; ++i) {
    for (int j = 0; j < SRTM_SIZE; ++j) {
      if (!file.read(reinterpret_cast<char *>(buffer), sizeof(buffer))) {
        std::string error =
            "Error reading file (" + std::to_string(file.tellg()) + ")";
        throw std::runtime_error(error);
      }
      // swap bytes as they are stored in big endian
      alt = (buffer[0] << 8) | buffer[1];
      heightmap[i][j] = alt;
    }
  }
}

short get_elevation(double lat, double lon, const std::string &prefix) {
  tilename_t tilename;
  get_tile_name(lat, lon, tilename);
  short xtile, ytile;
  get_tile_pixel(lat, lon, xtile, ytile);
  std::string filename = prefix + tilename;
  heightmap_t heightmap;
  read_tile(filename, heightmap);
  return heightmap[xtile][ytile];
}

trianglelist_t triangulate(const heightmap_t &heightmap,
                           interpolation_t interpolation) {
  trianglelist_t triangles;
  return triangles;
}

} // namespace srtm

int main(int argc, const char *argv[]) {
  std::cout << srtm::get_elevation(44.684413, 6.614422, "test/") << std::endl;

  return EXIT_SUCCESS;
}