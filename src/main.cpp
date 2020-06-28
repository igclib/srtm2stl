#include <cmath>
#include <fstream>
#include <iostream>

const int SRTM_SIZE = 1201;
const size_t TILENAME_SIZE = sizeof("N00E000.hgt");

using heightmap_t = short[SRTM_SIZE][SRTM_SIZE];
using tilename_t = char[TILENAME_SIZE];
using tilepixel_t = std::pair<short, short>;

void get_tile_name(double lat, double lon, tilename_t &tile_name) {
  char latchar = lat > 0.0 ? 'N' : 'S';
  char lonchar = lon > 0.0 ? 'E' : 'W';
  double intlat, intlon;
  modf(lat, &intlat);
  modf(lon, &intlon);
  snprintf(tile_name, TILENAME_SIZE, "%c%02d%c%03d.hgt", latchar, (int)intlat,
           lonchar, (int)intlon);
}

tilepixel_t get_tile_pixel(double lat, double lon) {
  double intlat, intlon;
  double latdec = modf(lat, &intlat);
  double londec = modf(lon, &intlon);
  short x = std::round(SRTM_SIZE - 1 - SRTM_SIZE * latdec);
  short y = std::round(SRTM_SIZE * londec);
  return std::make_pair(x, y);
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
      alt = (buffer[0] << 8) | buffer[1]; // big endian order
      heightmap[i][j] = alt;
    }
  }
}

int get_elevation(double lat, double lon, const std::string &prefix) {
  tilename_t tile;
  get_tile_name(lat, lon, tile);
  tilepixel_t pixel = get_tile_pixel(lat, lon);
  heightmap_t heightmap;
  std::string filename = prefix + tile;
  read_tile(filename, heightmap);
  return heightmap[pixel.first][pixel.second];
}

int main(int argc, const char *argv[]) {
  std::cout << get_elevation(44.684413, 6.614422, "test/") << std::endl;

  return EXIT_SUCCESS;
}