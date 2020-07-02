#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <hgt.h>
#include <triangle.h>
#include <vec3.h>

namespace srtm {
// constants

short get_elevation(double lat, double lon, const std::string &prefix) {
  std::string tilename = hgt::tilename(lat, lon);
  std::string filename = prefix + tilename;
  hgt heightmap(filename);
  return heightmap.at(lat, lon);
}

void tile2stl(const std::string &tile, const std::string &stl) {}

} // namespace srtm

int main(int argc, const char *argv[]) {

  std::cout << srtm::get_elevation(44.684413, 6.614422, "test/") << std::endl;

  hgt tile("test/N44E006.hgt");

  tile.toASCII("test/tile.stl");
  return EXIT_SUCCESS;
}