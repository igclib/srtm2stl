#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <getopt.h>

#include <terrain.h>
#include <triangle.h>
#include <vec3.h>

const std::string PROGRAM_NAME = "srtm2stl";
const std::string VERSION = "v0.1.0";

void usage() { std::cerr << PROGRAM_NAME << " " << VERSION << std::endl; }

int main(int argc, char *argv[]) {
  switch (getopt(argc, argv, "ho:")) {
  case '?':
  case 'h':
  default:
    usage();
    return EXIT_SUCCESS;

  case -1:
    break;
  }

  double lat = 44.6;
  double lon = 6.2;
  double fast_alt = terrain::elevation(lat, lon, "test/");

  return EXIT_SUCCESS;
}