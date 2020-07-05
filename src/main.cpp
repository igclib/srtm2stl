#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <getopt.h>

#include <hgt.h>
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

  std::cout << hgt::elevation(44.684413, 6.614422, "test/") << std::endl;

  // hgt tile("test/N44E006.hgt");
  // tile.toASCII("test/tilenosyncv2.stl");

  return EXIT_SUCCESS;
}