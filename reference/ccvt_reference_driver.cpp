// Runs Michael Balzer's original CCVT optimizer (GPL-3.0, github.com/michaelbalzer/ccvt)
// on an input file, using the plain Euclidean metric (not the torus from the example).
// Used only as a reference to validate the Python port.
//
// Input format (whitespace separated):
//   n_sites n_points
//   x y capacity        (n_sites lines)
//   x y                 (n_points lines)
//
// Output (stdout):
//   iteration <i> <energy>      one line per optimization pass
//   iterations <n>
//   energy <e>
//   seconds <t>
//   site <x> <y>                one line per site, in input order

#include <cstdio>
#include <fstream>
#include <time.h>
#include "ccvt_metric.h"
#include "ccvt_optimizer.h"
#include "ccvt_point.h"
#include "ccvt_site.h"

using namespace ccvt;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s input.txt\n", argv[0]);
    return 1;
  }
  std::ifstream in(argv[1]);
  int nSites, nPoints;
  in >> nSites >> nPoints;

  Site<Point2>::List sites;
  for (int i = 0; i < nSites; ++i) {
    double x, y;
    int capacity;
    in >> x >> y >> capacity;
    sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
  }
  Point2::List points;
  for (int i = 0; i < nPoints; ++i) {
    double x, y;
    in >> x >> y;
    points.push_back(Point2(x, y));
  }

  typedef Optimizer<Site<Point2>, Point2, MetricEuclidean2> Opt;
  Opt optimizer;
  MetricEuclidean2 metric;

  clock_t start = clock();
  optimizer.initialize(sites, points, metric);
  int iteration = 0;
  bool stable;
  do {
    stable = optimizer.optimize(true);
    ++iteration;
    std::printf("iteration %d %.17g\n", iteration, optimizer.energy());
  } while (!stable);
  clock_t end = clock();

  std::printf("iterations %d\n", iteration);
  std::printf("energy %.17g\n", optimizer.energy());
  std::printf("seconds %.6f\n", static_cast<double>(end - start) / CLOCKS_PER_SEC);
  const Site<Point2>::Vector& result = optimizer.sites();
  for (unsigned int i = 0; i < result.size(); ++i) {
    std::printf("site %.17g %.17g\n", result[i].location.x, result[i].location.y);
  }
  return 0;
}
