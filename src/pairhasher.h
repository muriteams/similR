#include <Rcpp.h>
using namespace Rcpp;

#ifndef PAIRHASHER_H
#define PAIRHASHER_H 1

// A pair of (signed) integers
struct Coords {
  int i;
  int j;
};

bool operator==(const Coords& p0, const Coords& p1) {
  return p0.i == p1.i && p0.j == p1.j;
}

struct PairHasher {
  std::size_t operator()(Coords const& p) const noexcept {
    
    std::hash<int> hasher;
    std::size_t hash = hasher(p.i);
    
    // ^ makes bitwise XOR
    // 0x9e3779b9 is a 32 bit constant (comes from the golden ratio)
    // << is a shift operator, something like lhs * 2^(rhs)
    return hash ^ (hasher(p.j) + 0x9e3779b9 + (hash << 6) + (hash >> 2));
    
  }
};

#endif