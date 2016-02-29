//
// Created by Michael Linderman on 2/28/16.
//

#pragma once

#include "aseq/util/flyweight.hpp"

namespace aseq {
namespace model {

struct contig_tag {};
typedef util::flyweight_string_no_track<contig_tag> Contig;

typedef int64_t Pos;

class HasRegion {
 public:
  HasRegion(const Contig &contig, Pos pos, Pos end) : contig_(contig), pos_(pos), end_(end) {}

  const Contig &contig() const { return contig_; }
  Pos pos() const { return pos_; }
  Pos end() const { return end_; }
  size_t size() const { return end_ - pos_ + 1; }

 protected:
  Contig contig_;
  Pos pos_, end_;
};
}
}