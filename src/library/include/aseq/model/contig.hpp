//
// Created by Michael Linderman on 12/17/15.
//

#pragma once

#include "aseq/util/flyweight.hpp"

namespace aseq {
namespace model {

struct contig_tag {};

typedef util::flyweight_string_no_track<contig_tag> Contig;

}  // namespace model
}  // namespace aseq