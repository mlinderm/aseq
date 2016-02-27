//
// Created by Michael Linderman on 12/17/15.
//

#pragma once

#include <string>

#include "aseq/io/reference.hpp"
#include "aseq/model/variant_context.hpp"

namespace aseq {
namespace algorithm {

std::string Consensus(io::ReferenceSource& ref, const model::VariantContext& cxt, uint64_t flank);
}
}