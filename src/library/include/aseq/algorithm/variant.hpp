//
// Created by Michael Linderman on 12/17/15.
//

#pragma once

#include <string>

#include "aseq/io/reference.hpp"
#include "aseq/model/variant_context.hpp"

namespace aseq {
namespace algorithm {

model::VariantContext UpdateREFAllele(io::ReferenceSource& ref, model::VariantContext&& cxt);
model::VariantContext LeftAlignAndTrimAlleles(io::ReferenceSource& ref,
                                              model::VariantContext&& cxt);
model::VariantContext Normalize(io::ReferenceSource& ref, model::VariantContext&& cxt);


std::string Consensus(io::ReferenceSource& ref, const model::VariantContext& cxt, size_t flank);
}
}