//
// Created by Michael Linderman on 3/6/16.
//

#include "aseq/util/exception.hpp"
#include "aseq/algorithm/variant.hpp"

namespace aseq {
namespace algorithm {

using model::Allele;
using model::VariantContext;

VariantContext UpdateREFAllele(io::ReferenceSource& ref, VariantContext&& cxt) {
  if (cxt.ref() == model::Allele::N) {
    // Update reference allele
    std::string seq = ref.Sequence(cxt.contig(), cxt.pos(), cxt.end());
    return VariantContext(std::move(cxt), model::Allele(seq));
  } else
    return std::move(cxt);
}
}
}