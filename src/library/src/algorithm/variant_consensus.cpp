//
// Created by Michael Linderman on 2/24/16.
//

#include <aseq/io/vcf.hpp>
#include "aseq/util/exception.hpp"
#include "aseq/algorithm/variant.hpp"

namespace aseq {
namespace algorithm {
std::string Consensus(io::ReferenceSource& ref, const model::VariantContext& cxt, uint64_t flank) {
  if (!cxt.IsBiallelic()) {
    throw util::invalid_argument()
        << util::error_message("Consensus is only valid for biallelic variants");
  }

  std::string result;
  result.reserve(2 * flank + cxt.size());

  // 1. Left-flank
  if (flank > 0) {
// TODO: Adjust flank when there is a padding base
      result = ref.Sequence(cxt.contig(), cxt.pos() - flank, cxt.pos() - 1);
  }
  // 2. Replacement
  const auto& alt = cxt.alt(0);
  if (!alt.IsSymbolic()) {
    result.append(alt);
  } else if (cxt.HasAttribute(io::VCFHeader::INFO::SVTYPE)) {
    const auto& type = cxt.GetAttribute<util::Attributes::String>(io::VCFHeader::INFO::SVTYPE);
    if (type == "DEL") {
      // For symbolic deletions we don't write out any alternate bases, just the padding reference
      // base
      const auto& ref = cxt.ref();
      if (ref.size() != 1) {
        throw util::invalid_argument() << util::error_message("Consensus assumes single REF base");
      }
      result.append(ref);
    } else {
      throw util::invalid_argument()
          << util::error_message("Unsupported symbolic variant for Consensus");
    }
  } else {
    throw util::invalid_argument()
        << util::error_message("Unsupported symbolic variant for Consensus");
  }
  // 3. Right-flank
  if (flank > 0) {
    result.append(ref.Sequence(cxt.contig(), cxt.end() + 1, cxt.end() + flank));
  }

  return result;
}
}
}
