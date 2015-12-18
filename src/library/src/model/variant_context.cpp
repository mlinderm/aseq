//
// Created by Michael Linderman on 12/17/15.
//

#include "aseq/model/variant_context.hpp"

namespace aseq {
namespace model {

VariantContext::VariantContext(const Contig &contig, int64_t pos, const Allele &ref)
    : VariantContext(contig, pos, ref, static_cast<Allele *>(nullptr),
                     static_cast<Allele *>(nullptr)) {}

VariantContext::VariantContext(const Contig &contig, int64_t pos, const Allele &ref,
                               const Allele &alt)
    : VariantContext(contig, pos, ref, &alt, &alt + 1) {}

VariantContext::VariantContext(const Contig &contig, int64_t pos, const Allele &ref,
                               std::initializer_list<Allele> alts)
    : VariantContext(contig, pos, ref, alts.begin(), alts.end()) {}

CompareVariants::result_type CompareVariants::operator()(const VariantContext &left,
                                                         const VariantContext &right) const {
  if (left.contig_ != right.contig_) {
    return (left.contig_ < right.contig_) ? result_type::BEFORE : result_type::AFTER;
  } else if (left.pos_ != right.pos_) {
    return (left.pos_ < right.pos_) ? result_type::BEFORE : result_type::AFTER;
  } else if (left.end_ != right.end_) {
    return (left.end_ < right.end_) ? result_type::BEFORE : result_type::AFTER;
  } else if (left.alts_ != right.alts_) {
    // Check for different kinds of overlap (including partial overlap)
    VariantContext::Alts left_alts(left.alts_), right_alts(right.alts_);
    std::sort(left_alts.begin(), left_alts.end());
    std::sort(right_alts.begin(), right_alts.end());
    if (left_alts == right_alts) {
      // Catch if they are differently sorted
      return result_type::EQUAL;
    } else if (std::includes(right_alts.begin(), right_alts.end(), left_alts.begin(),
                             left_alts.end())) {
      return result_type::SUBSET;
    } else if (std::includes(left_alts.begin(), left_alts.end(), right_alts.begin(),
                             right_alts.end())) {
      return result_type::SUPERSET;
    } else {
      assert(false);  // TODO: handle zero to partial overlap
    }
    return result_type::SUBSET;
  } else {
    return result_type::EQUAL;
  }
}

}  // namespace model
}  // namespace aseq