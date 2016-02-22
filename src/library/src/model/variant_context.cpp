//
// Created by Michael Linderman on 12/17/15.
//

#include "aseq/model/variant_context.hpp"

namespace aseq {
namespace model {

VariantContext::VariantContext(impl::VariantContextData &&data)
    : util::HasAttributes(std::move(data.attrs_)),
      contig_(std::move(data.contig_)),
      pos_(std::move(data.pos_)),
      end_(std::move(data.end_)),
      ref_(std::move(data.ref_)),
      alts_(std::move(data.alts_)),
      ids_(std::move(data.ids_)),
      qual_(std::move(data.qual_)),
      filters_(std::move(data.filters_)) {
  // TODO: Validate data
}

VariantContext::VariantContext(const Contig &contig, int64_t pos, const Allele &ref)
    : VariantContext(contig, pos, ref, {}) {}

VariantContext::VariantContext(const Contig &contig, int64_t pos, const Allele &ref,
                               const Allele &alt)
    : VariantContext(contig, pos, ref, {alt}) {}

VariantContext::VariantContext(const Contig &contig, Pos pos, const Allele &ref,
                               std::initializer_list<Allele> alts)
    : contig_(contig), pos_(pos), ref_(ref), alts_(alts) {
  end_ = pos_ + ref_.size() - 1;
  // TODO: Validate data
}

std::ostream &operator<<(std::ostream &os, const VariantContext &vc) {
  os << vc.contig_ << ":" << vc.pos_ << vc.ref_;
  return os;
}

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
    VariantContext::Alleles left_alts(left.alts_), right_alts(right.alts_);
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