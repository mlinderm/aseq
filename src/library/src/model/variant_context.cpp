//
// Created by Michael Linderman on 12/17/15.
//

#include "aseq/model/variant_context.hpp"

namespace aseq {
namespace model {

constexpr AlleleIndex VariantContext::kNoCallIdx, VariantContext::kNonRefIdx,
    VariantContext::kRefIdx, VariantContext::kFirstAltIdx;

VariantContext::VariantContext(const Contig &contig, int64_t pos, const Allele &ref)
    : VariantContext(contig, pos, ref, {}) {}

VariantContext::VariantContext(const Contig &contig, int64_t pos, const Allele &ref,
                               const Allele &alt)
    : VariantContext(contig, pos, ref, {alt}) {}

VariantContext::VariantContext(const Contig &contig, Pos pos, const Allele &ref,
                               std::initializer_list<Allele> alts)
    : VariantContext(contig, pos, pos + ref.size() - 1, ref, alts) {}

VariantContext::VariantContext(const Contig &contig, Pos pos, Pos end, const Allele &ref,
                               std::initializer_list<Allele> alts)
    : HasRegion(contig, pos, end), ref_(ref), alts_(alts) {
  // TODO: Validate data
}

VariantContext::VariantContext(VariantContext &&other)
    : util::HasAttributes(other),
      HasRegion(other),
      ids_(std::move(other.ids_)),
      qual_(std::move(other.qual_)),
      filters_(std::move(other.filters_)),
      ref_(std::move(other.ref_)),
      alts_(std::move(other.alts_))
{
  genotypes_.reserve(other.genotypes_.size());
  for (Genotype &g : other.genotypes_) {
    genotypes_.emplace_back(*this, std::move(g));
  }
}

VariantContext &VariantContext::operator=(VariantContext &&rhs) {
  this->HasRegion::operator=(std::move(rhs));
  this->util::HasAttributes::operator=(std::move(rhs));
  ids_ = std::move(rhs.ids_);
  qual_ = std::move(rhs.qual_);
  filters_ = std::move(rhs.filters_);
  ref_ = std::move(rhs.ref_);
  alts_ = std::move(rhs.alts_);
  genotypes_.clear();
  for (Genotype &g : rhs.genotypes_) {
    genotypes_.emplace_back(*this, std::move(g));
  }
  return *this;
}

VariantContext::VariantContext(VariantContext &&context, const Allele &ref)
    : VariantContext(std::move(context)) {
  ref_ = ref;
}

VariantContext &VariantContext::SetFlag(VariantContext::Flags flags) {
  flags_ |= std::underlying_type<Flags>::type(flags);
  return *this;
}

bool VariantContext::IsPASSing() const {
  static VariantContext::Filter kPASS("PASS");
  return filters_.size() == 1 && filters_.front() == kPASS;
}

const Genotype &VariantContext::GetGenotype(const Sample &sample) const {
  // TODO: Sort genotypes for quicker access
  auto r = std::find_if(genotypes_.begin(), genotypes_.end(),
                        [&](const Genotype &g) { return g.sample() == sample; });
  if (r == genotypes_.end()) throw util::no_such_sample();
  return *r;
}

Genotype VariantContext::GetGenotypeOrNoCall(const Sample &sample) const {
  // TODO: Sort genotypes for quicker access
  auto r = std::find_if(genotypes_.begin(), genotypes_.end(),
                        [&](const Genotype &g) { return g.sample() == sample; });
  return (r == genotypes_.end()) ? Genotype(*this, sample, Genotype::kNone, {}) : *r;
}

void VariantContext::MergeGenotypes(Genotypes &&genotypes) {
  // Append genotypes and remove duplicates
  for (auto &gt : genotypes) {
    genotypes_.emplace_back(*this, std::move(gt));
  }
  // stable sort preferences original genotypes from the same sample
  // TODO: preference call in new genotypes over NO_CALL in original genotypes
  std::stable_sort(genotypes_.begin(), genotypes_.end(),
                   [&](const Genotype &l, const Genotype &r) { return l.sample() < r.sample(); });
  genotypes_.erase(
      std::unique(genotypes_.begin(), genotypes_.end(),
                  [&](const Genotype &l, const Genotype &r) { return l.sample() == r.sample(); }),
      genotypes_.end());
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