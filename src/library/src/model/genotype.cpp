//
// Created by Michael Linderman on 3/2/16.
//

#include <algorithm>

#include <boost/functional/hash.hpp>

#include "aseq/model/genotype.hpp"
#include "aseq/model/variant_context.hpp"

namespace aseq {
namespace model {

static Genotype::Alleles::initializer genotype_alleles_init;

namespace impl {

PhasedIndices::PhasedIndices(bool phased, std::initializer_list<AlleleIndex> indices)
    : phased_(phased), indices_(indices) {
  if (!phased_) {
    // Sort alleles to permit linear equality comparison, set as phased if all alleles
    // are real (index >= 0) and the same
    std::sort(indices_.begin(), indices_.end());
    phased_ = std::adjacent_find(indices_.begin(), indices_.end(), [&](const AlleleIndex& a1,
                                                                       const AlleleIndex& a2) {
                return a1 < VariantContext::kRefIdx || a2 < VariantContext::kRefIdx || a1 != a2;
              }) == indices_.end();
  }
}

bool PhasedIndices::operator==(const PhasedIndices& rhs) const {
  return phased_ == rhs.phased_ && indices_ == rhs.indices_;
}

std::size_t hash_value(const PhasedIndices& indices) {
  std::size_t seed = 0;
  boost::hash_combine(seed, indices.phased_);
  boost::hash_combine(seed, indices.indices_);
  return seed;
}

std::ostream& operator<<(std::ostream& os, const PhasedIndices& v) { return os; }

}  // namespace impl

const Genotype::Alleles Genotype::kNone,
    Genotype::kNoCallNoCall(true, VariantContext::kNoCallIdx, VariantContext::kNoCallIdx),
    Genotype::kRefRef(true, VariantContext::kRefIdx, VariantContext::kRefIdx),
    Genotype::kRefAlt(false, VariantContext::kRefIdx, VariantContext::kFirstAltIdx),
    Genotype::kAltAlt(true, VariantContext::kFirstAltIdx, VariantContext::kFirstAltIdx),
    Genotype::kRefAltP(true, VariantContext::kRefIdx, VariantContext::kFirstAltIdx),
    Genotype::kAltRefP(true, VariantContext::kFirstAltIdx, VariantContext::kRefIdx);

Genotype::Genotype(const VariantContext& variant, const Sample& sample, const Alleles& alleles,
                   util::Attributes&& attr)
    : HasAttributes(std::move(attr)), variant_(variant), sample_(sample), alleles_(alleles) {}

Genotype::Genotype(const VariantContext& variant, Genotype&& other)
    : util::HasAttributes(other),
      variant_(variant),
      sample_(std::move(other.sample_)),
      alleles_(std::move(other.alleles_)) {}
}
}
