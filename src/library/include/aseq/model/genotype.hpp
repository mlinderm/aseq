//
// Created by Michael Linderman on 3/1/16.
//

#pragma once

#include <aseq/util/exception.hpp>
#include "aseq/util/flyweight.hpp"
#include "aseq/util/attributes.hpp"

namespace aseq {
namespace model {

class VariantContext;

typedef int AlleleIndex;

struct sample_tag {};
typedef util::flyweight_string_no_track<sample_tag> Sample;

namespace impl {

struct PhasedIndices {
  typedef std::vector<AlleleIndex> Indices;

  PhasedIndices() : phased_(false) {}
  explicit PhasedIndices(bool phased) : phased_(phased) {}
  explicit PhasedIndices(AlleleIndex a) : PhasedIndices(false, {a}) {}
  PhasedIndices(bool phased, AlleleIndex a1, AlleleIndex a2) : PhasedIndices(phased, {a1, a2}) {}
  PhasedIndices(bool phased, std::initializer_list<AlleleIndex> indices);

  size_t Ploidy() const { return indices_.size(); }
  bool operator==(const PhasedIndices& rhs) const;

  bool phased_;
  Indices indices_;
};

std::size_t hash_value(const PhasedIndices&);
std::ostream& operator<<(std::ostream& os, const PhasedIndices& v);

}  // namespace impl

class Genotype : public util::HasAttributes {
 public:
  typedef typename boost::flyweight<impl::PhasedIndices, boost::flyweights::no_tracking> Alleles;

  // Commonly observed genotypes (*P are phased)
  static const Alleles kNone, kNoCallNoCall, kRefRef, kRefAlt, kAltAlt, kRefAltP, kAltRefP;

  Genotype() = delete;
  Genotype(const VariantContext& variant, const Sample& sample)
      : Genotype(variant, sample, kNone, {}) {}
  Genotype(const VariantContext& variant, const Sample& sample, util::Attributes&& attr)
      : Genotype(variant, sample, kNone, std::move(attr)) {}
  Genotype(const VariantContext& variant, const Sample& sample, const Alleles& alleles)
      : Genotype(variant, sample, kNone, {}) {}
  Genotype(const VariantContext& variant, const Sample& sample, const Alleles& alleles,
           util::Attributes&& attr);
  Genotype(const VariantContext&, Genotype&&);

  const VariantContext& variant_context() const { return *variant_; }
  const Sample& sample() const { return sample_; }
  const Alleles& alleles() const { return alleles_; }

  bool Phased() const { return indices().phased_; }

  size_t Ploidy() const { return indices().Ploidy(); }
  bool Haploid() const { return Ploidy() == 1; }
  bool Diploid() const { return Ploidy() == 2; }

  friend void swap(Genotype& a, Genotype& b) {
    if (a.variant_ != b.variant_) {
      throw util::invalid_argument();
    }

    using std::swap;
    swap(static_cast<util::HasAttributes&>(a), static_cast<util::HasAttributes&>(b));
    swap(a.sample_, b.sample_);
    swap(a.alleles_, b.alleles_);
  }

 private:
  const impl::PhasedIndices& indices() const { return alleles_.get(); }

  const VariantContext* variant_;
  Sample sample_;
  Alleles alleles_;
};
}
}