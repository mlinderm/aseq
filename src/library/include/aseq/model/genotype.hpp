//
// Created by Michael Linderman on 3/1/16.
//

#pragma once

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
  PhasedIndices(bool phased, AlleleIndex a1, AlleleIndex a2) : PhasedIndices(phased, {a1, a2}) {}
  PhasedIndices(bool phased, std::initializer_list<AlleleIndex> indices);

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
  Genotype(const VariantContext& variant, const Sample& sample, const Alleles& alleles,
           util::Attributes&& attr);
  Genotype(const VariantContext&, Genotype&&);

  const VariantContext& variant_context() const { return variant_; }
  const Sample& sample() const { return sample_; }
  const Alleles& alleles() const { return alleles_; }

  bool Phased() const { return alleles_.get().phased_; }

  size_t Ploidy() const { return indices().size(); }
  bool Haploid() const { return Ploidy() == 1; }
  bool Diploid() const { return Ploidy() == 2; }

 private:
  const impl::PhasedIndices::Indices& indices() const { return alleles_.get().indices_; }

  const VariantContext& variant_;
  Sample sample_;
  Alleles alleles_;
};
}
}