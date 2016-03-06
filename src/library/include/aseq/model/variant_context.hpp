//
// Created by Michael Linderman on 12/13/15.
//

#pragma once

#include <vector>
#include <iosfwd>

#include <boost/optional.hpp>

#include "aseq/util/exception.hpp"
#include "aseq/util/attributes.hpp"
#include "aseq/model/region.hpp"
#include "aseq/model/allele.hpp"
#include "aseq/model/genotype.hpp"

namespace aseq {
namespace model {

namespace impl {
struct VariantContextData;
}  // namespace impl

class CompareVariants;

class VariantContext : public util::HasAttributes, public HasRegion {
 public:
  static_assert(std::is_signed<AlleleIndex>::value, "AlleleIndex must be signed type");
  static constexpr AlleleIndex kNoCallIdx = -1, kNonRefIdx = -2, kRefIdx = 0, kFirstAltIdx = 1;
  static constexpr size_t kMaxAlleles = std::numeric_limits<AlleleIndex>::max();

  typedef std::vector<Allele> Alleles;
  typedef std::vector<std::string> IDs;
  typedef float Qual;
  typedef util::Attributes::key_type Filter;
  typedef std::vector<Filter> Filters;

 public:
  VariantContext(impl::VariantContextData &&data);
  VariantContext(const Contig &contig, Pos pos, const Allele &ref);
  VariantContext(const Contig &contig, Pos pos, const Allele &ref, const Allele &alt);
  VariantContext(const Contig &contig, Pos pos, const Allele &ref,
                 std::initializer_list<Allele> alts);

  VariantContext(VariantContext &&);
  VariantContext &operator=(VariantContext &&);

  const Allele &ref() const { return ref_; }
  const Allele &alt(Alleles::size_type idx) const { return alts_.at(idx); }

  bool IsMonoallelic() const { return alts_.empty(); }
  bool IsBiallelic() const { return alts_.size() == 1; }
  bool IsMultiAllelic() const { return alts_.size() > 1; }

  bool HasQual() const { return static_cast<bool>(qual_); }

  size_t NumGenotypes() const { return genotypes_.size(); }
  const Genotype &GetGenotype(const Sample &sample) const;

  template <typename... Args>
  Genotype &AddGenotype(Args &&... args) {
    genotypes_.emplace_back(*this, std::forward<Args>(args)...);
    return genotypes_.back();
  }

  friend std::ostream &operator<<(std::ostream &, const VariantContext &);

 private:
  Allele ref_;
  Alleles alts_;

  // Context Fields
  IDs ids_;
  boost::optional<Qual> qual_;
  Filters filters_;

  std::vector<Genotype> genotypes_;

  friend class CompareVariants;
};

std::ostream &operator<<(std::ostream &, const VariantContext &);

class CompareVariants {
 public:
  typedef VariantContext first_argument_type;
  typedef VariantContext second_argument_type;
  enum class result_type : int { BEFORE = -2, SUBSET = -1, EQUAL = 0, SUPERSET = 1, AFTER = 2 };

  result_type operator()(const VariantContext &left, const VariantContext &right) const;
};

namespace impl {
struct VariantContextData {
  VariantContextData() : pos_(0), end_(0), ref_(Allele::MISSING) {}

  Contig contig_;
  Pos pos_, end_;
  Allele ref_;
  VariantContext::Alleles alts_;
  VariantContext::IDs ids_;
  boost::optional<VariantContext::Qual> qual_;
  VariantContext::Filters filters_;
  util::Attributes attrs_;
};
}

}  // namespace model
}  // namespace aseq
