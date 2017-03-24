//
// Created by Michael Linderman on 12/13/15.
//

#pragma once

#include <iosfwd>
#include <vector>

#include <boost/optional.hpp>

#include "aseq/model/allele.hpp"
#include "aseq/model/genotype.hpp"
#include "aseq/model/region.hpp"
#include "aseq/util/attributes.hpp"
#include "aseq/util/exception.hpp"

namespace aseq {
namespace model {

class CompareVariants;

class VariantContext : public util::HasAttributes, public HasRegion {
 public:
  static_assert(std::is_signed<AlleleIndex>::value, "AlleleIndex must be signed type");
  static constexpr AlleleIndex kNoCallIdx = -1, kNonRefIdx = -2, kRefIdx = 0, kFirstAltIdx = 1;
  static constexpr size_t kMaxAlleles = std::numeric_limits<AlleleIndex>::max();

  typedef std::vector<Allele> Alleles;
  typedef std::vector<std::string> IDs;
  typedef boost::optional<float> Qual;
  typedef util::Attributes::key_type Filter;
  typedef std::vector<Filter> Filters;
  typedef std::vector<Genotype> Genotypes;

  enum class Flags : unsigned int { kNone = 0, kSymbolic = 0x1 };

 public:
  VariantContext(const Contig &contig, Pos pos, const Allele &ref);
  VariantContext(const Contig &contig, Pos pos, const Allele &ref, const Allele &alt);
  VariantContext(const Contig &contig, Pos pos, const Allele &ref,
                 std::initializer_list<Allele> alts);
  VariantContext(const Contig &contig, Pos pos, Pos end, const Allele &ref,
                 std::initializer_list<Allele> alts);
  template <typename Iterator>
  VariantContext(const Contig &contig, Pos pos, Pos end, const Allele &ref, Iterator alt_begin,
                 Iterator alt_end);

  VariantContext(VariantContext &&);
  VariantContext &operator=(VariantContext &&);

  VariantContext(VariantContext &&, const Allele &ref);
  template <typename Iterator>
  VariantContext(VariantContext &&, const Contig &contig, Pos pos, const Allele &ref,
                 Iterator alt_begin, Iterator alt_end);

  // Flags
  bool flags(Flags flags) const { return flags_ & std::underlying_type<Flags>::type(flags); }
  VariantContext &SetFlag(Flags flags);

  const Allele &ref() const { return ref_; }
  const Alleles &alts() const { return alts_; }
  const Allele &alt(Alleles::size_type idx) const { return alts_.at(idx); }

  size_t NumAltAlleles() const { return alts_.size(); }

  bool IsMonoallelic() const { return alts_.empty(); }
  bool IsBiallelic() const { return alts_.size() == 1; }
  bool IsMultiAllelic() const { return alts_.size() > 1; }

  const IDs &ids() const { return ids_; }

  const Qual &qual() const { return qual_; }
  bool HasQual() const { return static_cast<bool>(qual_); }
  void SetQual(Qual::argument_type qual) { qual_ = qual; }
  void SetQual(const Qual &qual) { qual_ = qual; }

  const Filters &filters() const { return filters_; }
  bool HasFilter() const { return !filters_.empty(); }
  bool IsPASSing() const;

  const Genotypes &genotypes() const { return genotypes_; }
  Genotypes &&genotypes() { return std::move(genotypes_); }

  size_t NumGenotypes() const { return genotypes_.size(); }
  const Genotype &GetGenotype(const Sample &sample) const;
  Genotype GetGenotypeOrNoCall(const Sample &sample) const;

  template <typename... Args>
  Genotype &AddGenotype(Args &&... args) {
    genotypes_.emplace_back(*this, std::forward<Args>(args)...);
    return genotypes_.back();
  }
  void MergeGenotypes(Genotypes &&);

  friend std::ostream &operator<<(std::ostream &, const VariantContext &);

 public:
  // VCF-style Context Fields, including INFO attributes
  IDs ids_;
  Qual qual_;
  Filters filters_;

 private:
  std::underlying_type<Flags>::type flags_;

  Allele ref_;
  Alleles alts_;
  Genotypes genotypes_;

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

template <typename Iterator>
inline VariantContext::VariantContext(VariantContext &&context, const Contig &contig, Pos pos,
                                      const Allele &ref, Iterator alt_begin, Iterator alt_end)
    : VariantContext(std::move(context)) {
  this->HasRegion::operator=(HasRegion(contig, pos, pos + ref.size() - 1));
  ref_ = ref;
  alts_.assign(alt_begin, alt_end);
}

template <typename Iterator>
inline VariantContext::VariantContext(const Contig &contig, Pos pos, Pos end, const Allele &ref,
                                      Iterator alt_begin, Iterator alt_end)
    : HasRegion(contig, pos, end) {
  ref_ = ref;
  alts_.assign(alt_begin, alt_end);
}

inline VariantContext::Flags operator|(VariantContext::Flags lhs, VariantContext::Flags rhs) {
  return VariantContext::Flags(std::underlying_type<VariantContext::Flags>::type(lhs) |
                               std::underlying_type<VariantContext::Flags>::type(rhs));
}
}  // namespace model
}  // namespace aseq
