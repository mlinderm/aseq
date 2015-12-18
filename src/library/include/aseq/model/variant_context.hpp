//
// Created by Michael Linderman on 12/13/15.
//

#pragma once

#include "aseq/util/exception.hpp"
#include "aseq/model/contig.hpp"
#include "aseq/model/allele.hpp"

namespace aseq {
namespace model {

class CompareVariants;

class VariantContext {
 public:
  VariantContext(const Contig &contig, int64_t pos, const Allele &ref);
  VariantContext(const Contig &contig, int64_t pos, const Allele &ref, const Allele &alt);
  VariantContext(const Contig &contig, int64_t pos, const Allele &ref,
                 std::initializer_list<Allele> alts);

  template <typename Iterator>
  VariantContext(const Contig &contig, int64_t pos, const Allele &ref, Iterator alt_first,
                 Iterator alt_last);

  const Contig &contig() const { return contig_; }

  int64_t pos() const { return pos_; }
  int64_t end() const { return end_; }

 private:
  typedef std::vector<Allele> Alts;

  Contig contig_;
  int64_t pos_, end_;
  Allele ref_;
  Alts alts_;

  friend class CompareVariants;
};

template <typename Iterator>
VariantContext::VariantContext(const Contig &contig, int64_t pos, const Allele &ref,
                               Iterator alt_first, Iterator alt_last)
    : contig_(contig), pos_(pos), ref_(ref), alts_(alt_first, alt_last) {
  if (std::any_of(alts_.begin(), alts_.end(), [](const Allele &a) { return a.IsSymbolic(); })) {
    throw util::invalid_argument()
        << util::error_message("VariantContext with symbolic allele requires explicit END");
  }
  end_ = pos_ + ref_.size() - 1;  // Recall coordinates are 1-indexed
}

class CompareVariants {
 public:
  typedef VariantContext first_argument_type;
  typedef VariantContext second_argument_type;
  enum class result_type : int { BEFORE = -2, SUBSET = -1, EQUAL = 0, SUPERSET = 1, AFTER = 2 };

  result_type operator()(const VariantContext &left, const VariantContext &right) const;
};
}
}
