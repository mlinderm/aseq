//
// Created by Michael Linderman on 12/16/15.
//

#pragma once

#include <boost/range/iterator_range.hpp>

#include "aseq/util/flyweight.hpp"

namespace aseq {
namespace model {

struct allele_tag {};

class Allele : public util::flyweight_string_no_track<allele_tag> {
  typedef util::flyweight_string_no_track<allele_tag> base_type;

 public:
  typedef typename base_type::initializer initializer;

  static const Allele A, G, C, T, N;
  static const Allele MISSING;
  static const Allele NON_REF;

  Allele() : base_type(MISSING) {}
  Allele(const char* bases) : base_type(bases) {}
  Allele(const std::string& bases) : base_type(bases) {}
  template <class I>
  Allele(const boost::iterator_range<I>& range)
      : base_type(range.begin(), range.end()) {}
  template <typename I>
  Allele(I f, I l)
      : base_type(f, l) {}
  Allele(const Allele& allele) : base_type(allele) {}

  bool IsSymbolic() const;
};

}  // namespace model
}  // namespace aseq