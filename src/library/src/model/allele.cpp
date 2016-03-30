//
// Created by Michael Linderman on 12/16/15.
//

#include "aseq/model/allele.hpp"

namespace aseq {
namespace model {

static Allele::initializer allele_init;

const Allele Allele::A("A"), Allele::G("G"), Allele::C("C"), Allele::T("T"), Allele::N("N");
const Allele Allele::MISSING(".");
const Allele Allele::NON_REF("<NON_REF>");

bool Allele::IsSymbolic() const { return (front() == '<') && (back() == '>'); }

Allele Allele::SubAllele(size_t pos, size_t len) const {
  return Allele(substr(pos, len));
}

}  // namespace model
}  // namespace aseq
