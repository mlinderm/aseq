//
// Created by Michael Linderman on 3/6/16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "aseq/model/variant_context.hpp"
#include "aseq/algorithm/variant.hpp"
#include "aseq/io/reference-mock.hpp"

using namespace aseq::algorithm;

TEST(VariantUpdateRefAllele, UpdatesNAllele) {
  using aseq::model::Allele;
  using aseq::model::Contig;
  using aseq::model::VariantContext;

  aseq::io::testing::MockReferenceSource ref;
  EXPECT_CALL(ref, Sequence(Contig("1"), 2, 2)).WillRepeatedly(::testing::Return("A"));

  EXPECT_NO_THROW({
    VariantContext cxt("1", 2, Allele::N, "<DEL>");
    VariantContext new_cxt = UpdateREFAllele(ref, std::move(cxt));
    EXPECT_EQ(Allele::A, new_cxt.ref());
    EXPECT_EQ(Allele("<DEL>"), new_cxt.alt(0));
  });
}