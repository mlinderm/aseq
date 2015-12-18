//
// Created by Michael Linderman on 12/18/15.
//

#include <gtest/gtest.h>

#include "aseq/model/variant_context.hpp"

using namespace aseq::model;

TEST(VariantCompareTest, OrdersPreciseRefOnlyVariants) {
  CompareVariants cmp;
  VariantContext a("1", 100, Allele::A);
  {
    VariantContext b("1", 100, Allele::A);
    EXPECT_EQ(CompareVariants::result_type::EQUAL, cmp(a, b));
  }
  {
    VariantContext b("1", 101, Allele::A);
    EXPECT_EQ(CompareVariants::result_type::BEFORE, cmp(a, b));
  }
  {
    VariantContext b("1", 99, Allele::A);
    EXPECT_EQ(CompareVariants::result_type::AFTER, cmp(a, b));
  }
}

TEST(VariantCompareTest, OrdersPreciseBiAndMultiallelicVariants) {
  CompareVariants cmp;
  VariantContext a("1", 100, Allele::A, Allele::T);
  {
    VariantContext b("1", 100, Allele::A, Allele::T);
    EXPECT_EQ(CompareVariants::result_type::EQUAL, cmp(a, b));
  }
  {
    VariantContext b("1", 100, Allele::A, {Allele::T, Allele::C}),
        c("1", 100, Allele::A, {Allele::T, Allele::C});
    EXPECT_EQ(CompareVariants::result_type::EQUAL, cmp(b, c));
  }
  {
    VariantContext b("1", 100, Allele::A, {Allele::T, Allele::C}),
        c("1", 100, Allele::A, {Allele::C, Allele::T});
    EXPECT_EQ(CompareVariants::result_type::EQUAL, cmp(b, c));
  }
  {  // 'a' is a subset of 'b'
    VariantContext b("1", 100, Allele::A, {Allele::T, Allele::C});
    EXPECT_EQ(CompareVariants::result_type::SUBSET, cmp(a, b));
  }
  {  // 'b' is a superset of 'a'
    VariantContext b("1", 100, Allele::A, {Allele::T, Allele::C});
    EXPECT_EQ(CompareVariants::result_type::SUPERSET, cmp(b, a));
  }
}
