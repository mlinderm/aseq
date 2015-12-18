//
// Created by Michael Linderman on 12/16/15.
//

#include <gtest/gtest.h>

#include "aseq/model/allele.hpp"

using namespace aseq::model;

TEST(AllelePropertiesTest, IdentifiesSymbolicAllele) {
  EXPECT_TRUE(Allele("<DEL>").IsSymbolic());
  EXPECT_FALSE(Allele("A").IsSymbolic());
}
