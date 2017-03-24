//
// Created by Michael Linderman on 2/24/16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "aseq/io/vcf.hpp"
#include "aseq/model/variant_context.hpp"
#include "aseq/algorithm/variant.hpp"
#include "aseq/io/reference-mock.hpp"

using namespace aseq::util;
using namespace aseq::model;
using namespace aseq::io;
using namespace aseq::io::testing;
using namespace aseq::algorithm;

TEST(VariantConsensus, RejectsNotBiAllelic) {
  MockReferenceSource ref;
  {
    VariantContext cxt(VariantContext("1", 2, Allele::A, {Allele::T, Allele::C}));
    EXPECT_THROW(Consensus(ref, cxt, 0), invalid_argument);
  }

  {
    VariantContext cxt(VariantContext("1", 2, Allele::A));
    EXPECT_THROW(Consensus(ref, cxt, 0), invalid_argument);
  }
}

TEST(VariantConsensus, AppliesSNV) {
  MockReferenceSource ref;
  EXPECT_CALL(ref, Sequence(Contig("1"), 2, 2)).WillRepeatedly(::testing::Return("A"));
  EXPECT_CALL(ref, Sequence(Contig("1"), 1, 1)).WillRepeatedly(::testing::Return("C"));
  EXPECT_CALL(ref, Sequence(Contig("1"), 3, 3)).WillRepeatedly(::testing::Return("T"));

  {
    VariantContext cxt(VariantContext("1", 2, Allele::A, Allele::T));
    EXPECT_NO_THROW({
      std::string seq = Consensus(ref, cxt, 0);
      EXPECT_EQ(std::string("T"), seq);
    });
    EXPECT_NO_THROW({
      std::string seq = Consensus(ref, cxt, 1);
      EXPECT_EQ(std::string("CTT"), seq);
    });
  }
}

TEST(VariantConsensus, AppliesPreciseDeletion) {
  MockReferenceSource ref;
  EXPECT_CALL(ref, Sequence(Contig("1"), 990, 999)).WillRepeatedly(::testing::Return("AGCTAGCTAG"));
  EXPECT_CALL(ref, Sequence(Contig("1"), 1002, 1011))
      .WillRepeatedly(::testing::Return("AGCTAGCTAG"));

  {
    VariantContext cxt(VariantContext("1", 1000, Allele("AT"), Allele::A));
    EXPECT_NO_THROW({
      std::string seq = Consensus(ref, cxt, 10);
      EXPECT_EQ(std::string("AGCTAGCTAGAAGCTAGCTAG"), seq);
    });
  }
}

TEST(VariantConsensus, AppliesPreciseInsertion) {
  MockReferenceSource ref;
  EXPECT_CALL(ref, Sequence(Contig("1"), 990, 999)).WillRepeatedly(::testing::Return("AGCTAGCTAG"));
  EXPECT_CALL(ref, Sequence(Contig("1"), 1001, 1010))
      .WillRepeatedly(::testing::Return("AGCTAGCTAG"));

  {
    VariantContext cxt(VariantContext("1", 1000, Allele::A, Allele("AT")));
    EXPECT_NO_THROW({
      std::string seq = Consensus(ref, cxt, 10);
      EXPECT_EQ(std::string("AGCTAGCTAGATAGCTAGCTAG"), seq);
    });
  }
}

TEST(VariantConsensus, AppliesSymbolicDeletion) {
  MockReferenceSource ref;
  EXPECT_CALL(ref, Sequence(Contig("1"), 990, 999)).WillRepeatedly(::testing::Return("AGCTAGCTAG"));
  EXPECT_CALL(ref, Sequence(Contig("1"), 1011, 1020))
      .WillRepeatedly(::testing::Return("AGCTAGCTAG"));

  {
    VariantContext cxt("1", 1000, 1010, Allele::C, { "<DEL>" });
    cxt.SetAttribute(VCFHeader::INFO::SVTYPE, Attributes::String("DEL"));
    EXPECT_NO_THROW({
      std::string seq = Consensus(ref, cxt, 10);
      EXPECT_EQ(std::string("AGCTAGCTAGCAGCTAGCTAG"), seq);
    });
  }
}
