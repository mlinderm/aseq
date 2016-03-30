//
// Created by Michael Linderman on 3/6/16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "aseq/io/vcf.hpp"
#include "aseq/model/variant_context.hpp"
#include "aseq/algorithm/variant.hpp"
#include "aseq/io/reference-mock.hpp"

using namespace aseq::model;
using namespace aseq::algorithm;

TEST(VariantUpdateRefAllele, UpdatesNAllele) {
  using aseq::util::Attributes;
  using aseq::io::VCFHeader;
  aseq::io::testing::MockReferenceSource ref;
  EXPECT_CALL(ref, Sequence(Contig("1"), 2, 2)).WillRepeatedly(::testing::Return("A"));

  EXPECT_NO_THROW({
    aseq::model::impl::VariantContextData data;
    data.contig_ = "1";
    data.pos_ = 2;
    data.ref_ = Allele::N;
    data.alts_ = VariantContext::Alleles({"<DEL>"});
    data.attrs_.emplace(VCFHeader::INFO::END, Attributes::Integer(100));
    VariantContext cxt(std::move(data));
    VariantContext new_cxt = UpdateREFAllele(ref, std::move(cxt));
    EXPECT_EQ(Allele::A, new_cxt.ref());
    EXPECT_EQ(Allele("<DEL>"), new_cxt.alt(0));
  });
}

TEST(LeftAlignAndTrimVariantAllelesTest, LeftAlignAndTrimsAlleles) {
  aseq::io::testing::MockReferenceSource ref;
  EXPECT_CALL(ref, Sequence(Contig("1"), 1, 5)).Times(1).WillOnce(::testing::Return("GGGCA"));

  EXPECT_NO_THROW({
    VariantContext cxt("1", 6, "CAC", Allele::C);
    VariantContext new_cxt = LeftAlignAndTrimAlleles(ref, std::move(cxt));

    EXPECT_EQ(3, new_cxt.pos());
    EXPECT_EQ(Allele("GCA"), new_cxt.ref());
    EXPECT_EQ(Allele::G, new_cxt.alt(0));
  });
}