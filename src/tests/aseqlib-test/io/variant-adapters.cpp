//
// Created by Michael Linderman on 4/5/16.
//

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "aseq/io/variant-adapters.hpp"

using namespace aseq::io;
using namespace aseq::model;
using Attributes = aseq::util::Attributes;

TEST(VariantMergingSourceTest, MergesPreciseSitesOnlyVariants) {
  auto vcf =
      "##fileformat=VCFv4.2\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
      "1\t1\t.\tA\tT\t.\t.\t.";
  std::stringstream variants1(vcf), variants2(vcf);
  auto source = VariantMergeSourceInterface::MakeMergeVariantSource(
      VariantSourceInterface::MakeVariantSource(variants1),
      VariantSourceInterface::MakeVariantSource(variants2));
  ASSERT_TRUE(source);
  EXPECT_FALSE(source->IsIndexed());

  auto& vcf_header = dynamic_cast<const VCFHeader&>(source->header());
  EXPECT_TRUE(vcf_header.HasINFOField(source->source_key()));

  auto v = source->NextVariant();
  ASSERT_TRUE(v);

  EXPECT_NO_THROW({
    EXPECT_EQ(source->merged_label(),
              v->GetAttribute<VariantMergeSourceInterface::SourceLabel>(source->source_key()));
  });

  v = source->NextVariant();
  EXPECT_FALSE(v);
}

TEST(VariantMergingSourceTest, MergesSitesOnlyVariantsInOrder) {
  using SourceLabel = VariantMergeSourceInterface::SourceLabel;

  auto vcf1 =
      "##fileformat=VCFv4.2\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
      "1\t1\t.\tA\tT\t.\t.\t.";
  auto vcf2 =
      "##fileformat=VCFv4.2\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
      "1\t2\t.\tC\tG\t.\t.\t.";
  std::stringstream variants1(vcf1), variants2(vcf2);
  auto source = VariantMergeSourceInterface::MakeMergeVariantSource(
      VariantSourceInterface::MakeVariantSource(variants1),
      VariantSourceInterface::MakeVariantSource(variants2));
  ASSERT_TRUE(source);

  auto v = source->NextVariant();
  ASSERT_TRUE(v);

  EXPECT_EQ(1, v->pos());
  EXPECT_EQ(source->source1_label(), v->GetAttribute<SourceLabel>(source->source_key()));

  v = source->NextVariant();
  ASSERT_TRUE(v);

  EXPECT_EQ(2, v->pos());
  EXPECT_EQ(source->source2_label(), v->GetAttribute<SourceLabel>(source->source_key()));

  v = source->NextVariant();
  EXPECT_FALSE(v);
}

TEST(VariantMergingSourceTest, CombinesDifferentComponentsOfSitesOnlyVariants) {
  using SourceLabel = VariantMergeSourceInterface::SourceLabel;

  auto vcf1 =
      "##fileformat=VCFv4.2\n"
      "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
      "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
      "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
      "1\t1\t.\tA\tT\t.\t.\tAC=1;AN=2;DP=10";
  auto vcf2 =
      "##fileformat=VCFv4.2\n"
      "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
      "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
      "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
      "1\t1\t.\tA\tT\t100\t.\tAC=2;AN=4;DP=20";
  std::stringstream variants1(vcf1), variants2(vcf2);
  auto source = VariantMergeSourceInterface::MakeMergeVariantSource(
      VariantSourceInterface::MakeVariantSource(variants1),
      VariantSourceInterface::MakeVariantSource(variants2));
  ASSERT_TRUE(source);

  auto v = source->NextVariant();
  ASSERT_TRUE(v);
  EXPECT_NO_THROW({
    EXPECT_EQ(1, v->pos());
    EXPECT_EQ(VariantContext::Qual(100.0), v->qual());
    EXPECT_EQ(source->merged_label(), v->GetAttribute<SourceLabel>(source->source_key()));
    EXPECT_FALSE(v->HasAttribute(VCFHeader::INFO::DP)) << "Didn't drop mismatching INFO fields";
    // EXPECT_EQ(Attributes::Integers{3},
    // v->GetAttribute<Attributes::Integers>(VCFHeader::INFO::AC));
    // EXPECT_EQ(6, v->GetAttribute<Attributes::Integer>(VCFHeader::INFO::AN));
  });

  v = source->NextVariant();
  EXPECT_FALSE(v);
}

TEST(VariantMergingSourceTest, CombinesVariantsWithIdenticalSamples) {
  using SourceLabel = VariantMergeSourceInterface::SourceLabel;

  auto vcf1 =
      "##fileformat=VCFv4.2\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
      "1\t1\t.\tA\tT\t200\t.\t.\tGT\t0/0";
  auto vcf2 =
      "##fileformat=VCFv4.2\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
      "1\t1\t.\tA\tT\t100\t.\t.\tGT:GQ\t0/1:20";
  std::stringstream variants1(vcf1), variants2(vcf2);
  auto source = VariantMergeSourceInterface::MakeMergeVariantSource(
      VariantSourceInterface::MakeVariantSource(variants1),
      VariantSourceInterface::MakeVariantSource(variants2));
  ASSERT_TRUE(source);

  auto& header = dynamic_cast<const VCFHeader&>(source->header());
  ASSERT_EQ(1, header.NumSamples());

  auto v = source->NextVariant();
  ASSERT_TRUE(v);
  EXPECT_NO_THROW({
    EXPECT_EQ(VariantContext::Qual(100.0), v->qual());
    EXPECT_EQ(source->merged_label(), v->GetAttribute<SourceLabel>(source->source_key()));
    ASSERT_EQ(1, v->NumGenotypes());
    auto& gt = v->GetGenotype("Sample1");
    EXPECT_EQ(Genotype::kRefRef, gt.alleles());
    EXPECT_FALSE(gt.HasAttribute(VCFHeader::FORMAT::GQ))
        << "Shouldn't bring over attributes from unselected genotypes";
  });

  v = source->NextVariant();
  EXPECT_FALSE(v);
}

TEST(VariantMergingSourceTest, CombinesVariantsWithDifferentSamples) {
  using SourceLabel = VariantMergeSourceInterface::SourceLabel;

  auto vcf1 =
      "##fileformat=VCFv4.2\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
      "1\t1\t.\tA\tT\t100\t.\t.\tGT\t0/0";
  auto vcf2 =
      "##fileformat=VCFv4.2\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample2\n"
      "1\t1\t.\tA\tT\t.\t.\t.\tGT:GQ\t0/1:20";
  std::stringstream variants1(vcf1), variants2(vcf2);
  auto source = VariantMergeSourceInterface::MakeMergeVariantSource(
      VariantSourceInterface::MakeVariantSource(variants1),
      VariantSourceInterface::MakeVariantSource(variants2));
  ASSERT_TRUE(source);

  auto& header = dynamic_cast<const VCFHeader&>(source->header());
  ASSERT_EQ(2, header.NumSamples());

  auto v = source->NextVariant();
  ASSERT_TRUE(v);
  EXPECT_NO_THROW({
    EXPECT_EQ(VariantContext::Qual(100.0), v->qual());
    EXPECT_EQ(source->merged_label(), v->GetAttribute<SourceLabel>(source->source_key()));
    ASSERT_EQ(2, v->NumGenotypes());

    auto& gt1 = v->GetGenotype("Sample1");
    EXPECT_EQ(Genotype::kRefRef, gt1.alleles());
    EXPECT_FALSE(gt1.HasAttribute(VCFHeader::FORMAT::GQ));

    auto& gt2 = v->GetGenotype("Sample2");
    EXPECT_EQ(Genotype::kRefAlt, gt2.alleles());
    EXPECT_TRUE(gt2.HasAttribute(VCFHeader::FORMAT::GQ));
  });

  v = source->NextVariant();
  EXPECT_FALSE(v);
}

TEST(VariantMergingSourceTest, CombinesDifferentVariantsWithDifferentSamples) {
  using SourceLabel = VariantMergeSourceInterface::SourceLabel;

  auto vcf1 =
      "##fileformat=VCFv4.2\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
      "1\t1\t.\tA\tT\t.\t.\t.\tGT\t0/0";
  auto vcf2 =
      "##fileformat=VCFv4.2\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample2\n"
      "1\t2\t.\tA\tT\t.\t.\t.\tGT:GQ\t0/1:20";
  std::stringstream variants1(vcf1), variants2(vcf2);
  auto source = VariantMergeSourceInterface::MakeMergeVariantSource(
      VariantSourceInterface::MakeVariantSource(variants1),
      VariantSourceInterface::MakeVariantSource(variants2));
  ASSERT_TRUE(source);

  auto& header = dynamic_cast<const VCFHeader&>(source->header());
  ASSERT_EQ(2, header.NumSamples());

  auto v = source->NextVariant();
  ASSERT_TRUE(v);

  EXPECT_NO_THROW({
    EXPECT_EQ(source->source1_label(), v->GetAttribute<SourceLabel>(source->source_key()));
    ASSERT_EQ(2, v->NumGenotypes());

    auto& gt1 = v->GetGenotype("Sample1");
    EXPECT_EQ(Genotype::kRefRef, gt1.alleles());
    EXPECT_FALSE(gt1.HasAttribute(VCFHeader::FORMAT::GQ));

    auto& gt2 = v->GetGenotype("Sample2");
    EXPECT_EQ(Genotype::kNone, gt2.alleles());
    EXPECT_FALSE(gt2.HasAttribute(VCFHeader::FORMAT::GQ));
  });

  v = source->NextVariant();
  ASSERT_TRUE(v);

  EXPECT_NO_THROW({
    EXPECT_EQ(source->source2_label(), v->GetAttribute<SourceLabel>(source->source_key()));
    ASSERT_EQ(2, v->NumGenotypes());

    auto& gt1 = v->GetGenotype("Sample1");
    EXPECT_EQ(Genotype::kNone, gt1.alleles());
    EXPECT_FALSE(gt1.HasAttribute(VCFHeader::FORMAT::GQ));

    auto& gt2 = v->GetGenotype("Sample2");
    EXPECT_EQ(Genotype::kRefAlt, gt2.alleles());
    EXPECT_TRUE(gt2.HasAttribute(VCFHeader::FORMAT::GQ));
  });

  v = source->NextVariant();
  EXPECT_FALSE(v);
}

TEST(VariantMergingSourceTest, CombinesVariantsWithPartiallyOverlappingSamples) {
  using SourceLabel = VariantMergeSourceInterface::SourceLabel;

  auto vcf1 =
      "##fileformat=VCFv4.2\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n"
      "1\t1\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0";
  auto vcf2 =
      "##fileformat=VCFv4.2\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample2\tSample3\n"
      "1\t1\t.\tA\tT\t.\t.\t.\tGT:GQ\t0/1:20\t0/1";
  std::stringstream variants1(vcf1), variants2(vcf2);
  auto source = VariantMergeSourceInterface::MakeMergeVariantSource(
      VariantSourceInterface::MakeVariantSource(variants1),
      VariantSourceInterface::MakeVariantSource(variants2));
  ASSERT_TRUE(source);

  auto& header = dynamic_cast<const VCFHeader&>(source->header());
  ASSERT_EQ(3, header.NumSamples());

  auto v = source->NextVariant();
  ASSERT_TRUE(v);
  EXPECT_NO_THROW({
    EXPECT_EQ(source->merged_label(), v->GetAttribute<SourceLabel>(source->source_key()));
    ASSERT_EQ(3, v->NumGenotypes());

    auto& gt1 = v->GetGenotype("Sample1");
    EXPECT_EQ(Genotype::kRefRef, gt1.alleles());
    EXPECT_FALSE(gt1.HasAttribute(VCFHeader::FORMAT::GQ));

    auto& gt2 = v->GetGenotype("Sample2");
    EXPECT_EQ(Genotype::kRefRef, gt2.alleles());
    EXPECT_FALSE(gt2.HasAttribute(VCFHeader::FORMAT::GQ));

    auto& gt3 = v->GetGenotype("Sample3");
    EXPECT_EQ(Genotype::kRefAlt, gt3.alleles());
    EXPECT_FALSE(gt3.HasAttribute(VCFHeader::FORMAT::GQ));
  });

  v = source->NextVariant();
  EXPECT_FALSE(v);
}