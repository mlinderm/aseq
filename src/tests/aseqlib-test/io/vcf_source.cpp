//
// Created by Michael Linderman on 12/18/15.
//

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "aseq/io/vcf.hpp"
#include "aseq/model/variant_context.hpp"

using namespace aseq::io;
using namespace aseq::model;

namespace aseq {
namespace io {
namespace impl {

// Expose testing helper functions
VariantContext ParseVCFVariantLine(const std::string& line, VCFHeader& header);
}  // namespace impl
}  // namespace io
}  // namespace aseq

TEST(VCFHeaderParsingTest, ParsesMinimalVCFHeader) {
  std::stringstream content("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

  VCFSource source(FileFormat::VCF4_2, ASCIILineReaderInterface::MakeLineReader(content));
  EXPECT_EQ(FileFormat::VCF4_2, source.file_format());

  auto& header = source.header();
  EXPECT_TRUE(header.HasFILTERField(VCFHeader::FILTER::PASS));
  EXPECT_TRUE(header.HasFORMATField(VCFHeader::FORMAT::GT));
  EXPECT_EQ(0, header.NumSamples());
}

TEST(VCFHeaderParsingTest, ParsesMinimalVCFHeaderWithSample) {
  std::stringstream content(
      "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1");

  VCFSource source(FileFormat::VCF4_2, ASCIILineReaderInterface::MakeLineReader(content));
  EXPECT_EQ(FileFormat::VCF4_2, source.file_format());

  auto& header = source.header();
  EXPECT_EQ(1, header.NumSamples());
}

TEST(VCFHeaderParsingTest, ParsesKeyValueLines) {
  EXPECT_NO_THROW({
    std::stringstream content(
        "##fileformat=VCFv4.2\n##reference=\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    VCFSource source(FileFormat::VCF4_2, ASCIILineReaderInterface::MakeLineReader(content));
  });
}

TEST(VCFHeaderParsingTest, ParsesFixedLengthVectorHeaderLines) {
  EXPECT_NO_THROW({
    // clang-format off
    std::stringstream content(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    );
    // clang-format on
    VCFSource source(FileFormat::VCF4_2, ASCIILineReaderInterface::MakeLineReader(content));
    auto& header = source.header();
    ASSERT_TRUE(header.HasINFOField(VCFHeader::INFO::CIEND));
    auto& field = header.INFOField(VCFHeader::INFO::CIEND);
    EXPECT_EQ(2, field.nmbr_);
  });
}

class VCFVariantParsingTest : public ::testing::Test {
 public:
  VCFVariantParsingTest() : header_(FileFormat::VCF4_2) {}

 protected:
  virtual void SetUp() {
    for (auto f : {VCFHeader::INFO::END, VCFHeader::INFO::SVTYPE, VCFHeader::INFO::SVLEN,
                   VCFHeader::INFO::CIPOS, VCFHeader::INFO::CIEND, VCFHeader::INFO::HOMLEN,
                   VCFHeader::INFO::HOMSEQ})
      header_.AddINFOField(f);

    for (auto f : {VCFHeader::FORMAT::GT, VCFHeader::FORMAT::DP, VCFHeader::FORMAT::FT,
                   VCFHeader::FORMAT::GQ, VCFHeader::FORMAT::HQ})
      header_.AddFORMATField(f);
  }

  VCFHeader header_;
};

TEST_F(VCFVariantParsingTest, ParsesMinimalSitesOnlyVariant) {
  using aseq::io::impl::ParseVCFVariantLine;

  std::string line = "1\t1\t.\tA\t.\t.\t.\t.";
  VariantContext context = ParseVCFVariantLine(line, header_);
  EXPECT_EQ(Contig("1"), context.contig());
  EXPECT_EQ(1, context.pos());
  EXPECT_EQ(1, context.end());
  EXPECT_EQ(Allele::A, context.ref());
}

TEST_F(VCFVariantParsingTest, ParsesComplexSitesOnlyVariant) {
  using aseq::io::impl::ParseVCFVariantLine;
  using aseq::util::Attributes;

  EXPECT_NO_THROW({
    std::string line("1\t1\trs100\tA\tG,<NON_REF>\t100.0\tLowQual\tEND=101");
    VariantContext context = ParseVCFVariantLine(line, header_);
    EXPECT_EQ(101, context.GetAttribute<Attributes::Integer>("END"));
  });
  EXPECT_NO_THROW({
    std::string line(
        "1\t2827694\trs2376870\tCGTGGATGCGGGGAC\tC\t.\tPASS\tSVTYPE=DEL;END=2827708;HOMLEN=1;"
        "HOMSEQ=G;SVLEN=-14");
    VariantContext context = ParseVCFVariantLine(line, header_);
  });
  EXPECT_NO_THROW({
    std::string line(
        "2\t321682\t.\tT\t<DEL>\t6\tPASS\tSVTYPE=DEL;END=321887;SVLEN=-205;CIPOS=-56,20;CIEND=-10,"
        "62");
    VariantContext context = ParseVCFVariantLine(line, header_);
    EXPECT_EQ(321887, context.end());
  });
  EXPECT_NO_THROW({
    std::string line(
        "2\t14477084\t.\tC\t<DEL:ME:ALU>\t6\tPASS\tSVTYPE=DEL;END=14477381;SVLEN=-297;CIPOS=-22,18;"
        "CIEND=-12,32");
    VariantContext context = ParseVCFVariantLine(line, header_);
  });

  EXPECT_NO_THROW({
    std::string line("2\t321681\tbnd_W\tG\tG]17:198982]\t6\tPASS\tSVTYPE=BND");
    VariantContext context = ParseVCFVariantLine(line, header_);
  });

  EXPECT_NO_THROW({
    std::string line(
        "1\t83484909\t220\tN\t<DEL>\t196.12\tPASS\tCIEND=-5,4;END=83485005;SVLEN=-96;SVTYPE=DEL");
    VariantContext context = ParseVCFVariantLine(line, header_);
    EXPECT_EQ(Attributes::Integers({-5, 4}), context.GetAttribute<Attributes::Integers>("CIEND"));
  });
}

TEST_F(VCFVariantParsingTest, RejectsIllFormedSitesOnlyVariant) {
  using aseq::io::impl::ParseVCFVariantLine;
  using aseq::util::Attributes;
  using aseq::util::file_parse_error;

  EXPECT_THROW({
    std::string line("1\t1\t.\tA\t.\t.\t.\tEND");
    VariantContext context = ParseVCFVariantLine(line, header_);
  }, file_parse_error);

  EXPECT_THROW({
    std::string line("1\t1\t.\tA\t.\t.\t.\tEND=");
    VariantContext context = ParseVCFVariantLine(line, header_);
  }, file_parse_error);

  EXPECT_THROW({
    std::string line("1\t1\t.\tA\t.\t.\t.\tEND=A");
    VariantContext context = ParseVCFVariantLine(line, header_);
  }, file_parse_error);

  EXPECT_THROW({
    std::string line("1\t1\t.\tA\t.\t.\t.\tEND=2,3");
    VariantContext context = ParseVCFVariantLine(line, header_);
  }, file_parse_error);
}

TEST_F(VCFVariantParsingTest, ParsesUndefinedINFOfields) {
  using aseq::io::impl::ParseVCFVariantLine;
  using aseq::util::Attributes;

  EXPECT_NO_THROW({
    std::string line("1\t1\trs100\tA\tG\t100.0\tLowQual\tNEW_FLAG;NEW_STRINGS=1,2");
    VariantContext context = ParseVCFVariantLine(line, header_);
    ASSERT_TRUE(header_.HasINFOField("NEW_FLAG"));
    ASSERT_TRUE(header_.HasINFOField("NEW_STRINGS"));
  });
}

TEST_F(VCFVariantParsingTest, ParsesMinimalVariantWithSamples) {
  using aseq::io::impl::ParseVCFVariantLine;
  using aseq::util::Attributes;
  header_.SetSamples({"NA12878", "NA12891"});
  ASSERT_TRUE(header_.NumSamples() == 2);

  EXPECT_NO_THROW({
    std::string line("1\t1\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0|1");
    VariantContext context = ParseVCFVariantLine(line, header_);
    ASSERT_EQ(2, context.NumGenotypes());

    {
      auto& gt = context.GetGenotype("NA12878");
      ASSERT_EQ(&context, &gt.variant_context());
      EXPECT_FALSE(gt.HasAttribute(VCFHeader::FORMAT::GT));
      EXPECT_EQ(Genotype::kRefAlt, gt.alleles());
      EXPECT_FALSE(gt.Phased());
      EXPECT_TRUE(gt.Diploid());
    }

    {
      auto& gt = context.GetGenotype("NA12891");
      EXPECT_EQ(Genotype::kRefAltP, gt.alleles());
      EXPECT_TRUE(gt.Phased());
      EXPECT_TRUE(gt.Diploid());
    }

  });
}

TEST_F(VCFVariantParsingTest, ParsesComplexVariantWithSamples) {
  using aseq::io::impl::ParseVCFVariantLine;
  using aseq::util::Attributes;
  header_.SetSamples({"NA12878", "NA12891"});
  ASSERT_EQ(2, header_.NumSamples());

  EXPECT_NO_THROW({
    std::string line("1\t1\t.\tA\tG\t.\t.\t.\tGT:DP:FT:GQ:HQ\t0/1:20:dp;gq:30:4,5\t0|1:25:.:.:.,.");
    VariantContext context = ParseVCFVariantLine(line, header_);
    ASSERT_TRUE(context.NumGenotypes() == 2);

    {
      auto& gt = context.GetGenotype("NA12878");
      auto& ft = gt.GetAttribute<VariantContext::Filters>(VCFHeader::FORMAT::FT);
      EXPECT_EQ((VariantContext::Filters{"dp", "gq"}), ft);
    }

    {
      auto& gt = context.GetGenotype("NA12891");
      EXPECT_TRUE(gt.Phased());
      EXPECT_TRUE(gt.Diploid());
      EXPECT_EQ(25, gt.GetAttribute<Attributes::Integer>(VCFHeader::FORMAT::DP));
      EXPECT_FALSE(gt.HasAttribute(VCFHeader::FORMAT::FT));
      EXPECT_FALSE(gt.HasAttribute(VCFHeader::FORMAT::GQ));
      EXPECT_FALSE(gt.HasAttribute(VCFHeader::FORMAT::HQ));
    }

  });
}
