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
bool ParseVCFHeaderLine(const std::string& line, VCFHeader& header);
VariantContext ParseVCFVariant(const std::string& line, VCFHeader& header);
}  // namespace impl
}  // namespace io
}  // namespace aseq

TEST(VCFHeaderParsingTest, ParsesMinimalVCFHeader) {
  std::stringstream content("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

  VCFSource source(FileFormat::VCF4_2, ASCIILineReaderInterface::MakeLineReader(content));
  EXPECT_EQ(FileFormat::VCF4_2, source.file_format());

  auto& header = source.header();
  EXPECT_TRUE(header.FILTER_HasField("PASS"));
  // EXPECT_TRUE(header.HasField(VCFHeader::KeyType::FORMAT, "GT"));
  // EXPECT_EQ(0, header.NumSamples());
}

class VCFVariantParsingTest : public ::testing::Test {
 public:
  VCFVariantParsingTest() : header_(FileFormat::VCF4_2) {}

 protected:
  virtual void SetUp() {
    header_.INFO_AddField(VCFHeader::INFO::END);
    header_.INFO_AddField(VCFHeader::INFO::SVTYPE);
    header_.INFO_AddField(VCFHeader::INFO::SVLEN);
    header_.INFO_AddField(VCFHeader::INFO::CIPOS);
    header_.INFO_AddField(VCFHeader::INFO::CIEND);
    header_.INFO_AddField(VCFHeader::INFO::HOMLEN);
    header_.INFO_AddField(VCFHeader::INFO::HOMSEQ);
  }

  VCFHeader header_;
};

TEST_F(VCFVariantParsingTest, ParsesMinimalSitesOnlyVariant) {
  using aseq::io::impl::ParseVCFVariant;

  std::string line = "1\t1\t.\tA\t.\t.\t.\t.";
  VariantContext context = ParseVCFVariant(line, header_);
  EXPECT_EQ(Contig("1"), context.contig());
  EXPECT_EQ(1, context.pos());
  EXPECT_EQ(1, context.end());
  EXPECT_EQ(Allele::A, context.ref());
}

TEST_F(VCFVariantParsingTest, ParsesComplexSitesOnlyVariant) {
  using aseq::io::impl::ParseVCFVariant;
  using aseq::util::Attributes;

  EXPECT_NO_THROW({
    std::string line("1\t1\trs100\tA\tG,<NON_REF>\t100.0\tLowQual\tEND=101");
    VariantContext context = ParseVCFVariant(line, header_);
    EXPECT_EQ(101, context.GetAttribute<Attributes::Integer>("END"));
  });
  EXPECT_NO_THROW({
    std::string line(
        "1\t2827694\trs2376870\tCGTGGATGCGGGGAC\tC\t.\tPASS\tSVTYPE=DEL;END=2827762;HOMLEN=1;"
        "HOMSEQ=G;SVLEN=-68");
    VariantContext context = ParseVCFVariant(line, header_);
  });
  EXPECT_NO_THROW({
    std::string line(
        "2\t321682\t.\tT\t<DEL>\t6\tPASS\tSVTYPE=DEL;END=321887;SVLEN=-205;CIPOS=-56,20;CIEND=-10,"
        "62");
    VariantContext context = ParseVCFVariant(line, header_);
  });

  EXPECT_NO_THROW({
    std::string line(
        "2\t14477084\t.\tC\t<DEL:ME:ALU>\t6\tPASS\tSVTYPE=DEL;END=14477381;SVLEN=-297;CIPOS=-22,18;"
        "CIEND=-12,32");
    VariantContext context = ParseVCFVariant(line, header_);
  });

  EXPECT_NO_THROW({
    std::string line("2\t321681\tbnd_W\tG\tG]17:198982]\t6\tPASS\tSVTYPE=BND");
    VariantContext context = ParseVCFVariant(line, header_);
  });
}

TEST_F(VCFVariantParsingTest, ParsesUndefinedINFOfields) {
  using aseq::io::impl::ParseVCFVariant;
  using aseq::util::Attributes;

  EXPECT_NO_THROW({
    std::string line("1\t1\trs100\tA\tG\t100.0\tLowQual\tNEW_FLAG;NEW_STRINGS=1,2");
    VariantContext context = ParseVCFVariant(line, header_);
    ASSERT_TRUE(header_.INFO_HasField("NEW_FLAG"));
    ASSERT_TRUE(header_.INFO_HasField("NEW_STRINGS"));
  });
}
