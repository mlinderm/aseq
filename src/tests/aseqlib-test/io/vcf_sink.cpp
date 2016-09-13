//
// Created by Michael Linderman on 3/6/16.
//

#include <cppformat/format.h>
#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include "aseq/io/vcf.hpp"
#include "aseq/model/variant_context.hpp"

using namespace aseq::util;
using namespace aseq::io;
using namespace aseq::model;

namespace aseq {
namespace io {
namespace impl {
std::string GenerateVCFVariant(VCFHeader&, const VariantContext&);
}  // namespace impl
}  // namespace io
}  // namespace aseq

class VCFVariantGeneratingTest : public ::testing::Test {
 public:
  VCFVariantGeneratingTest() : header_(FileFormat::VCF4_2) {}

 protected:
  virtual void SetUp() {
    header_.AddINFOField(VCFHeader::INFO::AC);
    header_.AddFILTERField(VCFHeader::FILTER::PASS);
    header_.AddFORMATField(VCFHeader::FORMAT::GT);
    header_.AddFORMATField(VCFHeader::FORMAT::GQ);
  }

  VCFHeader header_;
};

TEST_F(VCFVariantGeneratingTest, GeneratesVCFHeaderWithSamples) {
  header_.SetSamples({"Sample0", "Sample1"});

  std::stringstream content;
  VCFSink sink(header_, ASCIILineWriterInterface::MakeLineWriter(content));
  auto source = VariantSourceInterface::MakeVariantSource(content);
  ASSERT_TRUE(source);

  EXPECT_EQ(FileFormat::VCF4_2, source->file_format());

  auto* vcf_source = dynamic_cast<VCFSource*>(source.get());
  auto& header = vcf_source->header();
  EXPECT_TRUE(header.HasINFOField(VCFHeader::INFO::AC));
  EXPECT_TRUE(header.HasFILTERField(VCFHeader::FILTER::PASS));
  EXPECT_TRUE(header.HasFORMATField(VCFHeader::FORMAT::GT));
  ASSERT_EQ(2, header.NumSamples());
  for (size_t i = 0; i < 2; i++) {
    EXPECT_EQ(Sample(fmt::format("Sample{}", i)), header.sample(i));
  }
}

TEST_F(VCFVariantGeneratingTest, GeneratesSitesOnlyVariants) {
  using aseq::io::impl::GenerateVCFVariant;
  EXPECT_NO_THROW({
    VariantContext cxt("1", 1, Allele::A, Allele::T);
    std::string line = GenerateVCFVariant(header_, cxt);
    EXPECT_EQ("1\t1\t.\tA\tT\t.\t.\t.", line);
  });

  EXPECT_NO_THROW({
    aseq::model::impl::VariantContextData data;
    data.contig_ = "1";
    data.pos_ = 1;
    data.ids_ = VariantContext::IDs({"rs100", "rs101"});
    data.ref_ = Allele::A;
    data.alts_ = VariantContext::Alleles({Allele::T, Allele::C});
    data.qual_ = 100.0;
    data.filters_ = VariantContext::Filters({VCFHeader::FILTER::PASS, "FAIL"});
    data.attrs_.emplace(VCFHeader::INFO::AC, Attributes::Integers({1}));

    header_.AddFILTERField(VCFHeader::Field("FAIL","Failed filter"));
    VariantContext cxt(std::move(data));
    std::string line = GenerateVCFVariant(header_, cxt);
    EXPECT_EQ("1\t1\trs100;rs101\tA\tT,C\t100.0\tPASS;FAIL\tAC=1", line);
  });

  EXPECT_NO_THROW({
    header_.AddINFOField(VCFHeader::INFO::CIEND);
    VariantContext cxt("1", 1, Allele::A, Allele::T);
    cxt.SetAttribute(VCFHeader::INFO::CIEND, Attributes::Integers({-5, 4}));
    std::string line = GenerateVCFVariant(header_, cxt);
    EXPECT_EQ("1\t1\t.\tA\tT\t.\t.\tCIEND=-5,4", line);
  });
}

TEST_F(VCFVariantGeneratingTest, GeneratesVariantsWithSamples) {
  header_.SetSamples({"Sample0", "Sample1"});

  using aseq::io::impl::GenerateVCFVariant;
  EXPECT_NO_THROW({
    VariantContext cxt("1", 1, Allele::A, Allele::T);
    cxt.AddGenotype("Sample0", Genotype::kRefAlt,
                    Attributes({{VCFHeader::FORMAT::GQ, Attributes::mapped_type(30)}}));
    cxt.AddGenotype("Sample1", Genotype::kRefAltP, Attributes());
    std::string line = GenerateVCFVariant(header_, cxt);
    EXPECT_EQ("1\t1\t.\tA\tT\t.\t.\t.\tGT:GQ\t0/1:30\t0|1:.", line);
  });
}

TEST_F(VCFVariantGeneratingTest, GeneratesVCFWithSamples) {
  header_.SetSamples({"Sample0", "Sample1"});
  std::stringstream content;
  ASSERT_NO_THROW({
    VCFSink sink(header_, ASCIILineWriterInterface::MakeLineWriter(content));
    {
      VariantContext cxt("1", 1, Allele::A, Allele::T);
      cxt.AddGenotype("Sample0", Genotype::kRefAlt, Attributes());
      cxt.AddGenotype("Sample1", Genotype::kRefAltP, Attributes());
      sink.PushVariant(cxt);
    }
  });

  std::string line;
  ASSERT_TRUE(std::getline(content, line));
  EXPECT_EQ("##fileformat=VCFv4.2", line);
  ASSERT_TRUE(std::getline(content, line));
  EXPECT_EQ(
      "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT "
      "allele, in the same order as listed\">",
      line);
  ASSERT_TRUE(std::getline(content, line));
  EXPECT_EQ("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", line);
  ASSERT_TRUE(std::getline(content, line));
  EXPECT_EQ("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">", line);
  ASSERT_TRUE(std::getline(content, line));
  EXPECT_EQ("##FILTER=<ID=PASS,Description=\"PASSing\">", line);
  ASSERT_TRUE(std::getline(content, line));
  EXPECT_EQ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample0\tSample1", line);
  ASSERT_TRUE(std::getline(content, line));
  EXPECT_EQ("1\t1\t.\tA\tT\t.\t.\t.\tGT:GQ\t0/1:.\t0|1:.", line);
  EXPECT_FALSE(std::getline(content, line));
}