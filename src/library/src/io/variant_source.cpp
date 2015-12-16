//
// Created by Michael Linderman on 12/13/15.
//

#include <iostream>

#include <boost/filesystem.hpp>
#include <aseq/io/vcf.hpp>

#include "aseq/util/exception.hpp"
#include "aseq/io/variant.hpp"

using namespace aseq::util;

namespace fs = boost::filesystem;

namespace aseq {
namespace io {

namespace {

static const std::map<std::string, FileFormat> kHeader2Format = {
    {"##fileformat=VCFv4.1", FileFormat::VCF4_1}, {"##fileformat=VCFv4.2", FileFormat::VCF4_2}};

FileFormat DetectSourceType(const ASCIILineReaderInterface::FactoryResult& reader) {
  auto line = reader->ReadNextLine();
  if (line) {
    auto f = kHeader2Format.find(boost::copy_range<std::string>(*line));
    return (f != kHeader2Format.end()) ? f->second : FileFormat::UNKNOWN;
  } else
    return FileFormat::UNKNOWN;
}

VariantSourceInterface::FactoryResult VariantSourcePicker(
    ASCIILineReaderInterface::FactoryResult reader) {
  FileFormat type = DetectSourceType(reader);
  switch (type) {
    default:
      throw file_parse_error() << error_message("unable to determine or unsupported file type");
    case FileFormat::VCF4_1:
    case FileFormat::VCF4_2:
      return std::make_unique<VCFSource>(type, std::move(reader));
  }
}

}  // anonymous namespace

VariantSourceInterface::FactoryResult VariantSourceInterface::MakeVariantSource(
    std::istream& istream) {
  auto reader = ASCIILineReaderInterface::MakeLineReader(istream);
  return VariantSourcePicker(std::move(reader));
}

VariantSourceInterface::FactoryResult VariantSourceInterface::MakeVariantSource(
    const boost::filesystem::path& path) {
  if (path == "-") return MakeVariantSource(std::cin);
  return VariantSourcePicker(ASCIILineReaderInterface::MakeLineReader(path));
}

}  // io namespace
}  // aseq namespac