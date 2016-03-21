//
// Created by Michael Linderman on 12/16/15.
//

#include <glog/logging.h>
#include <boost/filesystem.hpp>

#include "aseq/util/exception.hpp"
#include "aseq/io/variant.hpp"
#include "aseq/io/vcf.hpp"

using namespace aseq::util;
namespace fs = boost::filesystem;

namespace aseq {
namespace io {

VariantSinkInterface::FactoryResult VariantSinkInterface::MakeVariantSink(FileFormat format,
                                                                          std::ostream& ostream) {
  switch (format) {
    default:
      LOG(INFO) << "UNKNOWN or unsupported file format specified, defaulting to VCF 4.2 output";
      return VariantSinkInterface::MakeVariantSink(FileFormat::VCF4_2, ostream);
    case FileFormat::VCF4_1:
    case FileFormat::VCF4_2:
      VCFHeader header(format);
      auto writer = ASCIILineWriterInterface::MakeLineWriter(ostream);
      return std::make_unique<VCFSink>(header, std::move(writer));
  }
}

VariantSinkInterface::FactoryResult VariantSinkInterface::MakeVariantSink(FileFormat format,
                                                                          const fs::path& path) {
  switch (format) {
    default:
      LOG(INFO) << "UNKNOWN or unsupported file format specified, defaulting to VCF 4.2 output";
      return VariantSinkInterface::MakeVariantSink(FileFormat::VCF4_2, path);
    case FileFormat::VCF4_1:
    case FileFormat::VCF4_2:
      VCFHeader header(format);
      auto writer = ASCIILineWriterInterface::MakeLineWriter(path, format);
      return std::make_unique<VCFSink>(header, std::move(writer));
  }
}

VariantSinkInterface::FactoryResult VariantSinkInterface::MakeVariantSink(
    const VariantSourceInterface& source, std::ostream& ostream) {
  auto writer = ASCIILineWriterInterface::MakeLineWriter(ostream);
  switch (source.file_format()) {
    default:
      throw file_write_error() << error_message("unsupported variant sink file format");
    case FileFormat::VCF4_1:
    case FileFormat::VCF4_2:
      const VCFSource& vcf_source = dynamic_cast<const VCFSource&>(source);
      return std::make_unique<VCFSink>(vcf_source.header(), std::move(writer));
  }
}

}  // io namespace
}  // aseq namespace