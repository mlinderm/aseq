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

namespace impl {

VariantSinkInterface::FactoryResult MakeVariantSink(
    FileFormat format, ASCIILineWriterInterface::FactoryResult&& writer) {
  switch (format) {
    default:
      LOG(INFO) << "UNKNOWN or unsupported file format specified, defaulting to VCF 4.2 output";
      return MakeVariantSink(FileFormat::VCF4_2, std::move(writer));
    case FileFormat::VCF4_1:
    case FileFormat::VCF4_2:
      VCFHeader header(format);
      return std::make_unique<VCFSink>(header, std::move(writer));
  }
}

VariantSinkInterface::FactoryResult MakeVariantSink(
    const VariantSourceInterface& source, ASCIILineWriterInterface::FactoryResult&& writer,
    bool sites_only = false) {
  switch (source.file_format()) {
    default:
      throw file_write_error() << error_message("unsupported variant sink file format");
    case FileFormat::VCF4_1:
    case FileFormat::VCF4_2:
      const VCFHeader& vcf_header = dynamic_cast<const VCFHeader&>(source.header());
      if (sites_only) {
        VCFHeader new_header = vcf_header;
        new_header.SetSitesOnly();
        return std::make_unique<VCFSink>(new_header, std::move(writer));
      } else
        return std::make_unique<VCFSink>(vcf_header, std::move(writer));
  }
}

}  // impl namespace

VariantSinkInterface::FactoryResult VariantSinkInterface::MakeVariantSink(FileFormat format,
                                                                          std::ostream& ostream) {
  auto writer = ASCIILineWriterInterface::MakeLineWriter(ostream);
  return impl::MakeVariantSink(format, std::move(writer));
}

VariantSinkInterface::FactoryResult VariantSinkInterface::MakeVariantSink(FileFormat format,
                                                                          const fs::path& path) {
  auto writer = ASCIILineWriterInterface::MakeLineWriter(path, format);
  return impl::MakeVariantSink(format, std::move(writer));
}

VariantSinkInterface::FactoryResult VariantSinkInterface::MakeVariantSink(
    const VariantSourceInterface& source, std::ostream& ostream, bool sites_only) {
  auto writer = ASCIILineWriterInterface::MakeLineWriter(ostream);
  return impl::MakeVariantSink(source, std::move(writer), sites_only);
}

VariantSinkInterface::FactoryResult VariantSinkInterface::MakeVariantSink(
    const VariantSourceInterface& source, const boost::filesystem::path& path, bool sites_only) {
  auto writer = ASCIILineWriterInterface::MakeLineWriter(path, source.file_format());
  return impl::MakeVariantSink(source, std::move(writer), sites_only);
}

}  // io namespace
}  // aseq namespace