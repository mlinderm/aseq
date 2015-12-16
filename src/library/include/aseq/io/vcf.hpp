//
// Created by Michael Linderman on 12/13/15.
//

#pragma once

#include <type_traits>
#include <unordered_map>

#include "aseq/io/line.hpp"
#include "aseq/io/variant.hpp"

namespace aseq {
namespace io {

class VCFSource;
class VCFSink;

class VCFHeader {
 public:
  VCFHeader() : file_format_(FileFormat::UNKNOWN) {}
  VCFHeader(FileFormat type) : file_format_(type) {}

  FileFormat file_format() const { return file_format_; }
  void set_file_format(FileFormat format) { file_format_ = format; }

 private:
  FileFormat file_format_;

  friend class VCFSource;
};

namespace impl {

template <typename Line>
class VCFVariantLineParser;
}

class VCFSource : public VariantSourceInterface {
 public:
  typedef ASCIILineReaderInterface::FactoryResult Reader;

  VCFSource() = delete;
  VCFSource(FileFormat format, Reader&& reader);

  virtual FileFormat file_format() const override;
  virtual NextResult NextVariant() override;

  const VCFHeader& header() const { return header_; }

 private:
  VCFHeader header_;
  Reader reader_;
};

class VCFSink : public VariantSinkInterface {
 public:
  typedef ASCIILineWriterInterface::FactoryResult Writer;

 public:
  VCFSink() = delete;
  VCFSink(const VCFHeader& header, Writer&& writer);

  virtual void PushVariant(const model::VariantContext& context) override;

 private:
  VCFHeader header_;
  Writer writer_;
};
}
}