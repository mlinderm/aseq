//
// Created by Michael Linderman on 12/13/15.
//

#pragma once

#include <type_traits>
#include <unordered_map>

#include <boost/range/adaptor/map.hpp>

#include "aseq/io/line.hpp"
#include "aseq/io/variant.hpp"

namespace aseq {
namespace io {

class VCFSource;
class VCFSink;

class VCFHeader {
 public:
  class Field {
   public:
    typedef util::Attributes::key_type ID;
    enum class Type { FLAG, INTEGER, FLOAT, CHARACTER, STRING, GENOTYPE, FILTER };
    typedef int Number;

    static const Number R = -4, A = -3, G = -2, UNBOUNDED = -1;

    Field() : number_(0), type_(Type::FLAG) {}
    Field(const ID& id, const std::string& description = "")
        : id_(id), number_(0), type_(Type::FLAG), description_(description) {}
    Field(const ID& id, Number number, Type type, const std::string& description)
        : id_(id), number_(number), type_(type), description_(description) {}

    operator const ID&() const { return id_; }

    bool IsFlag() const { return number_ == 0 && type_ == Type::FLAG; }
    bool IsScalar() const { return number_ == 1; }

    ID id_;
    Number number_;
    Type type_;
    std::string description_;
  };

  struct INFO {
#define INFO_FIELD(id, number, type, description) static const Field id;
#include "vcf_fields.def"
#undef INFO_FIELD
  };

  struct FILTER {
#define FILTER_FIELD(id, description) static const Field id;
#include "vcf_fields.def"
#undef FILTER_FIELD
  };

  typedef std::unordered_map<util::Attributes::key_type, Field> Fields;

 public:
  VCFHeader() : file_format_(FileFormat::UNKNOWN) {}
  VCFHeader(FileFormat type) : file_format_(type) {}

  FileFormat file_format() const { return file_format_; }
  void set_file_format(const FileFormat& format) { file_format_ = format; }

  // Field methods
  static std::pair<const Field&, bool> AddField(Fields&, const Field&);
  static bool HasField(const Fields&, const Fields::key_type&);

#define FIELD(prefix)                                                   \
  Fields& prefix() { return prefix##_; }                                \
  boost::select_second_const_range<Fields> prefix##_Values() const {    \
    return boost::adaptors::values(prefix##_);                          \
  }                                                                     \
  std::pair<const Field&, bool> prefix##_AddField(const Field& field) { \
    return AddField(prefix##_, field);                                  \
  };                                                                    \
  bool prefix##_HasField(const Fields::key_type& key) const { return HasField(prefix##_, key); }

  FIELD(INFO)
  FIELD(FILTER)

#undef FIELD

 private:
  FileFormat file_format_;
  Fields INFO_, FILTER_;

  friend class VCFSource;
};

namespace impl {

template <typename Line>
class VCFVariantParser;
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

  typedef impl::VCFVariantParser<Line> Parser;
  struct ParserDeleter {
    void operator()(Parser*);
  };

  // Specify the deleter to enable the incomplete type
  std::unique_ptr<Parser, ParserDeleter> parser_;
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