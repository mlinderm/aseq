//
// Created by Michael Linderman on 12/13/15.
//

#pragma once

#include <type_traits>
#include <unordered_map>
#include <iosfwd>

#include <boost/range/adaptor/map.hpp>

#include "aseq/io/line.hpp"
#include "aseq/io/variant.hpp"
#include "aseq/model/genotype.hpp"

namespace aseq {
namespace io {

class VCFSource;
class VCFSink;

class VCFHeader : public VariantHeaderInterface {
 public:
  class Field {
   public:
    typedef util::Attributes::key_type ID;
    enum class Type { FLAG, INTEGER, FLOAT, CHARACTER, STRING, GENOTYPE, FILTER };
    typedef int Number;

    static const Number R = -4, A = -3, G = -2, UNBOUNDED = -1;

    Field() : nmbr_(0), type_(Type::FLAG) {}
    Field(const ID& id, const std::string& desc = "") : Field(id, 0, Type::FLAG, desc) {}
    Field(const ID& id, Number number, Type type, const std::string& desc)
        : id_(id), nmbr_(number), type_(type), desc_(desc) {}

    operator const ID&() const { return id_; }

    bool IsFlag() const { return nmbr_ == 0 && type_ == Type::FLAG; }
    bool IsScalar() const { return nmbr_ == 1; }

    ID id_;
    Number nmbr_;
    Type type_;
    std::string desc_;
  };

  struct FILTER {
#define FILTER_FIELD(id, description) static const Field id;
#include "vcf_fields.def"
#undef FILTER_FIELD
  };

  struct INFO {
#define INFO_FIELD(id, number, type, description) static const Field id;
#include "vcf_fields.def"
#undef INFO_FIELD
  };

  struct FORMAT {
#define FORMAT_FIELD(id, number, type, description) static const Field id;
#include "vcf_fields.def"
#undef FORMAT_FIELD
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

#define FIELD(prefix)                                                    \
  Fields& prefix() { return prefix##_; }                                 \
  boost::select_second_const_range<Fields> prefix##Values() const {      \
    return boost::adaptors::values(prefix##_);                           \
  }                                                                      \
  std::pair<const Field&, bool> Add##prefix##Field(const Field& field) { \
    return AddField(prefix##_, field);                                   \
  };                                                                     \
  bool Has##prefix##Field(const Fields::key_type& key) const { return HasField(prefix##_, key); }

  FIELD(FILTER)
  FIELD(INFO)
  FIELD(FORMAT)

#undef FIELD

  // Sample methods
  size_t NumSamples() const override { return samples_.size(); }
  const model::Sample& Sample(size_t idx) const override { return samples_.at(idx); }
  void SetSamples(std::initializer_list<model::Sample> samples) { samples_.assign(samples); }

 private:
  FileFormat file_format_;
  Fields FILTER_, INFO_, FORMAT_;
  std::vector<model::Sample> samples_;

  friend class VCFSource;
  friend class VCFSink;
};

std::ostream& operator<<(std::ostream&, const VCFHeader::Field&);

namespace impl {

template <typename Line>
class VCFVariantParser;

class VCFVariantGeneratorInterface {

};
}  // namespace impl

class VCFSource : public VariantSourceInterface {
 public:
  typedef ASCIILineReaderInterface::FactoryResult Reader;

  VCFSource() = delete;
  VCFSource(FileFormat format, Reader&& reader);

  virtual FileFormat file_format() const override;
  const VCFHeader& header() const override { return header_; }

  virtual NextResult NextVariant() override;

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

  // Specify the deleter to enable the incomplete type
  std::unique_ptr<impl::VCFVariantGeneratorInterface> generator_;
};
}
}