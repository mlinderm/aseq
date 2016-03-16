//
// Created by Michael Linderman on 12/16/15.
//
#define BOOST_SPIRIT_USE_PHOENIX_V3 1
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/repository/include/karma_confix.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/tuple/elem.hpp>

#include <glog/logging.h>
#include <cppformat/format.h>

#include "aseq/io/vcf.hpp"
#include "vcf_io.def"

namespace km = boost::spirit::karma;
namespace phx = boost::phoenix;
namespace repo = boost::spirit::repository;
namespace fusion = boost::fusion;

// Adapt structs for use with boost spirit
// clang-format off
BOOST_FUSION_ADAPT_STRUCT(
aseq::io::VCFHeader::Field,
(aseq::util::Attributes::key_type, id_)
(aseq::io::VCFHeader::Field::Number, nmbr_)
(aseq::io::VCFHeader::Field::Type, type_)
(std::string, desc_)
)

#define ENUM_THEN_STRING(r, data, elem) (BOOST_PP_TUPLE_ELEM(2, 0, elem), BOOST_PP_TUPLE_ELEM(2, 1,elem))

typedef fusion::vector<
aseq::model::Contig const &,
aseq::model::Pos,
aseq::model::VariantContext::IDs const &,
aseq::model::Allele const &,
aseq::model::VariantContext::Alleles const &,
aseq::model::VariantContext::Qual const &,
aseq::model::VariantContext::Filters const &,
aseq::util::Attributes const &
> VariantContext2VCF;
// clang-format on

namespace boost {
namespace spirit {
namespace traits {

namespace {

template <typename Transformed>
struct transform_hold_any {
  typedef Transformed const &type;
  static type pre(aseq::util::Attributes::mapped_type const &val) {
    return aseq::util::any_cast<Transformed>(val);
  }
};
}

#define TRANSFORM(type)                                                                   \
  template <>                                                                             \
  struct transform_attribute<aseq::util::Attributes::mapped_type const, type, km::domain> \
      : transform_hold_any<type> {};

TRANSFORM(aseq::util::Attributes::Integer)
TRANSFORM(aseq::util::Attributes::Integers)
TRANSFORM(aseq::util::Attributes::Float)
TRANSFORM(aseq::util::Attributes::Floats)
TRANSFORM(aseq::util::Attributes::Character)
TRANSFORM(aseq::util::Attributes::Characters)
TRANSFORM(aseq::util::Attributes::Strings)

#undef TRANSFORM

template <>
struct transform_attribute<aseq::model::VariantContext const, VariantContext2VCF, km::domain> {
  typedef VariantContext2VCF type;

  static type pre(aseq::model::VariantContext const &val) {
    return type(val.contig(), val.pos(), val.ids(), val.ref(), val.alts(), val.qual(),
                val.filters(), val.attributes());
  }
};
}  // namespace traits
}  // namespace spirit
}  // namespace boost

namespace aseq {
namespace io {
namespace impl {

template <typename Iterator>
struct VCFHeaderGenerator {
  VCFHeaderGenerator() {
    using km::lit;
    using namespace km::labels;

    // clang-format off
    file_formats_.add BOOST_PP_SEQ_FOR_EACH(ENUM_THEN_STRING, BOOST_PP_NIL, FILE_TYPES);
    field_number_.add BOOST_PP_SEQ_FOR_EACH(ENUM_THEN_STRING, BOOST_PP_NIL, FIELD_NUMBER);
    field_types_.add
      BOOST_PP_SEQ_FOR_EACH(ENUM_THEN_STRING, BOOST_PP_NIL, FIELD_TYPES)
      (VCFHeader::Field::Type::GENOTYPE, "String")
      ;

    format_ = lit("##fileformat=") << file_formats_ << km::eol;
    full_field_ = repo::confix('<','>')[(
      lit("ID=")           << km::string <<
      lit(",Number=")      << (field_number_ | km::int_) <<
      lit(",Type=")        << field_types_ <<
      lit(",Description=") << repo::confix('"','"')[km::string]
    )] << km::eol;
    id_and_desc_field_ = repo::confix('<','>')[(
      lit("ID=") << km::string <<
      km::skip[field_number_] << km::skip[field_types_] <<
      lit(",Description=") << repo::confix('"','"')[km::string]
    )] << km::eol;
    required_columns_ %= lit("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") << (
      km::eps(phx::size(_val) == 0)
      | (lit("\tFORMAT\t") << (km::stream % '\t'))
    );
    // clang-format on
  }

  km::rule<Iterator, FileFormat> format_;
  km::rule<Iterator, VCFHeader::Field> full_field_, id_and_desc_field_;
  km::rule<Iterator, std::vector<model::Sample> > required_columns_;

  km::symbols<FileFormat, const char *> file_formats_;
  km::symbols<VCFHeader::Field::Type, const char *> field_types_;
  km::symbols<VCFHeader::Field::Number, const char *> field_number_;
};

template <typename Iterator>
class VCFVariantGenerator : public VCFVariantGeneratorInterface {
 public:
  VCFVariantGenerator() = delete;
  VCFVariantGenerator(VCFHeader &header) : header_(header) {
    using namespace km::labels;
    using util::Attributes;

    for (auto f : header_.INFOValues()) {
      info_keys_.emplace(f, GetGenerator(f));
    }

    // clang-format off
    auto sep = km::lit('\t');
    auto missing = km::lit('.');

    variant_  = vcf_line_(_val);
    vcf_line_ =
      chrom_                 << sep <<
      pos_                   << sep <<
      (id_        | missing) << sep <<
      allele_                << sep <<
      alt_(_r1)              << sep <<
      (km::float_ | missing) << sep <<
      (filter_    | missing) << sep <<
      (info_      | missing);

    chrom_  = km::string;
    pos_    = km::int_generator<model::Pos,10,false>();
    alt_    %= (
      (allele_ % ',')
      | missing
    );

    id_     = km::string % ';';
    filter_ = km::string % ';';
    info_   = info_entry_ % ';';

    allele_ = km::string;

    auto get_info_rule = phx::bind(&VCFVariantGenerator::KeyToRule, this, phx::ref(info_keys_), phx::at_c<0>(_val));
    info_entry_ %= km::eps[_a = get_info_rule] <<
            km::string << (km::eps(!_a) | (km::lit('=') << km::lazy(*_a)));

    integer_value_    = km::attr_cast<Attributes::mapped_type const &, Attributes::Integer>(km::int_);
    integers_value_   = km::attr_cast<Attributes::mapped_type const &, Attributes::Integers>(km::int_ % ',');
    float_value_      = km::attr_cast<Attributes::mapped_type const &, Attributes::Float>(km::float_);
    floats_value_     = km::attr_cast<Attributes::mapped_type const &, Attributes::Floats>(km::float_ % ',');
    character_value_  = km::attr_cast<Attributes::mapped_type const &, Attributes::Character>(km::char_);
    characters_value_ = km::attr_cast<Attributes::mapped_type const &, Attributes::Characters>(km::char_ % ',');
    string_value_     = km::stream;
    strings_value_    = km::attr_cast<Attributes::mapped_type const &, Attributes::Strings>(km::string % ',');
    // clang-format on
  }

  bool Generate(Iterator &itr, const model::VariantContext &cxt) {
    return km::generate(itr, variant_, cxt);
  }

 private:
  typedef km::rule<Iterator, util::Attributes::mapped_type const &()> AttributeRule;
  typedef std::unordered_map<util::Attributes::key_type, AttributeRule const *> KeyToRuleMap;

  const AttributeRule *KeyToRule(KeyToRuleMap &map, const util::Attributes::key_type &key) {
    auto i = map.find(key);
    if (i == map.end()) {
      LOG(INFO) << fmt::format("Treating attribute {} as STRING", key);
      std::tie(i, std::ignore) = map.emplace(key, &string_value_);
    }
    return i->second;
  }

  const AttributeRule *GetGenerator(const VCFHeader::Field &field) {
#define TORULE(TYPE, SCALAR, VECTOR) \
  case VCFHeader::Field::Type::TYPE: \
    return (field.IsScalar()) ? SCALAR : VECTOR

    switch (field.type_) {
      default:
        throw util::file_write_error()
            << util::error_message(fmt::format("unsupported attribute type for {}", field));
      case VCFHeader::Field::Type::FLAG:
        return nullptr;
        TORULE(INTEGER, &integer_value_, &integers_value_);
        TORULE(FLOAT, &float_value_, &floats_value_);
        TORULE(CHARACTER, &character_value_, &characters_value_);
        TORULE(STRING, &string_value_, &strings_value_);
    }
#undef TORULE
  }

  km::rule<Iterator, model::VariantContext const &()> variant_;
  km::rule<Iterator, VariantContext2VCF(const model::VariantContext &)> vcf_line_;

  km::rule<Iterator, model::Contig> chrom_;
  km::rule<Iterator, model::Pos()> pos_;
  km::rule<Iterator, model::VariantContext::Alleles const &(const model::VariantContext &)> alt_;
  km::rule<Iterator, model::Allele> allele_;

  km::rule<Iterator, model::VariantContext::IDs const &()> id_;
  km::rule<Iterator, model::VariantContext::Filters const &()> filter_;
  km::rule<Iterator, util::Attributes const &()> info_;
  km::rule<Iterator, util::Attributes::value_type const &(), km::locals<AttributeRule const *> >
      info_entry_;

  AttributeRule integer_value_, integers_value_, float_value_, floats_value_, character_value_,
      characters_value_, string_value_, strings_value_;

  KeyToRuleMap info_keys_;
  VCFHeader &header_;
};

typedef VCFVariantGenerator<std::back_insert_iterator<std::string> > Generator;

std::string GenerateVCFVariant(VCFHeader &header, const model::VariantContext &context) {
  std::string line;
  auto itr = std::back_inserter(line);

  VCFVariantGenerator<decltype(itr)> gen(header);
  gen.Generate(itr, context);

  return line;
}

}  // namespace impl

VCFSink::VCFSink(const VCFHeader &header, Writer &&writer)
    : header_(header), writer_(std::move(writer)) {
  std::string line;
  auto itr = std::back_inserter(line);

  impl::VCFHeaderGenerator<decltype(itr)> gen;
  km::generate(itr, gen.format_, header_.file_format_);

#define FIELDS(KIND, GENERATOR)                                 \
  for (auto &f : header_.KIND##Values()) {                      \
    km::generate(itr, km::lit("##" #KIND "=") << GENERATOR, f); \
  }
  FIELDS(INFO, gen.full_field_);
  FIELDS(FORMAT, gen.full_field_);
  FIELDS(FILTER, gen.id_and_desc_field_);
#undef FIELDS

  km::generate(itr, gen.required_columns_, header_.samples_);

  writer_->WriteLine(line);

  // Create the variant generator for subsequent use
  generator_ = std::make_unique<impl::Generator>(header_);
}

void VCFSink::PushVariant(const model::VariantContext &cxt) {
  auto &gen = static_cast<impl::Generator &>(*generator_);

  std::string line;
  auto itr = std::back_inserter(line);
  gen.Generate(itr, cxt);
  writer_->WriteLine(line);
}

}  // namespace io
}  // namespace aseq