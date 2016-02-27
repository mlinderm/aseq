//
// Created by Michael Linderman on 12/13/15.
//

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <glog/logging.h>
#include <cppformat/format.h>

#include "aseq/util/exception.hpp"
#include "aseq/io/vcf.hpp"
#include "aseq/model/variant_context.hpp"

namespace x3 = boost::spirit::x3;

// Adapt structs for use with boost spirit
// clang-format off
BOOST_FUSION_ADAPT_STRUCT(
aseq::io::VCFHeader::Field,
(aseq::util::Attributes::key_type, id_)
(aseq::io::VCFHeader::Field::Number, number_)
(aseq::io::VCFHeader::Field::Type, type_)
(std::string, description_)
)
// clang-format on

namespace boost {
namespace algorithm {

namespace detail {
template <typename CharT>
struct is_charF : public predicate_facade<is_charF<CharT> > {
  typedef bool result_type;

  is_charF(CharT match) : match_(match) {}

  template <typename Char2T>
  bool operator()(Char2T value) const {
    return value == match_;
  }

  CharT match_;
};
}  // namespace detail

template <typename CharT>
inline detail::is_charF<CharT> is_char(CharT match) {
  return detail::is_charF<CharT>(match);
}

}  // namespace algorithm

// Pull names to the boost namespace
using algorithm::is_char;
}  // namespace boost

namespace boost {
namespace spirit {
namespace x3 {
namespace traits {
template <>
void move_to<const char*, aseq::util::Attributes::key_type>(const char* first, const char* last,
                                                            aseq::util::Attributes::key_type& key) {
  key = aseq::util::Attributes::key_type(first, last);
};
}  // namespace traits
}  // namespace x3
}  // namespace spirit
}  // namespace boost

namespace aseq {
namespace io {
namespace impl {

using model::VariantContext;

// Helper values and macros for parsing VCF fields
const auto kFieldSplitter = boost::is_char('\t');
const auto kSemicolonSplitter = boost::is_char(';');
const auto kCommaSplitter = boost::is_char(',');
const auto kFieldFinder = boost::token_finder(kFieldSplitter);
const auto kInfoFinder = boost::token_finder(kSemicolonSplitter);
const auto kFormatFinder = boost::token_finder(boost::is_char(':'));
const auto kEqualsFinder = boost::token_finder(boost::is_char('='));

template <typename Range>
bool IsDefined(const Range& range) {
  return range.size() != 1 || range.front() != '.';
}

template <typename Iterator>
bool IsDefined(const Iterator& begin, const Iterator& end) {
  return IsDefined(boost::make_iterator_range(begin, end));
}

template <typename Sequence, typename Range, typename Predicate>
Sequence& SplitOptionalRange(Sequence& container, const Range& input, Predicate predicate) {
  if (!IsDefined(input)) {
    return container;
  }
  return boost::split(container, input, predicate);
}

// Boost Spirit X3 infrastructure
namespace parser {

// Parser rules
using util::Attributes;
using model::VariantContext;
using x3::lit;
using x3::char_;

x3::rule<class header_line> const header_line = "VCF header line";
x3::rule<class header_field> const header_field = "VCF header field";
x3::rule<class field_attributes, VCFHeader::Field> const field_attributes =
    "VCF header field attributes";
x3::rule<class filter_or_alt_attributes, VCFHeader::Field> const filter_or_alt_attributes =
    "VCF header field with ID and description attributes";

struct header_tag {};
auto SetFormat = [&](auto& ctx) {
  VCFHeader& header = x3::get<header_tag>(ctx);
  header.set_file_format(x3::_attr(ctx));
};

#define FIELD_ADD(kind)                           \
  auto Add##kind##Field = [&](auto& ctx) {        \
    VCFHeader& header = x3::get<header_tag>(ctx); \
    header.kind##_AddField(x3::_attr(ctx));       \
  };

FIELD_ADD(INFO)
FIELD_ADD(FILTER)

#undef FIELD_ADD

const x3::symbols<FileFormat> kFileFormats{{"VCFv4.1", FileFormat::VCF4_1},
                                           {"VCFv4.2", FileFormat::VCF4_2}};

const x3::symbols<VCFHeader::Field::Type> kFieldTypes{
    {"String", VCFHeader::Field::Type::STRING},
    {"Character", VCFHeader::Field::Type::CHARACTER},
    {"Integer", VCFHeader::Field::Type::INTEGER},
    {"Float", VCFHeader::Field::Type::FLOAT},
    {"Flag", VCFHeader::Field::Type::FLAG}};

const x3::symbols<VCFHeader::Field::Number> kFieldNumbers{{"A", VCFHeader::Field::A},
                                                          {"G", VCFHeader::Field::G},
                                                          {"R", VCFHeader::Field::R},
                                                          {".", VCFHeader::Field::UNBOUNDED}};
// clang-format off
auto const id   = lit("ID=")          > x3::raw[+(char_ - ',')];
auto const nmbr = lit("Number=")      > (kFieldNumbers | x3::int_) ;
auto const type = lit("Type=")        > kFieldTypes;
auto const desc = lit("Description=") > x3::confix('"','"')[*(char_ - '"')];

auto const field_attributes_def =
    lit('=') > x3::confix('<','>')[id > ',' > nmbr > ',' > type > ',' > desc];
auto const filter_or_alt_attributes_def =
    lit('=') > x3::confix('<','>')[id > x3::attr(0) > x3::attr(VCFHeader::Field::Type::FLAG) > ',' > desc];

auto const header_field_def = (
    (lit("INFO") > field_attributes)[AddINFOField]
    | (lit("FILTER") > filter_or_alt_attributes)[AddFILTERField]
    | ((+x3::alnum) > '=' > (+(char_ - x3::eol)))
);

auto const header_line_def = lit("##") > (
    (lit("fileformat=") > kFileFormats[SetFormat])
    | header_field
) > -x3::eol;
// clang-format on

x3::rule<class ref_allele, VariantContext::Pos> const pos = "POS";
x3::rule<class qual, boost::optional<VariantContext::Qual> > const qual = "QUAL";

// clang-format off
auto const pos_def  = x3::long_long;
auto const qual_def = lit('.') | x3::float_;
auto const missing  = lit('.') % ',';
// clang-format on

BOOST_SPIRIT_DEFINE(header_line, header_field, field_attributes, filter_or_alt_attributes);
BOOST_SPIRIT_DEFINE(pos, qual);

template <typename Iterator>
struct AttributeParser {
  Attributes::key_type key_;

  AttributeParser() = delete;
  AttributeParser(const Attributes::key_type& key) : key_(key) {}

  virtual bool operator()(Iterator, const Iterator&, Attributes::mapped_type&) const {
    return true;
  }

  void operator()(Iterator begin, const Iterator& end, Attributes& attr) const {
    Attributes::mapped_type value;
    if (!this->operator()(begin, end, value)) {
      throw util::file_parse_error()
          << util::error_message(fmt::format("failed to parse value for {}", key_));
    }
    if (value.empty()) {
      return;
    }
    if (!attr.emplace(key_, std::move(value)).second) {
      throw util::file_parse_error()
          << util::error_message(fmt::format("attribute {} already exists", key_));
    }
  }
};

template <typename Iterator>
struct FlagAttributeParser : public AttributeParser<Iterator> {
  FlagAttributeParser(const Attributes::key_type key) : AttributeParser<Iterator>(key) {}
  bool operator()(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const {
    value = true;
    return begin == end;
  }
};

template <typename Iterator>
struct StringAttributeParser : public AttributeParser<Iterator> {
  StringAttributeParser(const Attributes::key_type key) : AttributeParser<Iterator>(key) {}
  bool operator()(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const {
    if (IsDefined(begin, end)) {
      value = std::move(Attributes::String(begin, end));
    }
    return true;
  }
};

template <typename Iterator>
struct StringsAttributeParser : public AttributeParser<Iterator> {
  StringsAttributeParser(const Attributes::key_type key) : AttributeParser<Iterator>(key) {}
  bool operator()(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const {
    Attributes::Strings strings;
    SplitOptionalRange(strings, boost::make_iterator_range(begin, end), kCommaSplitter);
    if (!strings.empty()) {
      value = std::move(strings);
    }
    return true;
  }
};

#define PARSER(NAME, TYPE, EXPR)                                                                 \
  template <typename Iterator>                                                                   \
  struct NAME : public AttributeParser<Iterator> {                                               \
    NAME(const Attributes::key_type key) : AttributeParser<Iterator>(key) {}                     \
    bool operator()(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const { \
      TYPE actual_value;                                                                         \
      bool r = x3::parse(begin, end, EXPR, actual_value) && begin == end;                        \
      if (r) {                                                                                   \
        value = std::move(actual_value);                                                         \
      }                                                                                          \
      return r;                                                                                  \
    }                                                                                            \
  };

PARSER(IntegerAttributeParser, Attributes::Integer, missing | x3::int_);
PARSER(IntegersAttributeParser, Attributes::Integers, missing | (x3::int_ % ','));
PARSER(FloatAttributeParser, Attributes::Float, missing | x3::float_);
PARSER(FloatsAttributeParser, Attributes::Floats, missing | (x3::float_ % ','));
PARSER(CharacterAttributeParser, Attributes::Character, missing | x3::char_);
PARSER(CharactersAttributeParser, Attributes::Characters, missing | (x3::char_ % ','));

#undef PARSER

template <typename Range, typename Parser, typename Attribute>
void ParseRange(const Range& range, const Parser& parser, Attribute& attr) {
  typename Range::iterator begin = range.begin();
  if (!x3::parse(begin, range.end(), parser, attr) || begin != range.end()) {
    throw util::file_parse_error()
        << util::error_message(fmt::format("failed to parse {}", parser.name));
  }
};

}  // namespace parser

template <typename Line>
class VCFVariantParser {
  typedef typename Line::const_iterator Iterator;
  typedef parser::AttributeParser<Iterator> AttributeParser;
  using Attributes = util::Attributes;

 public:
  VCFVariantParser(VCFHeader& header) : header_(header) {
    // "Missing" INFO field
    info_keys_.emplace(".", std::make_unique<AttributeParser>("."));
    for (auto f : header.INFO_Values()) {
      info_keys_.emplace(f, GetParser(f));
    }
  }

  VariantContext ParseVCFVariant(const Line& line) {
    using util::file_parse_error;
    using util::error_message;

    // Split VCF line into component fields
    auto fields_itr = boost::make_split_iterator(line, kFieldFinder);
    static const auto kSplitEnd = boost::split_iterator<Iterator>();

    std::array<boost::iterator_range<Iterator>, 8> fields;
    for (size_t i = 0; i < 8; i++, ++fields_itr) {
      if (fields_itr == kSplitEnd)
        throw file_parse_error() << error_message("fewer than 8 required fields");
      fields[i] = *fields_itr;
    }

    model::impl::VariantContextData context;

    // Core variant fields (CHROM, POS, REF, ALT)
    context.contig_ = fields[0];
    parser::ParseRange(fields[1], parser::pos, context.pos_);
    context.ref_ = fields[3];
    SplitOptionalRange(context.alts_, fields[4], kCommaSplitter);

    // Context fields (ID, QUAL, FILTER)
    SplitOptionalRange(context.ids_, fields[2], kSemicolonSplitter);
    parser::ParseRange(fields[5], parser::qual, context.qual_);
    SplitOptionalRange(context.filters_, fields[6], kSemicolonSplitter);

    // INFO field
    for (auto i = boost::make_split_iterator(fields[7], kInfoFinder); i != kSplitEnd; ++i) {
      auto equals = boost::find(*i, kEqualsFinder);
      Attributes::key_type key(i->begin(), equals.begin());

      // Look up attribute key and parse value with returned parser, detect '=' if present
      // in case the attibute is not yet defined
      auto& parser = FindOrAddParser(info_keys_, key, header_.INFO(), equals.begin() == i->end());
      parser(equals.end(), i->end(), context.attrs_);
    }

    // Fix up variant based on INFO fields
    context.end_ = context.attrs_.at_or<util::Attributes::Integer>(
        VCFHeader::INFO::END, context.pos_ + context.ref_.size() - 1);

    return VariantContext(std::move(context));
  }

 private:
  VCFHeader& header_;

  typedef std::unordered_map<Attributes::key_type, std::unique_ptr<AttributeParser> > AttrTypes;
  AttrTypes info_keys_, format_keys_;

  static typename AttrTypes::mapped_type GetParser(const VCFHeader::Field& field) {
#define CASE(TYPE, SCALAR, VECTOR)                               \
  case VCFHeader::Field::Type::TYPE:                             \
    if (field.IsScalar()) {                                      \
      return std::make_unique<parser::SCALAR<Iterator> >(field); \
    } else {                                                     \
      return std::make_unique<parser::VECTOR<Iterator> >(field); \
    }

    switch (field.type_) {
      default:
        throw util::file_parse_error() << util::error_message("unsupported attribute type");
      case VCFHeader::Field::Type::FLAG:
        return std::make_unique<parser::FlagAttributeParser<Iterator> >(field);

        CASE(STRING, StringAttributeParser, StringsAttributeParser);
        CASE(INTEGER, IntegerAttributeParser, IntegersAttributeParser);
        CASE(FLOAT, FloatAttributeParser, FloatsAttributeParser);
        CASE(CHARACTER, CharacterAttributeParser, CharactersAttributeParser);
    }
#undef CASE
  }

  const AttributeParser& FindOrAddParser(AttrTypes& keys, util::Attributes::key_type key,
                                         VCFHeader::Fields& fields, bool as_flag) {
    typename AttrTypes::const_iterator i = keys.find(key);
    if (i == keys.end()) {
      VCFHeader::Field field(key, "Auto-generated");
      if (as_flag) {
        LOG(INFO) << fmt::format("Treating attribute {} as FLAG", key);
      } else {
        LOG(INFO) << fmt::format("Treating attribute {} as unbounded STRING vector", key);
        field.number_ = VCFHeader::Field::UNBOUNDED;
        field.type_ = VCFHeader::Field::Type::STRING;
      }
      VCFHeader::AddField(fields, field);
      i = keys.emplace(key, GetParser(field)).first;
    }
    return *(i->second);
  }
};

VariantContext ParseVCFVariant(const std::string& line, VCFHeader& header) {
  VCFVariantParser<const std::string> parser(header);
  return parser.ParseVCFVariant(line);
}

}  // namespace impl

void VCFSource::ParserDeleter::operator()(VCFSource::Parser* p) { delete p; }

VCFSource::VCFSource(FileFormat format, VCFSource::Reader&& reader)
    : header_(format), reader_(std::forward<Reader>(reader)) {
  auto line_parser =
      x3::with<impl::parser::header_tag>(std::ref(header_))[impl::parser::header_line];
  while (auto line = reader_->ReadNextLine()) {
    if (boost::starts_with(*line, "##")) {
      x3::parse(line->begin(), line->end(), line_parser);
      continue;
    } else if (boost::starts_with(*line, "#CHROM")) {
      // Parse required line with optional sample identifiers
      std::vector<Line> tokens;
      boost::split(tokens, *line, impl::kFieldSplitter);
      if (tokens.size() < 8 || tokens.size() == 9)
        throw util::file_parse_error() << util::error_message(
            fmt::format("unexpected number of fields ({}) on #CHROM line", tokens.size()));
      // if (tokens.size() > 9) header_.samples_.assign(tokens.begin() + 9,
      // tokens.end());
      break;
    } else {
      throw util::file_parse_error() << util::error_message("unexpected line in the header");
    }
  }

  // Add or update standard header fields
  header_.FILTER_AddField(VCFHeader::FILTER::PASS);

  // Create the variant parser for subsequent use
  parser_.reset(new impl::VCFVariantParser<Line>(header_));
}

FileFormat VCFSource::file_format() const { return header_.file_format(); }

VariantSourceInterface::NextResult VCFSource::NextVariant() {
  auto line = reader_->ReadNextLine();
  return line ? NextResult(std::move(parser_->ParseVCFVariant(*line))) : NextResult();
}
}  // namespace io
}  // namespace aseq