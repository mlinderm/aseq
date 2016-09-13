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
(aseq::io::VCFHeader::Field::Number, nmbr_)
(aseq::io::VCFHeader::Field::Type, type_)
(std::string, desc_)
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
static model::Genotype::Alleles::initializer genotype_alleles_init;

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
  return !IsDefined(input) ? container : boost::split(container, input, predicate);
}

// Boost Spirit X3 infrastructure
namespace parser {

// "as" directive needed to parse into underlying type of attribute, but then convert
// that typed attribute into the generic "any" used as the attribute mapped type
template <typename Subject, typename T>
struct as_directive : x3::unary_parser<Subject, as_directive<Subject, T> > {
  typedef x3::unary_parser<Subject, as_directive<Subject, T> > base_type;
  typedef T attribute_type;
  typedef Subject subject_type;
  static bool const handles_container = Subject::handles_container;

  as_directive(const Subject& subject) : base_type(subject) {}

  template <typename Iterator, typename Context, typename RContext, typename Attribute>
  bool parse(Iterator& first, Iterator const& last, Context const& context, RContext& rcontext,
             Attribute& attr) const {
    T as_attr;
    if (this->subject.parse(first, last, context, rcontext, as_attr)) {
      attr = std::move(as_attr);
      return true;
    } else {
      return false;
    }
  }
};

template <typename T>
struct as {
  template <typename Subject>
  as_directive<typename x3::extension::as_parser<Subject>::value_type, T> operator[](
      const Subject& subject) const {
    return {x3::as_parser(subject)};
  }
};

// Parser rules
using util::Attributes;
using model::VariantContext;
using model::Genotype;
using x3::lit;
using x3::char_;

x3::rule<class header_line> const header_line = "VCF header line";
x3::rule<class header_field> const header_field = "VCF header field";
x3::rule<class field_attributes, VCFHeader::Field> const field_attributes =
    "VCF header field attributes";
x3::rule<class filter_or_alt_attributes, VCFHeader::Field> const filter_or_alt_attributes =
    "VCF header field with ID and description attributes";

struct header_tag {};
auto SetFormat = [](auto& ctx) {
  VCFHeader& header = x3::get<header_tag>(ctx);
  header.set_file_format(x3::_attr(ctx));
};

#define FIELD_ADD(kind)                           \
  auto Add##kind##Field = [](auto& ctx) {        \
    VCFHeader& header = x3::get<header_tag>(ctx); \
    header.Add##kind##Field(x3::_attr(ctx));      \
  };

FIELD_ADD(FILTER)
FIELD_ADD(INFO)
FIELD_ADD(FORMAT)

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
    | (lit("FORMAT") > field_attributes)[AddFORMATField]
    | (lit("FILTER") > filter_or_alt_attributes)[AddFILTERField]
    | ((+x3::alnum) > '=' > (*(char_ - x3::eol)))
);

auto const header_line_def = lit("##") > (
    (lit("fileformat=") > kFileFormats[SetFormat])
    | header_field
) > -x3::eol;
// clang-format on

x3::rule<class ref_allele, model::Pos> const pos = "POS";
x3::rule<class qual, VariantContext::Qual> const qual = "QUAL";

#define ATTR_RULE(name, desc) x3::rule<class name, Attributes::mapped_type> const name = desc
ATTR_RULE(integer_value, "Integer value");
ATTR_RULE(integers_value, "Integers value");
ATTR_RULE(float_value, "Float value");
ATTR_RULE(floats_value, "Floats value");
ATTR_RULE(char_value, "Character value");
ATTR_RULE(chars_value, "Characters value");
#undef ATTR_RULE

x3::rule<class genotype_allele, int> const genotype_allele = "GT allele";
x3::rule<class genotype_alleles, model::impl::PhasedIndices> const genotype_alleles = "GT";
x3::rule<class genotype_value, Attributes::mapped_type> const genotype_value = "GT attribute";

// Can't use static values in Genotype due to initialization issues
const x3::symbols<Genotype::Alleles> kGenotypeStrings{
    {".", Genotype::Alleles(VariantContext::kNoCallIdx)},
    {"0", Genotype::Alleles(0)},
    {"1", Genotype::Alleles(1)},
    {"./.", Genotype::Alleles(true, VariantContext::kNoCallIdx, VariantContext::kNoCallIdx)},
    {"0/0", Genotype::Alleles(true, 0, 0)},
    {"0/1", Genotype::Alleles(false, 0, 1)},
    {"1/0", Genotype::Alleles(false, 0, 1)},
    {"1/1", Genotype::Alleles(true, 1, 1)},
    {"0|0", Genotype::Alleles(true, 0, 0)},
    {"0|1", Genotype::Alleles(true, 0, 1)},
    {"1|0", Genotype::Alleles(true, 1, 0)},
    {"1|1", Genotype::Alleles(true, 1, 1)}};

// clang-format off
auto const pos_def  = x3::long_long;
auto const qual_def = lit('.') | x3::float_;
auto const missing  = (lit('.') % ',');

auto const integer_value_def  = missing | as<Attributes::Integer>()[x3::int_];
auto const integers_value_def = missing | as<Attributes::Integers>()[x3::int_ % ','];
auto const float_value_def    = missing | as<Attributes::Float>()[x3::float_];
auto const floats_value_def   = missing | as<Attributes::Floats>()[x3::float_ % ','];
auto const char_value_def     = missing | as<Attributes::Character>()[x3::char_];
auto const chars_value_def    = missing | as<Attributes::Characters>()[x3::char_ % ','];

auto missing_allele = [](auto& ctx) { x3::_val(ctx) = VariantContext::kNoCallIdx; };
auto init_phase     = [](auto& ctx) { x3::_val(ctx).phased_ = true; };
auto update_phase   = [](auto& ctx) { x3::_val(ctx).phased_ &= (x3::_attr(ctx) == '|'); };
auto add_allele     = [](auto& ctx) { x3::_val(ctx).indices_.push_back(x3::_attr(ctx)); };

// Use symbol table to catch common genotypes
auto const genotype_value_def  =
    (kGenotypeStrings >> x3::eoi) | as<Genotype::Alleles>()[genotype_alleles];
auto const genotype_alleles_def =
    x3::eps[init_phase] >> (genotype_allele[add_allele] % x3::char_("/|")[update_phase]);
auto const genotype_allele_def = x3::rule<class genotype_allele, int> {}
    %= x3::lit('.')[missing_allele] | x3::int_;
// clang-format on

BOOST_SPIRIT_DEFINE(header_line, header_field, field_attributes, filter_or_alt_attributes);
BOOST_SPIRIT_DEFINE(pos, qual);
BOOST_SPIRIT_DEFINE(integer_value, integers_value, float_value, floats_value, char_value,
                    chars_value);
BOOST_SPIRIT_DEFINE(genotype_allele, genotype_alleles, genotype_value);

template <typename Iterator>
struct AttributeParser {
  Attributes::key_type key_;

  AttributeParser() = delete;
  AttributeParser(const Attributes::key_type& key) : key_(key) {}

  virtual bool Parse(Iterator, const Iterator&, Attributes::mapped_type&) const { return true; };

  void Parse(Iterator begin, const Iterator& end, Attributes& attr) const {
    Attributes::mapped_type value;
    if (!this->Parse(begin, end, value)) {
      throw util::file_parse_error()
          << util::error_message(fmt::format("failed to parse value for attribute {}", key_));
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
  bool Parse(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const {
    return (begin == end) ? value = true, true : false;
  }
};

template <typename Iterator>
struct StringAttributeParser : public AttributeParser<Iterator> {
  StringAttributeParser(const Attributes::key_type key) : AttributeParser<Iterator>(key) {}
  bool Parse(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const {
    if (IsDefined(begin, end)) {
      value = std::move(Attributes::String(begin, end));
    }
    return true;
  }
};

template <typename Iterator>
struct StringsAttributeParser : public AttributeParser<Iterator> {
  StringsAttributeParser(const Attributes::key_type key) : AttributeParser<Iterator>(key) {}
  bool Parse(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const {
    Attributes::Strings strings;
    SplitOptionalRange(strings, boost::make_iterator_range(begin, end), kCommaSplitter);
    if (!strings.empty()) {
      value = std::move(strings);
    }
    return true;
  }
};

#define PARSER(NAME, EXPR)                                                                  \
  template <typename Iterator>                                                              \
  struct NAME : public AttributeParser<Iterator> {                                          \
    NAME(const Attributes::key_type key) : AttributeParser<Iterator>(key) {}                \
    bool Parse(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const { \
      return x3::parse(begin, end, EXPR, value) && begin == end;                            \
    }                                                                                       \
  };

PARSER(IntegerAttributeParser, parser::integer_value);
PARSER(IntegersAttributeParser, parser::integers_value);
PARSER(FloatAttributeParser, parser::float_value);
PARSER(FloatsAttributeParser, parser::floats_value);
PARSER(CharacterAttributeParser, parser::char_value);
PARSER(CharactersAttributeParser, parser::chars_value);
PARSER(GenotypeAttributeParser, parser::genotype_value);
#undef PARSER

template <typename Iterator>
struct FilterAttributeParser : public AttributeParser<Iterator> {
  FilterAttributeParser(const Attributes::key_type key) : AttributeParser<Iterator>(key) {}
  bool Parse(Iterator begin, const Iterator& end, Attributes::mapped_type& value) const {
    VariantContext::Filters filters;
    SplitOptionalRange(filters, boost::make_iterator_range(begin, end), kSemicolonSplitter);
    if (!filters.empty()) {
      value = std::move(filters);
    }
    return true;
  }
};

}  // namespace parser

template <typename Range, typename Parser, typename Attribute>
void ParseRange(const Range& range, const Parser& parser, Attribute& attr) {
  typename Range::iterator begin = range.begin();
  if (!x3::parse(begin, range.end(), parser, attr) || begin != range.end()) {
    throw util::file_parse_error()
        << util::error_message(fmt::format("failed to parse {}", parser.name));
  }
};

template <typename Line>
class VCFVariantParser {
  typedef typename Line::const_iterator Iterator;
  typedef parser::AttributeParser<Iterator> AttributeParser;
  using Attributes = util::Attributes;

 public:
  VCFVariantParser(VCFHeader& header) : header_(header) {
    // "Missing" INFO field
    info_keys_.emplace(".", std::make_unique<AttributeParser>("."));
    for (auto f : header.INFOValues()) {
      info_keys_.emplace(f, GetParser(f));
    }
    for (auto f : header.FORMATValues()) {
      format_keys_.emplace(f, GetParser(f));
    }
  }

  VariantContext ParseVCFVariant(const Line& line) {
    using util::file_parse_error;
    using util::error_message;
    using model::Genotype;

    // Split VCF line into component fields
    auto fields_itr = boost::make_split_iterator(line, kFieldFinder);
    static const auto kSplitEnd = boost::split_iterator<Iterator>();

    std::array<boost::iterator_range<Iterator>, 8> fields;
    for (size_t i = 0; i < 8; i++, ++fields_itr) {
      if (fields_itr == kSplitEnd)
        throw file_parse_error() << error_message(fmt::format("{} fields found, 8 required", i));
      fields[i] = *fields_itr;
    }

    model::impl::VariantContextData data;

    // Core variant fields (CHROM, POS, REF, ALT)
    data.contig_ = fields[0];
    ParseRange(fields[1], parser::pos, data.pos_);
    data.ref_ = fields[3];
    SplitOptionalRange(data.alts_, fields[4], kCommaSplitter);

    // Context fields (ID, QUAL, FILTER)
    SplitOptionalRange(data.ids_, fields[2], kSemicolonSplitter);
    ParseRange(fields[5], parser::qual, data.qual_);
    SplitOptionalRange(data.filters_, fields[6], kSemicolonSplitter);

    // INFO field
    for (auto i = boost::make_split_iterator(fields[7], kInfoFinder); i != kSplitEnd; ++i) {
      auto equals = boost::find(*i, kEqualsFinder);
      Attributes::key_type key(i->begin(), equals.begin());

      // Look up attribute key and parse value with returned parser, detect if flag (or valued)
      // in case the attribute is not yet defined
      auto& parser = FindOrAddParser(info_keys_, key, header_.INFO(), equals.begin() == i->end());
      parser.Parse(equals.end(), i->end(), data.attrs_);
    }

    // Fix up variant based on INFO fields
    data.end_ = data.attrs_.at_or<util::Attributes::Integer>(VCFHeader::INFO::END,
                                                             (data.pos_ + data.ref_.size()) - 1);

    VariantContext context(std::move(data));

    // FORMAT and genotypes (if relevant)
    if (header_.NumSamples() > 0) {
      if (fields_itr == kSplitEnd) {
        throw file_parse_error() << error_message("missing FORMAT field");
      }

      // Parse FORMAT field
      std::vector<const AttributeParser*> format;
      for (auto f = boost::make_split_iterator(*fields_itr, kFormatFinder); f != kSplitEnd; ++f) {
        format.push_back(&FindOrAddParser(format_keys_, *f, header_.FORMAT()));
      }
      ++fields_itr;

      // Parse individual sample entries
      for (size_t g = 0; g < header_.NumSamples(); g++, ++fields_itr) {
        if (fields_itr == kSplitEnd) {
          throw file_parse_error() << error_message(
              fmt::format("expected {} sample entries, found {}", header_.NumSamples(), g));
        }

        Attributes sample_attr;

        auto f = boost::make_split_iterator(*fields_itr, kFormatFinder);
        for (auto p = format.begin(); p != format.end() && f != kSplitEnd; ++p, ++f) {
          auto& parser = *(*p);
          parser.Parse(f->begin(), f->end(), sample_attr);
        }
        if (f != kSplitEnd) {
          throw file_parse_error() << error_message("more attributes than specified in FORMAT");
        }

        Genotype::Alleles alleles;
        {  // Get and remove GT from attributes if present
          auto gt = sample_attr.find(VCFHeader::FORMAT::GT);
          if (gt != sample_attr.end()) {
            alleles = Attributes::at<Genotype::Alleles>(gt);
            sample_attr.erase(gt);
          }
        }

        context.AddGenotype(header_.sample(g), alleles, std::move(sample_attr));
      }
      if (fields_itr != kSplitEnd) {
        throw file_parse_error() << error_message("more samples than specified in header");
      }
    }

    return context;
  }

 private:
  VCFHeader& header_;

  typedef std::unordered_map<Attributes::key_type, std::unique_ptr<AttributeParser> > AttrTypes;
  AttrTypes info_keys_, format_keys_;

  static typename AttrTypes::mapped_type GetParser(const VCFHeader::Field& field) {
#define CASE1(TYPE, PARSER)          \
  case VCFHeader::Field::Type::TYPE: \
    return std::make_unique<parser::PARSER<Iterator> >(field);

#define CASE2(TYPE, SCALAR, VECTOR)                              \
  case VCFHeader::Field::Type::TYPE:                             \
    if (field.IsScalar()) {                                      \
      return std::make_unique<parser::SCALAR<Iterator> >(field); \
    } else {                                                     \
      return std::make_unique<parser::VECTOR<Iterator> >(field); \
    }

    switch (field.type_) {
      CASE1(FLAG, FlagAttributeParser);
      CASE1(GENOTYPE, GenotypeAttributeParser);
      CASE1(FILTER, FilterAttributeParser);
      CASE2(STRING, StringAttributeParser, StringsAttributeParser);
      CASE2(INTEGER, IntegerAttributeParser, IntegersAttributeParser);
      CASE2(FLOAT, FloatAttributeParser, FloatsAttributeParser);
      CASE2(CHARACTER, CharacterAttributeParser, CharactersAttributeParser);
      default:
        throw util::file_parse_error()
            << util::error_message(fmt::format("unsupported attribute type for {}", field));
    }
#undef CASE2
#undef CASE1
  }

  const AttributeParser& FindOrAddParser(AttrTypes& keys, util::Attributes::key_type key,
                                         VCFHeader::Fields& fields, bool as_flag = false) {
    typename AttrTypes::const_iterator i = keys.find(key);
    if (i == keys.end()) {
      VCFHeader::Field field(key, "Auto-generated");
      if (as_flag) {
        LOG(INFO) << fmt::format("Treating attribute {} as FLAG", key);
      } else {
        LOG(INFO) << fmt::format("Treating attribute {} as unbounded STRING vector", key);
        field.nmbr_ = VCFHeader::Field::UNBOUNDED;
        field.type_ = VCFHeader::Field::Type::STRING;
      }
      VCFHeader::AddField(fields, field);
      std::tie(i, std::ignore) = keys.emplace(key, GetParser(field));
    }
    return *(i->second);
  }
};

VariantContext ParseVCFVariantLine(const std::string& line, VCFHeader& header) {
  VCFVariantParser<const std::string> parser(header);
  return parser.ParseVCFVariant(line);
}

}  // namespace impl

void VCFSource::ParserDeleter::operator()(VCFSource::Parser* p) { delete p; }

VCFSource::VCFSource(FileFormat format, VCFSource::Reader&& reader)
    : header_(format), reader_(std::forward<Reader>(reader)) {
  using util::file_parse_error;

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
      if (tokens.size() < 8 || tokens.size() == 9) {
        throw file_parse_error() << util::error_message(
            fmt::format("unexpected number of fields ({}) on #CHROM line", tokens.size()));
      } else if (tokens.size() > 9) {
        // File has samples
        header_.samples_.assign(tokens.begin() + 9, tokens.end());
      }
      break;
    } else {
      throw file_parse_error() << util::error_message("unexpected line in the header");
    }
  }

  // Add or update standard header fields
  header_.AddFILTERField(VCFHeader::FILTER::PASS);
  header_.FORMAT_[VCFHeader::FORMAT::GT] = VCFHeader::FORMAT::GT;
  header_.FORMAT_[VCFHeader::FORMAT::FT] = VCFHeader::FORMAT::FT;

  // Create the variant parser for subsequent use
  parser_.reset(new impl::VCFVariantParser<Line>(header_));
}

FileFormat VCFSource::file_format() const { return header_.file_format(); }



void VCFSource::SetRegion(model::Contig contig, model::Pos pos, model::Pos end) {
  reader_->SetRegion(contig, pos, end);
}

VariantSourceInterface::NextResult VCFSource::NextVariant() {
  if (auto line = reader_->ReadNextLine())
    return NextResult(parser_->ParseVCFVariant(*line));
  else
    return NextResult();
}
}  // namespace io
}  // namespace aseq