//
// Created by Michael Linderman on 3/6/16.
//

#include <iostream>
#include <algorithm>

#include <docopt.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include <cppformat/format.h>

#include "aseq-version.h"
#include "aseq/io/variant.hpp"
#include "aseq/io/reference.hpp"
#include "aseq/io/fasta.hpp"
#include "aseq/algorithm/variant.hpp"

#include "commands.hpp"

namespace fs = boost::filesystem;

namespace {
const char USAGE[] = R"(aseq Variant analysis commands

Usage:
  aseq variants (-h | --help)
  aseq variants split [-l <N>] <file>
  aseq variants normalize -R <ref> [--minimal] <file>
  aseq variants intervals [--flank <F>] <file>
  aseq variants consensus [--flank <F>] [--noREF | --noALT] -R <ref> <file>
  aseq variants table [--F <field>...] [--GF <field>...] <file>

Options:
  -h --help                  Show this screen.
  -R <ref>, --ref <ref>      Reference fasta (indexed)
  --flank <F>                Length of flanks [default: 1000]
  --noREF                    Don't emit reference consensus sequence
  --noALT                    Don't emit alternate consensus sequence
  --F <field>                INFO field
  --GF <field>               Genotype (FORMAT) field
  -l <N>                     Number of variants per split [default: 1]
  --minimal                  Sites-only output
)";

int SplitMain(std::map<std::string, docopt::value>& args) {
  using namespace aseq::io;

  long n = args["-l"].asLong();
  if (n < 0) {
    std::cerr << "split count must be > 0" << std::endl;
    std::cerr << USAGE;
    return 1;
  }

  std::string pattern = (fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.vcf").native();

  auto source = VariantSourceInterface::MakeVariantSource(args["<file>"].asString());
  for (auto v = source->NextVariant(); v;) {
    auto path = fs::unique_path(pattern);
    auto sink = VariantSinkInterface::MakeVariantSink(*source, path);
    for (size_t i = 0; i < n && v; i++) {
      sink->PushVariant(*v);
      v = source->NextVariant();
    }
    std::cout << path.native() << std::endl;
  }

  return 0;
}

int NormalizeMain(std::map<std::string, docopt::value>& args) {
  using namespace aseq::io;
  using namespace aseq::algorithm;

  ReferenceSource ref(args["--ref"].asString());
  auto source = VariantSourceInterface::MakeVariantSource(args["<file>"].asString());
  auto sink = VariantSinkInterface::MakeVariantSink(*source, std::cout, args["--minimal"].asBool());
  while (auto v = source->NextVariant()) {
    sink->PushVariant(Normalize(ref, std::move(*v)));
  }
  return 0;
}

int IntervalsMain(std::map<std::string, docopt::value>& args) {
  using namespace aseq::io;
  using aseq::model::Pos;

  Pos flank = args["--flank"].asLong();
  if (flank < 0) {
    std::cerr << "--flank argument must be > 0" << std::endl;
    std::cerr << USAGE;
    return 1;
  }

  auto source = VariantSourceInterface::MakeVariantSource(args["<file>"].asString());
  while (auto v = source->NextVariant()) {
    fmt::print(std::cout, "{}:{}-{}\n", v->contig(),
               std::max(v->pos() - flank, static_cast<Pos>(1)), v->end() + flank);
  }
  return 0;
}

int ConsensusMain(std::map<std::string, docopt::value>& args) {
  using namespace aseq::io;
  using namespace aseq::algorithm;

  aseq::model::Pos flank = args["--flank"].asLong();
  if (flank < 0) {
    std::cerr << "--flank argument must be > 0" << std::endl;
    std::cerr << USAGE;
    return 1;
  }

  ReferenceSource ref(args["--ref"].asString());
  auto source = VariantSourceInterface::MakeVariantSource(args["<file>"].asString());
  FastaSink sink(std::cout);

  auto v = source->NextVariant();
  if (!v) {
    LOG(ERROR) << "No variants found in input";
    return 1;
  }

  if (!args["--noREF"].asBool()) {
    sink.PushSequence("ref", ref.Sequence(v->contig(), v->pos() - flank, v->end() + flank));
  }
  if (!args["--noALT"].asBool()) {
    sink.PushSequence("alt", Consensus(ref, *v, flank));
  }

  if (source->NextVariant()) {
    LOG(WARNING) << "Multiple variants in input, but FASTA only emitted for the first";
  }
  return 0;
}

int VariantsToTableMain(std::map<std::string, docopt::value>& args) {
  using namespace aseq::io;
  using aseq::util::Attributes;

  auto missing = Attributes::mapped_type(std::string("NA"));
  auto& Fs = args["--F"].asStringList();
  auto& GFs = args["--GF"].asStringList();

  auto source = VariantSourceInterface::MakeVariantSource(args["<file>"].asString());
  auto& header = source->header();

  for (size_t i = 0; i < Fs.size(); i++) {
    fmt::print(std::cout, ((i == 0) ? "{}" : "\t{}"), Fs[i]);
  }
  for (size_t i = 0; i < header.NumSamples(); i++) {
    for (size_t j = 0; j < GFs.size(); j++) {
      fmt::print(std::cout, ((j == 0 && Fs.empty()) ? "{}.{}" : "\t{}.{}"), header.Sample(i),
                 GFs[j]);
    }
  }
  std::cout << std::endl;

  while (auto v = source->NextVariant()) {
    for (size_t i = 0; i < Fs.size(); i++) {
      fmt::print(std::cout, ((i == 0) ? "{}" : "\t{}"),
                 v->GetAttributeOr<Attributes::mapped_type>(Fs[i], missing));
    }
    for (size_t i = 0; i < header.NumSamples(); i++) {
      auto& gt = v->GetGenotype(header.Sample(i));
      for (size_t j = 0; j < GFs.size(); j++) {
        fmt::print(std::cout, ((j == 0 && Fs.empty()) ? "{}" : "\t{}"),
                   gt.GetAttributeOr<Attributes::mapped_type>(GFs[j], missing));
      }
    }
    std::cout << std::endl;
  }
  return 0;
}

}  // anonymous namespace

int VariantsMain(const std::vector<std::string>& argv) {
  std::map<std::string, docopt::value> args = docopt::docopt(USAGE, argv, true, ASEQ_VERSION);

  try {
    if (args["split"].asBool()) {
      return SplitMain(args);
    } else if (args["normalize"].asBool()) {
      return NormalizeMain(args);
    } else if (args["consensus"].asBool()) {
      return ConsensusMain(args);
    } else if (args["intervals"].asBool()) {
      return IntervalsMain(args);
    } else if (args["table"].asBool()) {
      return VariantsToTableMain(args);
    }
  } catch (aseq::util::exception_base& e) {
    LOG(ERROR) << e.what();
  }

  return 1;
}