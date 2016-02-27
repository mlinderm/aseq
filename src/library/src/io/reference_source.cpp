//
// Created by Michael Linderman on 2/22/16.
//

#include <memory>

#include <boost/filesystem.hpp>
#include <htslib/faidx.h>
#include <cppformat/format.h>
#include <string>
#include <bitset>

#include "aseq/util/exception.hpp"
#include "aseq/io/reference.hpp"

namespace aseq {
namespace io {

namespace fs = boost::filesystem;
using namespace aseq::util;

ReferenceSource::ReferenceSource(const fs::path& file)
    : faidx_(fai_load(file.c_str()), &fai_destroy) {
  if (!faidx_) {
    throw file_parse_error() << error_message(
        fmt::format("Could not open indexed FASTA file {}", file));
  }
}

ReferenceSource::ReferenceSource(const std::string& path) : ReferenceSource(fs::path(path)) {}
ReferenceSource::ReferenceSource() : faidx_(0, &fai_destroy) {}

std::string ReferenceSource::Sequence(const model::Contig& ctg, int64_t pos, int64_t end) {
  int length = 0;
  std::unique_ptr<char, void (*)(void*)> seq(
      faidx_fetch_seq(faidx_.get(), ctg.c_str(), pos - 1, end - 1, &length), &std::free);
  if (!seq) {
    // TODO: Improve error message: length=-1 no contig by that name, length=-2 seeking failed
    throw file_parse_error() << error_message("Error reading FASTA file");
  }

  // TODO: Get rid of this extra copy by adding function to htslib that takes a buffer
  return std::string(seq.get());
}
}
}