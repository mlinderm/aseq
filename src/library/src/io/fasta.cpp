//
// Created by Michael Linderman on 2/28/16.
//

#include <ostream>
#include <fstream>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/utility/string_ref.hpp>

#include "aseq/io/fasta.hpp"

namespace aseq {
namespace io {

FastaSink::FastaSink(const boost::filesystem::path& path, size_t columns)
    : owned_ostream_(std::make_unique<std::ofstream>(path.native())),
      ostream_(*owned_ostream_),
      columns_(80) {}

FastaSink& FastaSink::PushSequence(const std::string& name, const std::string& sequence) {
  ostream_ << '>' << name << std::endl;
  auto sequence_ref = boost::string_ref(sequence);
  for (size_t i = 0; i < sequence.size(); i += columns_) {
    ostream_ << sequence_ref.substr(i, std::min(columns_, sequence.size() - i)) << std::endl;
  }
  return *this;
}
}
}