//
// Created by Michael Linderman on 12/16/15.
//

#include <ostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <glog/logging.h>

#include "aseq/util/exception.hpp"
#include "aseq/io/line.hpp"

using namespace aseq::util;
namespace fs = boost::filesystem;

namespace aseq {
namespace io {

namespace impl {

class ASCIIStreamLineWriter : public ASCIILineWriterInterface {
 public:
  ASCIIStreamLineWriter() = delete;

  explicit ASCIIStreamLineWriter(std::ostream &ostream) : ostream_(ostream) {}

  explicit ASCIIStreamLineWriter(const fs::path &path)
      : owned_ostream_(std::make_unique<std::ofstream>(path.native())), ostream_(*owned_ostream_) {}

  virtual void WriteLine(const Line &line) { ostream_ << line << std::endl; }

 private:
  std::unique_ptr<std::ostream> owned_ostream_;
  std::ostream &ostream_;
};
}  // namespace impl

ASCIILineWriterInterface::FactoryResult ASCIILineWriterInterface::MakeLineWriter(
    std::ostream &ostream) {
  return FactoryResult(new impl::ASCIIStreamLineWriter(ostream));
}

ASCIILineWriterInterface::FactoryResult ASCIILineWriterInterface::MakeLineWriter(
    const fs::path &path) {
  LOG_IF(INFO, fs::exists(path)) << "Overwriting existing file";

  return std::make_unique<impl::ASCIIStreamLineWriter>(path);
}
}  // namespace io
}  // namespace aseq