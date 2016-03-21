//
// Created by Michael Linderman on 12/16/15.
//

#include <ostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <glog/logging.h>
#include <cppformat/format.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>

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

  virtual void Write(const Line &line) { ostream_ << line; }

 private:
  std::unique_ptr<std::ostream> owned_ostream_;
  std::ostream &ostream_;
};

class BGZipLineWriter : public ASCIILineWriterInterface {
 public:
  BGZipLineWriter() = delete;

  BGZipLineWriter(const fs::path &path, FileFormat format)
      : file_(nullptr, &hts_close), format_(format) {
    file_.reset(hts_open(path.c_str(), "wz"));
    if (!file_) {
      throw file_write_error() << error_message("could not open tabix file for writing");
    }
  }
  ~BGZipLineWriter() {
    tbx_conf_t const *conf = nullptr;
    switch (format_) {
      default:
        break;
      case FileFormat::VCF4_1:
      case FileFormat::VCF4_2:
        conf = &tbx_conf_vcf;
    }
    if (conf) {
      fs::path path(file_->fn);
      file_.reset(nullptr);  // Force closure of file
      if (!fs::exists(path) || tbx_index_build(path.c_str(), 0, conf) != 0) {
        throw file_write_error() << error_message("could not generate tabix index");
      }
    }
  }

  virtual void Write(const Line &line) {
    bgzf_write(file_->fp.bgzf, line.begin(), line.end() - line.begin());
  }

 private:
  std::unique_ptr<htsFile, decltype(&hts_close)> file_;
  FileFormat format_;
};

}  // namespace impl

ASCIILineWriterInterface::FactoryResult ASCIILineWriterInterface::MakeLineWriter(
    std::ostream &ostream) {
  return FactoryResult(new impl::ASCIIStreamLineWriter(ostream));
}

ASCIILineWriterInterface::FactoryResult ASCIILineWriterInterface::MakeLineWriter(
    const fs::path &path, FileFormat format) {
  LOG_IF(INFO, fs::exists(path)) << fmt::format("Overwriting existing file {}", path.native());

  if (path.extension() == ".gz")
    return std::make_unique<impl::BGZipLineWriter>(path, format);
  else
    return std::make_unique<impl::ASCIIStreamLineWriter>(path);
}
}  // namespace io
}  // namespace aseq