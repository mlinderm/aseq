//
// Created by Michael Linderman on 12/12/15.
//

#include <istream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>

#include "aseq/util/exception.hpp"
#include "aseq/io/line.hpp"

extern "C" {

long hts_utell(htsFile *fp);
}

using namespace aseq::util;
namespace fs = boost::filesystem;

namespace aseq {
namespace io {

namespace impl {

class ASCIIStreamLineReader : public ASCIILineReaderInterface {
 public:
  ASCIIStreamLineReader() = delete;

  explicit ASCIIStreamLineReader(std::istream &istream) : istream_(istream) {}

  explicit ASCIIStreamLineReader(const fs::path &path)
      : owned_istream_(std::make_unique<std::ifstream>(path.native())), istream_(*owned_istream_) {}

  virtual NextResult ReadNextLine() override {
    if (!std::getline(istream_, line_))
      return NextResult();
    else
      return NextResult(
          boost::make_iterator_range<Line::iterator>(line_.data(), line_.data() + line_.size()));
  }

 private:
  std::unique_ptr<std::istream> owned_istream_;
  std::istream &istream_;
  std::string line_;
};

namespace {
// Simple callback function to read bgzip compressed line
int TabixReadLine(BGZF *fp, void *tbxv, void *sv, int *tid, int *beg, int *end) {
  kstring_t *s = (kstring_t *)sv;
  *tid = 0;
  *beg = 0;
  *end = 0;
  return bgzf_getline(fp, '\n', s);
}
}

class TabixLineReader : public ASCIILineReaderInterface {
 public:
  TabixLineReader() = delete;

  explicit TabixLineReader(const fs::path &path)
      : file_(nullptr, &hts_close),
        index_(nullptr, &tbx_destroy),
        iter_(nullptr, &hts_itr_destroy),
        line_({}) {
    file_.reset(hts_open(path.c_str(), "r"));
    if (!file_) {
      throw file_parse_error() << error_message("could not open tabix file for reading");
    }
    index_.reset(tbx_index_load(path.c_str()));
    if (!index_) {
      throw file_parse_error() << error_message("Could not open .tbi index");
    }

    // Create synthetic iterator pointing to the start of the file
    iter_.reset((hts_itr_t *)calloc(1, sizeof(hts_itr_t)));
    iter_->read_rest = 1;
    iter_->curr_off = hts_utell(file_.get());
    iter_->readrec = TabixReadLine;
  }

  ~TabixLineReader() { free(ks_release(&line_)); }

  virtual NextResult ReadNextLine() override {
    // If we have a live iterator, use that to read from the file
    if (iter_ && tbx_itr_next(file_.get(), index_.get(), iter_.get(), &line_) < 0) {
      iter_.reset(nullptr);
      return NextResult();
    } else
      return NextResult(boost::make_iterator_range(line_.s, line_.s + line_.l));
  }

 private:
  std::unique_ptr<htsFile, decltype(&hts_close)> file_;
  std::unique_ptr<tbx_t, decltype(&tbx_destroy)> index_;
  std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iter_;
  kstring_t line_;
};

}  // impl namespace

ASCIILineReaderInterface::FactoryResult ASCIILineReaderInterface::MakeLineReader(
    std::istream &istream) {
  return std::make_unique<impl::ASCIIStreamLineReader>(istream);
}

ASCIILineReaderInterface::FactoryResult ASCIILineReaderInterface::MakeLineReader(
    const fs::path &path) {
  if (!fs::exists(path)) {  // Catch a common error case
    throw file_parse_error() << error_message("input file does not exist");
  }
  if (path.extension() == ".gz")
    return std::make_unique<impl::TabixLineReader>(path);
  else
    return std::make_unique<impl::ASCIIStreamLineReader>(path);
}
} // io namespace
} // aseq namespace