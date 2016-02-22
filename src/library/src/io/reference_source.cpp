//
// Created by Michael Linderman on 2/22/16.
//

#include <memory>

#include <boost/filesystem.hpp>
//#include <htslib/faidx.h>

#include "aseq/util/exception.hpp"
#include "aseq/io/reference.hpp"

namespace aseq {
namespace io {

namespace fs = boost::filesystem;
using namespace aseq::util;

    /*
ReferenceSource::ReferenceSource(const fs::path& file)
    : faidx_(fai_load(file.c_str()), &fai_destroy) {
  if (!faidx_) throw file_parse_error() << error_message("Could not open indexed FASTA file");
}
     */
}
}