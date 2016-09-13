//
// Created by Michael Linderman on 12/13/15.
//

#pragma once

namespace aseq {
namespace io {

enum class FileFormat { UNKNOWN, VCF4_1, VCF4_2 };

inline bool IsVCF(FileFormat format) {
  return format == FileFormat::VCF4_1 || format == FileFormat::VCF4_2;
}
}
}