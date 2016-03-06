//
// Created by Michael Linderman on 12/12/15.
//

#pragma once

#include <stdexcept>
#include <typeinfo>

#include <boost/exception/all.hpp>

namespace aseq {
namespace util {

// Common Exception base_type
typedef boost::error_info<struct error_tag_msg, std::string> error_message;
typedef boost::error_info<struct error_tag_col_number, size_t> error_column_number;
typedef boost::error_info<struct error_tag_line_number, size_t> error_line_number;

struct exception_base : virtual std::exception, virtual boost::exception {
 public:
  const char* what() const throw() {
    if (std::string const* m = boost::get_error_info<error_message>(*this)) {
      return m->c_str();
    } else
      return typeid(*this).name();
  }
};

// Exception types
class invalid_argument : public exception_base {};
class file_parse_error : public exception_base {};
class file_write_error : public exception_base {};

class no_such_sample : public exception_base {};
}
}