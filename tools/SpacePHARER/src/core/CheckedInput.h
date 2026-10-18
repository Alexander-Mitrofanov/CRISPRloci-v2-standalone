// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "Debug.h"
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
namespace spacepharer {
inline void require_columns(std::size_t actual, std::size_t minimum) {
    if (actual < minimum || actual > 255) {
        Debug(Debug::ERROR) << "Malformed native record: invalid column count\n";
        EXIT(EXIT_FAILURE);
    }
}
inline double finite_number(const char *value) {
    char *end = nullptr;
    errno = 0;
    const double number = std::strtod(value, &end);
    if (end == value || errno == ERANGE || !std::isfinite(number) ||
        (*end != '\0' && *end != '\t' && *end != '\n' && *end != ' ')) {
        Debug(Debug::ERROR) << "Malformed native record: invalid numeric field\n";
        EXIT(EXIT_FAILURE);
    }
    return number;
}
// Database readers require explicit close before destruction. This noncopyable
// guard also closes them on future early-return paths.
template <class Reader> class ReaderGuard {
    Reader &reader_;

  public:
    explicit ReaderGuard(Reader &reader) : reader_(reader) {}
    ReaderGuard(const ReaderGuard &) = delete;
    ReaderGuard &operator=(const ReaderGuard &) = delete;
    ~ReaderGuard() { reader_.close(); }
};
} // namespace spacepharer
