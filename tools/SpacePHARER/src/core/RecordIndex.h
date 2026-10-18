// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <unordered_map>
#include <vector>

namespace spacepharer {
// Borrow database rows for one query. Each target retains its original row order;
// callers own parsing, reader lifetime and scientific last-match/tie policies.
template <class KeyOf, class NextLine>
inline std::unordered_map<unsigned int, std::vector<char *>>
index_records(char *data, KeyOf key_of, NextLine next_line) {
    std::unordered_map<unsigned int, std::vector<char *>> index;
    for (char *row = data; *row != '\0'; row = next_line(row)) {
        index[key_of(row)].push_back(row);
    }
    return index;
}
} // namespace spacepharer
