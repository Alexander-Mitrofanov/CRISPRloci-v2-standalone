#include "core/RecordIndex.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

int main() {
    const auto key = [](const char *row) { return static_cast<unsigned int>(std::strtoul(row, nullptr, 10)); };
    const auto next = [](char *row) { return std::strchr(row, '\n') + 1; };
    std::string data;
    for (unsigned int i = 0; i < 10000; ++i) {
        data += std::to_string((i * 13) % 101) + "\t" + std::to_string(i) + "\n";
    }
    std::size_t indexedRows = 0;
    const auto index = spacepharer::index_records(data.data(), key, [&](char *row) {
        ++indexedRows;
        return next(row);
    });
    assert(indexedRows == 10000);
    for (unsigned int target = 0; target < 101; ++target) {
        // Independently retain the original full-scan algorithm as the oracle.
        std::vector<char *> expected;
        for (char *row = data.data(); *row; row = next(row)) {
            if (key(row) == target) expected.push_back(row);
        }
        assert(index.at(target) == expected);
        assert(index.at(target).back() == expected.back());
    }
    assert(index.find(101) == index.end());
    char empty[] = "";
    assert(spacepharer::index_records(empty, key, next).empty());
}
