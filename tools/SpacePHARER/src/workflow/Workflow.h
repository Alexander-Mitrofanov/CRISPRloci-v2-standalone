// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include <string>
#include <vector>
namespace spacepharer {
inline const char *enabled(bool value) { return value ? "TRUE" : nullptr; }
inline int execute_workflow(CommandCaller &caller, const std::string &program, const unsigned char *source,
                            std::size_t length, const std::vector<std::string> &arguments) {
    if (program.empty() || arguments.empty()) {
        Debug(Debug::ERROR) << "Workflow requires a path and input arguments\n";
        return EXIT_FAILURE;
    }
    FileUtil::writeFile(program, source, length);
    caller.execProgram(program.c_str(), arguments);
    return EXIT_FAILURE; // exec must not return successfully without execution.
}
} // namespace spacepharer
