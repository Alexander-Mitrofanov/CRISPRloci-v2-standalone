#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
native_dir="${project_dir}/components/native"
source_file="${native_dir}/literal_fuzzy.c"
path_feature_source="${native_dir}/path_features.c"
output_file="${native_dir}/libliteral_fuzzy.so"
compiler="${CC:-x86_64-conda-linux-gnu-cc}"

if ! command -v "${compiler}" >/dev/null 2>&1; then
    echo "Error: the CRISPRidentify Conda C compiler is unavailable." >&2
    echo "Run scripts/setup_local.sh to install and activate the pinned toolchain." >&2
    exit 1
fi

compiler_path="$(command -v "${compiler}")"
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    case "${compiler_path}" in
        "${CONDA_PREFIX}"/bin/*) ;;
        *)
            echo "Error: compiler resolved outside the CRISPRidentify environment: ${compiler_path}" >&2
            exit 1
            ;;
    esac
fi

temporary_file="$(mktemp --tmpdir="${native_dir}" .libliteral_fuzzy.XXXXXX.so)"
cleanup() {
    rm -f "${temporary_file}"
}
trap cleanup EXIT

"${compiler_path}" \
    -O3 -std=c99 -fPIC -shared -Wall -Wextra -Werror \
    -fno-fast-math -ffp-contract=off \
    "${source_file}" "${path_feature_source}" -lm -o "${temporary_file}"
chmod 0755 "${temporary_file}"
mv "${temporary_file}" "${output_file}"
trap - EXIT

echo "Built native fuzzy matcher: ${output_file}"
