#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
profile="core"
use_lock=1

usage() {
    cat <<'EOF'
Usage: scripts/setup_local.sh [--profile core|is-elements] [--no-lock]

The core profile supports array detection, classification, portable strand
prediction, CasAndra Cas annotation, and every standard report.
The is-elements profile additionally installs fasta36 and the IS-element model
for --is_element True.

Existing environments are reconciled with the selected profile before the
vendored, checksummed CasAndra package is installed without network access.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --profile)
            if [[ $# -lt 2 ]]; then
                echo "Error: --profile requires a value." >&2
                exit 2
            fi
            profile="$2"
            shift 2
            ;;
        --no-lock)
            use_lock=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown setup option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

host_os="$(uname -s)"
host_arch="$(uname -m)"
if [[ "${host_os}" != "Linux" || "${host_arch}" != "x86_64" ]]; then
    echo "Error: CRISPRidentify v2 supports Linux x86_64 only; detected ${host_os} ${host_arch}." >&2
    exit 1
fi

case "${profile}" in
    core)
        environment_file="${project_dir}/environment.yml"
        default_env_dir="${project_dir}/.envs/crispridentify-core"
        ;;
    is-elements)
        environment_file="${project_dir}/environments/is-elements.yml"
        default_env_dir="${project_dir}/.envs/crispridentify-is-elements"
        ;;
    *)
        echo "Error: profile must be 'core' or 'is-elements'." >&2
        exit 2
        ;;
esac

env_dir="${CRISPRIDENTIFY_ENV:-${default_env_dir}}"
case "${env_dir}" in
    /*) ;;
    *) env_dir="${PWD}/${env_dir}" ;;
esac

lock_file="${project_dir}/locks/conda-linux-64-core.lock"
casandra_dir="${project_dir}/vendor/casandra"

if [[ ! -f "${casandra_dir}/pyproject.toml" ]]; then
    echo "Error: vendored CasAndra source is missing at ${casandra_dir}." >&2
    exit 1
fi
casandra_source_manifest="${casandra_dir}/SOURCE_MANIFEST.sha256"
if [[ ! -f "${casandra_source_manifest}" ]]; then
    echo "Error: vendored CasAndra source manifest is missing." >&2
    exit 1
fi
if ! command -v sha256sum >/dev/null 2>&1; then
    echo "Error: sha256sum is required to verify vendored CasAndra." >&2
    exit 1
fi
echo "Verifying vendored CasAndra source manifest"
(
    cd "${casandra_dir}"
    sha256sum --check --strict --quiet "${casandra_source_manifest}"
)


if command -v mamba >/dev/null 2>&1; then
    solver="$(command -v mamba)"
elif command -v micromamba >/dev/null 2>&1; then
    solver="$(command -v micromamba)"
elif command -v conda >/dev/null 2>&1; then
    solver="$(command -v conda)"
else
    echo "Error: install Conda, Mamba, or Micromamba before running this script." >&2
    exit 1
fi

locked_core=0
if [[ \
        "${profile}" == "core" \
        && "${use_lock}" -eq 1 \
        && -f "${lock_file}" ]]; then
    locked_core=1
fi

if [[ "${locked_core}" -eq 1 ]]; then
    if [[ ! -x "${env_dir}/bin/python" ]]; then
        echo "Creating locked core environment: ${env_dir}"
        "${solver}" create \
            --prefix "${env_dir}" \
            --file "${lock_file}" \
            --yes
    else
        echo "Reconciling locked core environment: ${env_dir}"
        "${solver}" install \
            --prefix "${env_dir}" \
            --file "${lock_file}" \
            --yes
    fi
elif [[ ! -x "${env_dir}/bin/python" ]]; then
    echo "Solving ${profile} environment: ${env_dir}"
    "${solver}" env create \
        --prefix "${env_dir}" \
        --file "${environment_file}" \
        --yes
else
    echo "Reconciling ${profile} environment: ${env_dir}"
    "${solver}" env update \
        --prefix "${env_dir}" \
        --file "${environment_file}" \
        --prune \
        --yes
fi

export CONDA_PREFIX="${env_dir}"
export PATH="${env_dir}/bin:${PATH}"
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1

# Activation hooks configure the pinned Conda compiler and Vmatch resources.
shopt -s nullglob
for hook in "${env_dir}"/etc/conda/activate.d/*.sh; do
    # shellcheck source=/dev/null
    source "${hook}"
done
shopt -u nullglob

casandra_staging="$(mktemp -d "${TMPDIR:-/tmp}/crispridentify-casandra.XXXXXX")"
cleanup_casandra_staging() {
    rm -rf -- "${casandra_staging}"
}
trap cleanup_casandra_staging EXIT
cp -a "${casandra_dir}/." "${casandra_staging}/"

echo "Installing vendored CasAndra into ${env_dir}"
"${env_dir}/bin/python" -m pip install \
    --disable-pip-version-check \
    --no-index \
    --no-deps \
    --no-build-isolation \
    --force-reinstall \
    "${casandra_staging}"
"${env_dir}/bin/python" -m pip check

"${project_dir}/scripts/build_native_matcher.sh"

"${env_dir}/bin/python" \
    "${project_dir}/scripts/check_environment.py" \
    --profile "${profile}"

echo
echo "Setup complete."
echo "Run the smoke test with:"
if [[ "${env_dir}" != "${default_env_dir}" ]]; then
    echo "  CRISPRIDENTIFY_ENV=${env_dir} \\"
fi
echo "  ${project_dir}/scripts/run_crispridentify.sh \\"
echo "    --file examples/NC_017040.1.fasta \\"
echo "    --result_folder runs/smoke \\"
echo "    --parallel False --cas True --html_report True"
