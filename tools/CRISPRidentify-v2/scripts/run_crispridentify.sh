#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
caller_dir="$(pwd -P)"
env_dir="${CRISPRIDENTIFY_ENV:-${project_dir}/.envs/crispridentify-core}"
case "${env_dir}" in
    /*) ;;
    *) env_dir="${PWD}/${env_dir}" ;;
esac

caller_absolute_path() {
    local value="$1"
    case "${value}" in
        ""|/*) printf '%s' "${value}" ;;
        *) printf '%s/%s' "${caller_dir}" "${value}" ;;
    esac
}

# The Python entry point relies on the repository as its working directory for
# bundled models and databases. Resolve every user-facing filesystem argument
# first so paths continue to mean "relative to where the launcher was called".
rewritten_arguments=()
result_folder_seen=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --file|--input_folder|--input_folder_multifasta|\
        --result_folder|--pickle_report|--json_report|\
        --sequence_topology_manifest|--_folder_queue)
            option="$1"
            [[ "${option}" == "--result_folder" ]] && result_folder_seen=1
            rewritten_arguments+=("${option}")
            shift
            if [[ $# -gt 0 ]]; then
                rewritten_arguments+=("$(caller_absolute_path "$1")")
                shift
            fi
            ;;
        --file=*|--input_folder=*|--input_folder_multifasta=*|\
        --result_folder=*|--pickle_report=*|--json_report=*|\
        --sequence_topology_manifest=*|--_folder_queue=*)
            option="${1%%=*}"
            value="${1#*=}"
            [[ "${option}" == "--result_folder" ]] && result_folder_seen=1
            rewritten_arguments+=(
                "${option}=$(caller_absolute_path "${value}")"
            )
            shift
            ;;
        *)
            rewritten_arguments+=("$1")
            shift
            ;;
    esac
done

if [[ "${result_folder_seen}" -eq 0 ]]; then
    rewritten_arguments+=("--result_folder" "${caller_dir}/Results")
fi

if [[ ! -x "${env_dir}/bin/python" ]]; then
    echo "Error: CRISPRidentify environment not found at ${env_dir}." >&2
    echo "Run ./install.sh first." >&2
    exit 1
fi

export CONDA_PREFIX="${env_dir}"
export PATH="${env_dir}/bin:${PATH}"
export PYTHONHASHSEED=0
export PYTHONNOUSERSITE=1
export LC_ALL=C

# Reproduce Conda activation hooks, including Vmatch's symbol-map path.
shopt -s nullglob
for hook in "${env_dir}"/etc/conda/activate.d/*.sh; do
    # shellcheck source=/dev/null
    source "${hook}"
done
shopt -u nullglob

cd "${project_dir}"
exec "${env_dir}/bin/python" \
    "${project_dir}/CRISPRidentify.py" \
    "${rewritten_arguments[@]}"
