#!/usr/bin/env bash
set -euo pipefail
umask 0077

if [[ $# -ne 2 || $1 != /* || $2 != https://github.com/*/*.git ]]; then
    echo "Usage: $0 /absolute/frontend-checkout https://github.com/OWNER/FRONTEND.git" >&2
    exit 2
fi

script_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
source_path="${script_root}/frontend"
target_root="$1"
while [[ ${target_root} != / && ${target_root} == */ ]]; do
    target_root="${target_root%/}"
done
target_root="$(realpath -ms -- "${target_root}")"

if [[ ! -d ${source_path} || -L ${source_path} ]]; then
    echo "Public frontend export source must be a normal directory." >&2
    exit 1
fi
source_root="$(cd -- "${source_path}" && pwd -P)"
python3 "${script_root}/scripts/check_private_data.py" export-source --export "${source_root}"

if [[ -L ${target_root} || ! -d ${target_root}/.git || -L ${target_root}/.git ]]; then
    echo "Target must be a normal Git checkout directory." >&2
    exit 2
fi
resolved_target_root="$(cd -- "${target_root}" && pwd -P)"
if [[ ${resolved_target_root} != "${target_root}" ]]; then
    echo "Target checkout path cannot contain symbolic links." >&2
    exit 2
fi

remote_url="$(git -C "${target_root}" config --get remote.origin.url)"
if [[ "${remote_url}" != "$2" ]]; then
    echo "Target origin does not match the explicit frontend repository." >&2
    exit 2
fi
if [[ "${target_root}" == "${script_root}" || "${target_root}" == "${script_root}/"* || "${script_root}" == "${target_root}/"* ]]; then
    echo "Export target must be separate from the private project." >&2
    exit 2
fi

if [[ -n $(git -C "${target_root}" status --porcelain) ]]; then
    echo "Target checkout must be clean before export." >&2
    exit 1
fi

source_symlinks="$(find "${source_root}" \( -path "${source_root}/node_modules" -o -path "${source_root}/dist" -o -path "${source_root}/coverage" \) -prune -o -type l -print)"
if [[ -n ${source_symlinks} ]]; then
    echo "Symbolic links are forbidden in the public frontend export source:" >&2
    printf '%s\n' "${source_symlinks}" >&2
    exit 1
fi

rsync -a --no-links --delete \
    --include=.env.example \
    --exclude='.env*' \
    --exclude='*.local' \
    --exclude='*.log' \
    --exclude='*.pem' \
    --exclude='*.key' \
    --exclude=.DS_Store \
    --exclude=.git \
    --exclude=node_modules \
    --exclude=dist \
    --exclude=coverage \
    "${source_root}/" "${target_root}/"

unexpected_local_files="$(find "${target_root}" \( -path "${target_root}/.git" -o -path "${target_root}/node_modules" -o -path "${target_root}/dist" \) -prune -o -type f \( -name '.env*' -o -name '*.local' -o -name '*.log' -o -name '*.pem' -o -name '*.key' \) ! -name '.env.example' -print)"
if [[ -n ${unexpected_local_files} ]]; then
    echo "Private or local-only files escaped into public export:" >&2
    printf '%s\n' "${unexpected_local_files}" >&2
    exit 1
fi

unexpected_symlinks="$(find "${target_root}" \( -path "${target_root}/.git" -o -path "${target_root}/node_modules" -o -path "${target_root}/dist" \) -prune -o -type l -print)"
if [[ -n ${unexpected_symlinks} ]]; then
    echo "Symbolic links escaped into public export:" >&2
    printf '%s\n' "${unexpected_symlinks}" >&2
    exit 1
fi

for private_path in backend deploy deployment reports CRISPRidentify SpacerPlacer CRISPR-evOr Integration; do
    if [[ -e ${target_root}/${private_path} ]]; then
        echo "Private workspace path escaped into public export: ${private_path}" >&2
        exit 1
    fi
done

python3 "${script_root}/scripts/check_private_data.py" export --export "${target_root}"
git -C "${target_root}" diff --check
git -C "${target_root}" status --short
