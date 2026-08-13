#!/usr/bin/env bash
set -euo pipefail

source ../../config.env

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

echo "Checking paths..."

for path in \
    "$DATADIR" \
    "$UV_PATH" \
    "$CELLRANGER_PATH" \
    "$STAR_PATH" \
    "$CELLBENDER_PATH"
do
    if [[ ! -e "$path" ]]; then
        echo "ERROR: Path does not exist: $path" >&2
        exit 1
    fi
done

echo "All paths exist."


# ---------------------------------------------------------------------------
# Run containers
# ---------------------------------------------------------------------------

validate_container() {
    local image="$1"
    local name="$2"
    local image_name

    echo
    echo "Checking: $image"

    image_name="$(singularity exec --containall "$image" \
        bash -c 'echo "$SINGULARITY_NAME"')"

    if [[ "$name" != "$image_name" ]]; then
        echo "ERROR: Container name mismatch." >&2
        echo "  Expected: $name" >&2
        echo "  Actual:   ${image_name:-<unset>}" >&2
        exit 1
    fi

    echo "  Container name: $image_name"
}

validate_container "$UV_PATH" "uv_r_cuda_0.2.1.sif"
validate_container "$CELLRANGER_PATH" "cellranger_9.0.0.sif"
validate_container "$STAR_PATH" "STAR_2.7.11b.sif"
validate_container "$CELLBENDER_PATH" "cellbender_0.3.0.sif"

echo
echo "All containers validated successfully."