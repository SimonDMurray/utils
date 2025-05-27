#!/bin/bash

# Can't set -e as if query fails due to permissions I want it to continue
set -uo pipefail

# Manifest hsa structure
# <file-name>\t<s3-path>
MANIFEST="${1:?'/path/to/manifest-files'}"
PARALLEL_NUMBER="${2:-'100'}"

s3_copy() {
  local line=$1

  # Set these variables to what you need
  AWS_PROFILE="<aws-profile-to-use>"
  OUTDIR="<s3://output-bucket/output-prefix"

  #timestamp is here to check parallelisation
  timestamp=$(date +%F_%T)

  file_name=$(echo ${line} | awk '{print $1}')
  file_path=$(echo ${line} | awk '{print $2}')

  outfile="${OUTDIR}/${file_name}"

  # stdout for sanity checks
  echo "${timestamp} ${file_name} ${file_path}"

  aws s3 cp --profile ${AWS_PROFILE} ${file_path} ${outfile}
}

# Export function for use with xargs
export -f s3_copy

index=1
# Only get the file length value from `wc -l`
MANIFEST_LENGTH=$(wc -l ${MANIFEST} | awk '{print $1}')
# Throttle number of jobs in parallel by doing a sequential iteration in chunks and then parallise each chunk
while [ "${index}" -le "${MANIFEST_LENGTH}" ]; do
    next_index=$(( index + PARALLEL_NUMBER ))

    # Create a subshell for isolation
    (
        count=0
        awk "NR >= ${index} && NR < ${next_index}" "${MANIFEST}" | while IFS= read -r line; do
            grep_s3_copy_info "$line" &  # launch in background
            count=$((count + 1))

            # Wait after reaching PARALLEL_NUMBER jobs in this chunk
            if [ "$count" -ge "${PARALLEL_NUMBER}" ]; then
                wait
                count=0
            fi
        done
        wait  # wait for remaining jobs in the chunk
    )

    index=${next_index}
done
