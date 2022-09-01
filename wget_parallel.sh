#!/bin/bash

set -euo pipefail

#This script takes a list of URLs and downloads 4 at a time in parallel.

sf=URLs_list.txt

cat $sf | xargs -P 4 -I{} wget {}

## wget --basement is an alternative if you want to leave it running in background and continue script
