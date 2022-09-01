#!/bin/bash

set -euo pipefail

#This takes all .bz2 files in a directory and converts them to .gz in parallel

parallel ‘bzcat {} | gzip -c > {.}.gz’ ::: *bz2
