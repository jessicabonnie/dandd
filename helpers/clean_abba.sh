#!/bin/bash

set -ex
filename=$1

#if [[ ! -f trio_progu_clean.csv ]] ; then
  sed -e 's/\/scratch16\/blangme2\/jessica\/data\/HVSVC2\/consensus\/allvar_//g' < ${filename} \
    | sed -e 's/\.fasta\.gz//g' \
    | sed 's/[" ]//g' | sed "s/[']//g" \
    | sed 's/\[//' | sed -e 's/\]//'
#fi