#!/bin/bash

DIRNAME=$(dirname "$0")
echo "The script is located in: $DIRNAME"
cd "$DIRNAME"
for file in *-out.mae; do
    base_name="${file%.*}-mini"; # This strips the extension from the file name and add some modification
    /shared/shared/schrodinger/2022-1/utilities/structconvert "$file" -split-nstructures 1 "$base_name.mae";
done
