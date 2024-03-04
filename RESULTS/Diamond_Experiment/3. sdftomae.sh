#!/bin/bash

DIRNAME=$(dirname "$0")
echo "The script is located in: $DIRNAME"
cd "$DIRNAME"
for file in *.sdf; do
    base_name="${file%.*}"; # This strips the extension from the file name
    /shared/shared/schrodinger/2022-1/utilities/structconvert "$file" "$base_name.mae";
done
