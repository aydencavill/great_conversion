#!/bin/bash

example_series=$1
new_src_dir=$2
temp_file=$(mktemp)
jq ".series.src_dir = \"$new_src_dir\"" "$example_series" > "$temp_file"
mv $temp_file $example_series
echo $1
echo $new_src_dir
