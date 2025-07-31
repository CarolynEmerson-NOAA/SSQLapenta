#!/bin/bash


for i in $(seq 1 18); do
        dir_name="ENS_$i"
        mkdir -p "$dir_name"
        echo "Created directory: $dir_name"
done
