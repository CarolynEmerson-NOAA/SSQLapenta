#!/bin/bash

PYTHON_SCRIPT="generate_urls.py"

for dir in ENS_*; do
    if [ -d "$dir" ]; then
        num=$(echo "$dir" | grep -oE '[0-9]+')
        echo "🔍 Processing $dir (ensemble number $num)"

        # Run Python script with ensemble number and output folder
        python "$PYTHON_SCRIPT" "$num" "$dir"

        links_file="$dir/wofs_file_list.txt"

        if [ -f "$links_file" ]; then
            echo "⬇️  Downloading files listed in $links_file"
            wget -i "$links_file" -P "$dir" --user-agent=Carolyn.Emerson

            echo "Deleting $links_file"
            rm "$links_file"
        else
            echo "⚠️  No links file found for $dir"
        fi
    fi
done
