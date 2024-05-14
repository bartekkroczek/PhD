#!/bin/bash

# Temporary file suffixes
tmp_suffix=".utf8"
tmp_cleaned=".cleaned"

# Loop through all .csv files in the current directory
for file in *.csv; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Convert Windows CRLF line endings to Unix LF line endings
        dos2unix "$file"

        # Convert file to UTF-8 and output to a temporary file
        iconv -f $(file -bi "$file" | sed -e 's/.*charset=//') -t utf-8 "$file" > "${file}${tmp_suffix}"
        
        # Remove empty lines and lines with only whitespace from the converted file and output to another temporary file
        sed '/^\s*$/d' "${file}${tmp_suffix}" > "${file}${tmp_cleaned}"
        
        # Replace the original file with the cleaned, converted file
        mv "${file}${tmp_cleaned}" "$file"
        
        # Clean up the intermediate UTF-8 conversion file
        rm -f "${file}${tmp_suffix}"
    fi
done

