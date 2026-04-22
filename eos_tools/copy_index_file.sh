#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <source_file> <target_directory>"
    exit 1
fi

SOURCE_FILE="$1"
TARGET_DIR="$2"

# Check if source file exists
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file '$SOURCE_FILE' does not exist."
    exit 1
fi

# Check if target directory exists
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Target directory '$TARGET_DIR' does not exist."
    exit 1
fi

# Find all subdirectories and copy the file into them
find "$TARGET_DIR" -type d -exec cp "$SOURCE_FILE" {} \;

echo "File '$SOURCE_FILE' copied to all subdirectories of '$TARGET_DIR'"