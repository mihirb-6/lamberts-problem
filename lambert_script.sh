#!/bin/bash

# Check if both arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <rust_file_search_term> <input_file_to_process>"
    exit 1
fi

SEARCH_TERM=$1
INPUT_FILE=$2

# Find the first Rust file matching the search term in the current directory or subdirectories
RUST_FILE=$(find . -type f -name "${SEARCH_TERM}" | head -n 1)

if [ -z "$RUST_FILE" ]; then
    echo "Error: No Rust file found matching '${SEARCH_TERM}'"
    exit 1
fi

echo "Found Rust file: $RUST_FILE"
echo "Compiling..."


# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Executing with input file: $INPUT_FILE"
    echo "----------------------------------------"
    # Run the compiled binary and pass the input file
    "$RUST_FILE" "$INPUT_FILE"
    echo "----------------------------------------"
else
    echo "Error: Compilation failed."
    exit 1
fi
