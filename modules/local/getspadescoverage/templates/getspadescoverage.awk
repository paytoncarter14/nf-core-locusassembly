#!/usr/bin/env awk
BEGIN {
    FS = ":" # Set field separator to colon for initial split
    print "locus,length,coverage" # Print CSV header
}
/^>/ {
    # Extract the name before the first colon
    name = \$1
    sub(/^>/, "", name) # Remove the leading '>'

    # Extract length
    # Look for "length_DIGITS" anywhere in the header
    if (match(\$0, /length_([0-9]+)/)) {
        length_val = substr(\$0, RSTART + 7, RLENGTH - 7)
    } else {
        length_val = "N/A"
    }

    # Extract coverage
    # Look for "cov_DIGITS.DIGITS" anywhere in the header
    if (match(\$0, /cov_([0-9.]+)/)) {
        coverage_val = substr(\$0, RSTART + 4, RLENGTH - 4)
    } else {
        coverage_val = "N/A"
    }

    print name "," length_val "," coverage_val
}