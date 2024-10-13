#!/bin/bash

filename=$1

# sort lines that do not contain underscore
grep -v '_' $filename | sort > /tmp/sorted_no_underscore

# sort lines that contain underscore
grep '_' $filename | sort > /tmp/sorted_underscore

# concatenate the two sorted files
cat /tmp/sorted_no_underscore /tmp/sorted_underscore > $filename

# remove the temporary files
rm /tmp/sorted_no_underscore /tmp/sorted_underscore

