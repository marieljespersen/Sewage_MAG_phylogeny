#!/bin/bash

in=$1
out=$2

perl -ne '$s = $_ =~ s/S(\w\w\w\d+)C\d+/$1/g; print $_, "\n"' $in > $out

