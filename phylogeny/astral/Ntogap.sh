#§/bin/bash

in=$1
out=$2

cat $in | tr 'N' '-' > $out
