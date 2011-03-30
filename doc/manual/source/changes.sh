#!/bin/sh

files=`find . -type f`
for f in $files; do
    [ "$f" == "changes.sh" ] && continue
    [ "$f" == "movefiles.sh" ] && continue
    revs=`hg log --template '{rev} ' $f`
    revs="$revs 1"
    for r in $revs; do
	echo " --> $r"
    done
done