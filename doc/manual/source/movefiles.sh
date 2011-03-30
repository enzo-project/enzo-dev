#!/bin/sh

dir=doc/manual/source

src=/Users/jwise/codes/EnzoWOC/enzo-gc
dest=/Users/jwise/codes/EnzoWOC/copy

cd ${src}/${dir}
find . -name '*~' -exec rm {} \;
src_dirs=`find . -type d`
src_files=`find . -type f`

cd ${dest}/${dir}
find . -name '*~' -exec rm {} \;
find . -type d > /tmp/dest_dirs.txt
find . -type f > /tmp/dest_files.txt

for d in $src_dirs; do
    count=`grep -c $d /tmp/dest_dirs.txt`
    if [ $count -eq 0 ]; then
	mkdir -vp ${dest}/${dir}/$d
    fi
done

cd ${dest}/${dir}
for f in $src_files; do
    count=`grep -c $f /tmp/dest_files.txt`
    if [ $count -eq 0 ]; then
	cp -v ${src}/${dir}/$f ${dest}/${dir}/$f
	hg add $f
    fi
done

