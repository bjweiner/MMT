#!/bin/sh
#
# argument = dirname
#
# 
case $# in
0)
  echo 'Usage: prune_files.sh <dirnames>'
  exit 0 ;;
# 1)
#   dir=$1 ;;
esac

dirs=$*

for dir in $dirs
do
cd $dir
if [ ! -d tmp_skycam ]; then
  mkdir tmp_skycam
  mv *skycam.fits tmp_skycam
fi

ls *.fits > list_all_fts.txt
awk '{print $1"[1]"}' list_all_fts.txt > list_all_fts_ext1.txt

pyraf <<EOF
hselect @list_all_fts_ext1.txt \$I,disperse yes > list_all_fts_ext1.grating.txt
.exit
EOF
sed 's/\[1\]//' < list_all_fts_ext1.grating.txt > list_all_fts.grating.txt

# cd lists
# for i in arc.list bias.list cal.list dark.list dflat.list sflat.list
for i in list_all_fts.grating.txt
do 
  grep 600_gpm $i | awk '{print $1}' > ztmp1
  if [ -s ztmp1 ]; then
    if [ ! -d tmp_600grating ]; then mkdir tmp_600grating ; fi
    mv `cat ztmp1` tmp_600grating
    cp $i $i.with600
    grep -v 600_gpm $i.with600 > $i
    echo Removed 600_gpm files from  $dir  $i
    cat ztmp1 >> list.with600line
  fi
done

# This won't catch sky flats since they don't get grating
# recorded in header. Or do they now?

cd ..
# finish loop over dirs
done


