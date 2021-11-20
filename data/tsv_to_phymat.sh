#$1 tsv
#$2 phy mat
#works for fastme. not tested in other softwares.

tail -n +2 $1 | wc -l > $2
tail -n +2 $1 >> $2
sed -i "s/\t/  /g" $2

