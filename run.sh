#!/bin/sh
# wyong 2010.5.10
if [ $# -lt 2 ]; then
	echo "$0 1ank.pdb 214"
	exit
else
	pdb=$1
	numres=$2
fi

segid=A
lastres=`expr 1 + $numres`

# generate monomer contact list
./gen_csu_contacts_list.pl -p $pdb -i 1 -j $lastres -c $segid | awk "(\$2 < \$5 && (\$5-\$2) >= 4)" > .csu
awk "{ print \$2, \$5 }" .csu > list.dat

# calculate the alpha-carbon distances between residues in contact
./native_contacts_dist.pl -p $pdb -l list.dat -s $segid -m $lastres > contact_csu.dat
