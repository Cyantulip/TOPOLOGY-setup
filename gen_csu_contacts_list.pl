#!/usr/bin/perl -w

use strict;
use Getopt::Std;

# make sure command line argments okay
if ($#ARGV == -1) {
    &usage() and exit -1;
}

# defaults and global variables
my $exec = "./resc";
my $pdb_file;
my $frst_res;
my $last_res;
my $chain = " ";

# evaluate any flags that may change the default values
&getopts('p:i:j:c:', \my %opts);


# pdb filename
if (defined $opts{'p'}) {
    $pdb_file = $opts{'p'};
} else {
    &usage() and exit;
}

if (defined $opts{'i'}) {
    $frst_res = $opts{'i'};
} else {
    &usage() and exit;
}

if (defined $opts{'j'}) {
    $last_res = $opts{'j'};
    &usage() and exit 
	unless ($last_res > $frst_res);
}
$last_res = $frst_res unless (defined $last_res);

if (defined $opts{'c'}) {
    $chain = $opts{'c'};
}


for (my $a = $frst_res; $a <= $last_res; $a++) {

        #printf("$exec $pdb_file $a $chain \n");
    # run CSU program
    open(CSU, "$exec $pdb_file $a $chain |")
	or die "Unable to run $exec.";

    # scan in the lines of output
    while (my $cline = <CSU>) {
	last if ($cline =~ /^Table II /);
    }

    while (my $cline = <CSU>) {
	next if ($cline =~ /\-{2,}/);
	next if ($cline =~ /Specific/);
	next if ($cline =~ /Residue/);
	next if ($cline =~ /^$/);
	last if ($cline =~ /^Table III /);
	print "$cline";
    }
    close (CSU);
}


exit 0; # end of "main"

sub usage {
    print STDERR << "EOF";

usage: $0 -p pdb_file -i first_resnum -j second_resnum [-c segid]

 -p       : input pdb file
 -i       : first residue number
 -j       : second residue number
 -c       : segment id

example: $0 -p 1I0C.pdb -i 6 -j 65 -c A

EOF
}

