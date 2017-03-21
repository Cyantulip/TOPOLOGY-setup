#!/usr/bin/perl -w

use strict;
use Getopt::Std;

# make sure command line argments okay
if ($#ARGV == -1) {
    &usage() and exit -1;
}

# defaults and global variables
my $pdb_file;
my $midres;
my $segid = " ";
my $nc_list_file;

# evaluate any flags that may change the default values
&getopts('tp:s:l:m:', \my %opts);

# input target pdb file name
if (defined $opts{'p'}) {
    $pdb_file = $opts{'p'};
} else {
    &usage() and exit;
}

if (defined $opts{'m'}) {
    $midres = $opts{'m'};
} else {
    &usage() and exit;
}

# segid
if (defined $opts{'s'}) {
    $segid = $opts{'s'};
}

# input native contact list
if (defined $opts{'l'}) {
    $nc_list_file = $opts{'l'};
} else {
    &usage() and exit;
}

# open pdb file to read in structure
open(PDB, "$pdb_file")
    or die "Unable to open $pdb_file.";

my $x = 0.0;
my $y = 0.0;
my $z = 0.0;
my @seginfo = ();
my $curr_seg;
my $i = 0;
my $resnum1="XXX";

# scan in the pdb file info for each segment
SCAN_PDB: while (my $line = <PDB>) {

    if ($line =~ /^ATOM/) {

        # use only natural proteins
	    my $atomtype = substr($line, 13, 3);
	    my $restype = substr($line, 17, 3);
	    my $resnum = substr($line, 23, 3);
	    $curr_seg = substr($line, 21, 1);
            if($resnum ne $resnum1){
               $i++;
            }

	    $x = substr($line, 30, 8);
	    $y = substr($line, 38, 8);
	    $z = substr($line, 46, 8);
#            print "ATOM",$resnum,$i,"\n";
            

	    if ($segid =~ /$curr_seg/i) {
		push(@seginfo, "$atomtype:$restype:$i:$x:$y:$z");
	    }

            $resnum1=$resnum;
    }
}

close(PDB);

# open native contact list file
open(NCL, "$nc_list_file")
    or die "Unable to open $nc_list_file.";

my @nc_list = ();
# scan in the native contact list file info
while (my $line = <NCL>) {
    unless ($line =~ /^#/) {
	(my $resnum1, my $resnum2) = split (/\s+/, $line);
#	print "$resnum1 $resnum2\n";
	$nc_list[$resnum1][$resnum2] = 1;
    }
}

my $num_contacts = 0;

my @first_select = ();
my @second_select = ();

# find contacts within monomer
# print "# Contacts within Monomer $segid\n";
@first_select = @seginfo;
@second_select = @seginfo;
$num_contacts = &find_contacts(\@first_select, \@second_select);
# print "# Number of contacts: $num_contacts\n";

exit 0; # end of "main"

sub calc_dist {
    my $ax = shift;
    my $ay = shift;
    my $az = shift;

    my $bx = shift;
    my $by = shift;
    my $bz = shift;

    my $dist;

    my $mx = $bx - $ax;
    my $my = $by - $ay;
    my $mz = $bz - $az;

    $dist = sqrt($mx*$mx + $my*$my + $mz*$mz);

    return $dist;
}

sub code_aa {
    my $aacode = shift;

    if ($aacode eq  "ALA") {
	return (1);
    }
    elsif ($aacode eq  "ARG") {
	return (2);
    }
    elsif ($aacode eq  "ASN") {
	return (3);
    }
    elsif ($aacode eq  "ASP") {
	return (4);
    }
    elsif ($aacode eq  "CYS") {
	return (5);
    }
    elsif ($aacode eq  "GLN") {
	return (6);
    }
    elsif ($aacode eq  "GLU") {
	return (7);
    }
    elsif ($aacode eq  "GLY") {
	return (8);
    }
    elsif ($aacode eq  "HIS") {
	return (9);
    }
    elsif ($aacode eq  "ILE") {
	return (10);
    }
    elsif ($aacode eq  "LEU") {
	return (11);
    }
    elsif ($aacode eq  "LYS") {
	return (12);
    }
    elsif ($aacode eq  "MET") {
	return (13);
    }
    elsif ($aacode eq  "PHE") {
	return (14);
    }
    elsif ($aacode eq  "PRO") {
	return (15);
    }
    elsif ($aacode eq  "SER") {
	return (16);
    }
    elsif ($aacode eq  "THR") {
	return (17);
    }
    elsif ($aacode eq  "TRP") {
	return (18);
    }
    elsif ($aacode eq  "TYR") {
	return (19);
    }
    elsif ($aacode eq  "VAL") {
	return (20);
    }
    else {
	return(0); # otherwise unknown
    }
}

sub is_calpha {
    my $atomtype = shift;

    if ($atomtype eq "CA ") {
	return (1);
    } else {
	return (0);
    }
}

sub is_sidechain {
    my $atomtype = shift;

    unless ($atomtype eq "N  "
	or $atomtype eq "CA "
	or $atomtype eq "C  "
	or $atomtype eq "O  ") {
	return (1);
    } else {
	return (0);
    }
}

sub find_contacts {
    my $first_listref = shift;
    my @first_list = @$first_listref;
    my $second_listref = shift;
    my @second_list = @$second_listref;

    my $num_contacts = 0;

    my $first_start = 0;
    my $first_end = @first_list;
    my $second_start = 0;
    my $second_end = @second_list;

    my $chainone=1;
    my $chaintwo=1;

    # check to see if same list
    my $same = 0;
    if (@first_list == @second_list) {
	$same = 1;
	for (my $i = 0; $i < @first_list; $i++) {
	    if ($first_list[$i] ne $second_list[$i]) {
		$same = 0;
		last;
	    }
	}
    }

    $first_end-- if ($same == 1);

    # calculate the ca pair distances
    my %ca_contact_dist = ();
    for (my $i = $first_start; $i < $first_end; $i++) {
	$second_start = ($i + 1) if ($same == 1);
	for (my $j = $second_start; $j < $second_end; $j++) {

	    (my $at_one, my $restype_one, my $resnum_one, 
	     my $x_one, my $y_one, my $z_one) = split(":", $first_list[$i]);

	    (my $at_two, my $restype_two, my $resnum_two,
	     my $x_two, my $y_two, my $z_two) = split(":", $second_list[$j]);

	    $resnum_one =~ s/\s//g;
	    $resnum_two =~ s/\s//g;
	    unless (defined ($nc_list[$resnum_one][$resnum_two])) {
		next;
	    }

	    my $dist = &calc_dist($x_one, $y_one, $z_one, $x_two, $y_two, $z_two);

#            print "$resnum_one $resnum_two $dist\n";
	    if (&is_calpha($at_one) and &is_calpha($at_two)) {
		$ca_contact_dist{"$restype_one$resnum_one:$restype_two$resnum_two"} = $dist;
		$nc_list[$resnum_one][$resnum_two] = $dist;
		$num_contacts++;
#	   	printf("%4d %4d %.4f\n", $resnum_one,$resnum_two,$dist);
		if ( $resnum_one <= $midres ) {
			$chainone=1; 
		} else {
			$chainone=2; 
		}
                if ( $resnum_two <= $midres ) {
                        $chaintwo=1;
                } else {
                        $chaintwo=2;
                }
                printf("%1d %4d %1d %4d %.4f\n", $chainone,$resnum_one,$chaintwo,$resnum_two,$dist);
	    }
	}
    }

    return $num_contacts;
}

sub usage {
    print STDERR << "EOF";

usage: $0 -p pdb_file [-s segid] -l native_contact_list -m midres

 -p       : protein PDB file
 -s       : segment id of interest
 -m       : the mid-residue number 
 -l       : native contact list

example: $0 -p SUC1.pdb -s A -l native_contact_list.dat -m 192

EOF
}

