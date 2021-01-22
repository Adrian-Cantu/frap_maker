#!/usr/bin/perl
#
use Getopt::Long;
use strict;
use warnings;

my $id;
my $hits_only;
my $norm_only;
my %desc;
my %Tj;
my $tj_file;
my $millon;
my $lmean=50000;

GetOptions ("f=s"     => \$id,
	    "h"       => \$hits_only,	
	    "n"       => \$norm_only,    
	    "t=s"       => \$tj_file,
 	    "m"       => \$millon,
		"l=i"	=> \$lmean,
);  # flag

my @glob_ARGV=@ARGV;
my %list;
my %length;
my %sum;
my %norm;
my %samples;
my %sum_row;
my @order;
open IDS , '<' , $id;

#start hashs
foreach (@ARGV){
	$list{$_}={};
	$norm{$_}={};
	$sum{$_}=0;
}

#read the fasta and extract id, length and description for each sequence
my $ll;
my $ll_name;
while(<IDS>) {
	if (/^>/) {
		chomp;
		my @f=split;
		my $id = shift @f;
		$id =~ s/^>//;
		$desc{$id}=join(' ',@f);
		push @order,$id;

		$length{$ll_name}=length($ll) if $ll_name;
		$ll='';
		$ll_name=$id;
		foreach (@ARGV){
        		$list{$_}->{$id}=0;
		}
	}
	else {
		chomp;
		$ll=$ll.$_;
	}
}
$length{$ll_name}=length($ll);


my %filehandlers;
foreach(@ARGV) {
	chomp;
	open $filehandlers{$_}, '<' , $_  or die "Can't open $_ for output: $!";
	my $filename = (split(/\//,$_))[-1];
	my $sample_id = (split(/\./,$filename))[1];
	$samples{$_} = $sample_id;	
}


open TJ , '<', $tj_file;
while(<TJ>) {
	chomp;
	my @f=split;
	$Tj{$f[0]}=$f[1];
}


foreach(@ARGV) {
	my $f_id=$_;
	while(readline($filehandlers{$f_id})) {
	my @fields=split;
	$list{$f_id}->{$fields[1]}=+$fields[0];
	}
}
foreach(@glob_ARGV) {
	my $f_id=$_;
	my $total=0;
	foreach(@order) {
		$total=$total+$list{$f_id}->{$_};
	}
	$sum{$f_id}=$total;
}
		
#foreach(@glob_ARGV) {
#	my $f_id=$_;
#    my $local_sum=0;
#	foreach my $g_id (@order) {
#        my $normalized_score=($list{$f_id}->{$_}/$Tj{$samples{$f_id}})*($lmean/$length{$_});
#		$norm{$f_id}->{$_}=$normalized_score;
#        $local_sum+=$normalized_score;
#	}
#    $sum_row{$g_id}=$local_sum
#}

foreach my $g_id (@order) {
    my $local_sum=0;
    foreach my $f_id (@glob_ARGV) {
        my $normalized_score=($list{$f_id}->{$g_id}/$Tj{$samples{$f_id}})*($lmean/$length{$g_id});
        $norm{$f_id}->{$g_id}=$normalized_score;
        $local_sum+=$normalized_score;
    }
    $sum_row{$g_id}=$local_sum
}    

print "\t";
foreach(@glob_ARGV) {
	print "$samples{$_}";
	unless(  \$_ == \$glob_ARGV[-1]  ) {
        print "\t";
    	}
}
print "\tid\n";

#foreach(@order) {
foreach my $g_id (sort { $sum_row{$b} <=> $sum_row{$a} } keys %sum_row) {
 #   my $g_id=$_;
    next unless ($sum_row{$g_id}>0);
	print "$desc{$g_id}";
	foreach my $f_id (@glob_ARGV) {
		if ($hits_only) {
			print "\t$list{$f_id}->{$g_id}";
		} elsif ($norm_only) {
			print "\t$norm{$f_id}->{$g_id}";
		} elsif ($millon) {
			my $t_kk = $norm{$f_id}->{$g_id} * 1000000;
		       print "\t$t_kk";
	       }	       
		
	}
	print "\t$g_id\n";
}


