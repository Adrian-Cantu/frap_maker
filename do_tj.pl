
#$ARGV[0]= path to folder with datasets
#$ARGV[1]= path to results folder that would be created

$run="mkdir -p $ARGV[1]";
system $run;

$run="ls $ARGV[0] | grep 'fasta' | cut -d '.' -f1 > $ARGV[1]/IDS.txt";
system $run;

my $x = "/IDS.txt";
my $path_to_ids = $ARGV[1].$x;

open my $handle, '<', $path_to_ids;
chomp(my @IDS = <$handle>);
close $handle;

my $tj='';

foreach(@IDS){
        $tj.=`echo -n "$_ "; grep -c ">" $ARGV[0]/$_.fasta`;
}

open OUT, '>' , "$ARGV[1]/tj.txt";
print OUT $tj;
close OUT;


