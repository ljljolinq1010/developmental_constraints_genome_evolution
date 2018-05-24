use strict;
use warnings;
## use it to generate elegans_GeneID_sequenceName.txt
## perl ID_conversion.pl  elegans_GeneID_transcriptID.txt elegans_GeneID_sequenceName.txt
my $infile = shift @ARGV;
my $outPutFile = $ARGV[0];
#open a file handle IN to save this file
open (IN,"<",$infile) || die "can not open a file: $!\n";
open (OUT, ">", $outPutFile) or die "can not open a file: $!\n";
while (<IN>) {
  #delete \n in each line of this file
  chomp;
  #split this file according to space
  my @array1=split(/\t/, $_);
  my @array2=split(/\./, $array1[1]);
  print $array1[1],"\n";
  print $array2[1],"\n";

  # remove the english letter for $array2[1]
  my $lastLetter= substr $array2[1],-1;
  print $lastLetter,"\n";
  if ($lastLetter =~ /^[a-zA-Z]+$/ ) {
    chop $array2[1];
  }
  print $array2[1],"\n";

  my $transcriptID=join('.',$array2[0],$array2[1]);
  print $transcriptID,"\n";
  print OUT "$array1[0]","\t","$transcriptID","\n";
}
close IN;
close OUT;
