#!perl -w
use strict;

my @arr;
my @parents;
my @kids;
my $line;
my $pheno;
my $geno;
my $pedid;

$line = <>;
while (defined ($line)) {
  @arr = split (" ", $line);
  @arr = @arr[9,10,11];
  $pheno = shift (@arr);
  $geno = join ("", sort(@arr));
  @parents = ($geno . $pheno);

  $line = <>;
  @arr = split (" ", $line);
  $pedid = $arr[0];
  @arr = @arr[9,10,11];
  $pheno = shift (@arr);
  $geno = join ("", sort(@arr));
  push (@parents, $geno . $pheno);

  @kids = ();
  while ($line = <>) {
      @arr = split (" ", $line);
      ($arr[0] != $pedid) and last;
      
      @arr = @arr[9,10,11];
      $pheno = shift (@arr);
      $geno = join ("", sort(@arr));
      push (@kids, $geno . $pheno);
  }

  print (join (' ', sort (@parents)), " ", join (' ', sort (@kids)), "\n");
}
