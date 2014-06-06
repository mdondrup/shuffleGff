#!/usr/bin/env perl

use strict;
use warnings;
use RangeTree;
use Set::IntervalTree;
use Bio::Tools::GFF;
use Bio::SeqFeatureI;
use Getopt::Long;

my $excludef   = ''; ## gff3 file 
my $chromf = 'chromsizes.txt'; ## whitespace separated file
my $verbose = 0;
my $progress = 0;
my $outf = "";
my $chr_minlength = 1000; ## filter chromosmes

my $opts = GetOptions ("exclude=s"   => \$excludef,
		      "chromf=s" => \$chromf,
		       "minlength=i" => \$chr_minlength,
		       "out=s" => \$outf,
		       "progress"  => \$progress,
		       "verbose"  => \$verbose);

my @order = (); ## maintain the feature order of the original file


my $max_iter = 10000; ## maximum number of re-draws per feature
my $seed = srand();
print ("rand seed: $seed\n") if ($verbose);
print ("parsing chromosomes\n") if ($verbose);
my %chrom_sizes = parse_chroms($chromf); ## chrom sizes file
print ("read ". ( scalar (keys %chrom_sizes)). " chromosomes\n") if ($verbose);
#my %chrom_prob = _weight_to_dist(%chrom_sizes); ## we nee the probabilities for length bias draws
my ($chromtree, $max) = weight_to_interval(%chrom_sizes); ## we can do that faster

print ("parsing excludes \n") if ($verbose);
my $excludes =  ($excludef) ? parse_excludes($excludef) : undef; ## a RangeTree
print ("done parsing excludes\n") if ($verbose);

print ("parsing GFF \n") if ($verbose);
my %store = parse_gff($ARGV[0]); # store each to be processed ggf feature, with id as key, 
##  to avoid having to parse twice
print ("read". ( scalar keys %store). " GFF features \n") if ($verbose);

#### do the shuffling ######
print ("shuffling...\n") if ($verbose);
%store = shuffle_features(%store); 

write_gff($outf);

sub parse_chroms {
  my $file = shift;
  open my $fh, '<', $file or die "couldn't open $file: $!\n";  
  my %ret = ();
  while (<$fh>) {
    chomp;
    next if /^\s*$/; 
    my ($chrom, $size) = split /\t/;
    ## do not use chroms that are too small to host most genes: 
    next if $size < $chr_minlength;
    $ret{$chrom} = $size;
  }
  return %ret;
}

sub parse_excludes {
  my $tree = new RangeTree();
  my $gff = Bio::Tools::GFF->new(-file => $_[0], -gff_version => 3);
  while(my $f = $gff->next_feature()) {
    $tree->addOne($f->seq_id(), $f->start, $f->end); 
  }
  $gff->close;
  return $tree;
}

sub parse_gff {
  my %ret;
  my $gff = Bio::Tools::GFF->new(-file => $_[0], -gff_version => 3);
  while(my $f = $gff->next_feature()) {
    my $id = $f->primary_id || _new_id($f);
    $ret{$id} = $f;
    push @order, $id;
  }
  $gff->close;
  return %ret;
}

sub write_gff{
  my $out = Bio::Tools::GFF->new(-file => '>'.$_[0], -gff_version => 3);
  foreach (@order) {
    $out->write_feature($store{$_});
 }
}

sub shuffle_features{
  my %store = @_;
  my %shuffled = (); 
  my %offset = (); # store each processed parent **offset only**, with id as key;
  my $n = 0;
  my $stime = time;
  my $qcount = scalar (keys %store);

 ITEMSLEFT: while (scalar (keys %store) > 0) {
  ITER: while (my ($id, $f) = each %store) {
      if ($progress) {
	$n++;
	if (($n % 100) == 0 || $n == $qcount) {
	  my $perc = 100*$n/$qcount;
	  my $usec = time - $stime;
	  my $eta = ($qcount - $n)*($usec/$n)     /(60*60);
	  printf STDERR "procesing %d ( %.2f percent) of $qcount in %d sec, ETA %.2f h \r", $n, $perc, time - $stime, $eta ;
	}
      }




      if ($f->has_tag('Parent')) {
	### feature has a parent
	###  multi-parent features do not exist
	my ($parent) = $f->get_tag_values('Parent');
	### most likely, the parent feature has already been seen:
	if (exists $offset{$parent}) {
	  ### we can store the feature immediately
	  ### and if it has an ID itself store the offset in the offset hash
	  ### the offset can be applied directly to the feature
	  if (defined $offset{$parent}) {
	    ### do not do this for unplaced parent feature (indicated by undef), 
	    ### instead only remove this feature
	    $f->seq_id($shuffled{$parent}->seq_id());
	    $f->start($f->start + $offset{$parent});
	    $f->end($f->end + $offset{$parent});	   
	    $shuffled{$id} = $f;
	  }
	  $offset{$id} = $offset{$parent};
	  delete $store{$id};
	} else {
	  ### if the GFF file wasn't sorted by parent rank, we have to shuffle the parent first
	  ### do not do anything now, the parent will be picked up in the next iteration
	  ### many GFFs will only need one iteration
	  next ITER;
	}
      } else {
	#### this is a Top-level feature, shuffle coordinates:
	#### truts random function to return appropriate new start position
	my ($newseq_id, $newstart) = random_location_no_overlap_sparse($f);
	if (defined ($newseq_id)) { # feature could be placed
	  $f->seq_id($newseq_id);
	  $offset{$id} = $newstart - $f->start;
	  $f->start($f->start + $offset{$id});
	  $f->end($f->end + $offset{$id});
	  $shuffled{$id} = $f;
	 } else {
	   # feature couldn't be placed, indicated by undef hash value
	   $offset{$id} = undef;
	 }
	delete $store{$id};
	
      }
    } ### ITER: 
    print ("\n left with ". ( scalar keys %store). " features after iteration \n") if ($verbose);
    
  }

  return %shuffled;
}


### draw a random location, works best for sparse genomes, where |excludes| << |genome|
### for denser genomes it might be better to 
sub random_location_no_overlap_sparse {
  my $f = shift;
  my $seq_id = "";
  my $start = 0;
  ### draw sequence_id
  ### try max $max_iter times to draw a valid position
  ### total max tries is max_iter * max_chromiter
  my $c = 0;
 REPEAT:
  for ( $c = 0; $c < $max_iter; $c++) {
    while (1) {
      $seq_id = weighted_rand_interval($chromtree, $max);
      die "no overlap in drawing interval" unless defined $seq_id;
      last if $chrom_sizes{$seq_id} >= $f->length;
      printf "redraw chromosome, because %s is shorter (%d) than feature %s (%d)\n", 
	$seq_id,  $chrom_sizes{$seq_id}, $f->primary_id, $f->length
	if $verbose;
    }
    #### try max M*length times to place the start without overlap on this chromosome before drawing a new one
    for (my $d = 0; $d < 100; $d++ ) {
      
      $start = int (rand ($chrom_sizes{$seq_id} - $f->length - 1)) + 1;
      if (defined $excludes) {
	last REPEAT if ( $excludes->countOverlaps($seq_id, $start, $start+$f->length) == 0 );
	print "need to redraw ".$f->primary_id." (re-try $d-$c), because position $seq_id:$start overlaps exclude region\n" if $verbose;
      }    
    }
  }
  if ($c == $max_iter) {
    warn "tried $max_iter times to place feature ".$f->primary_id.", feature will be removed!\n";
    return;
  }
  return ($seq_id, $start);
}





### surrogate feature id to keep features in store that do not have an id
sub _new_id {
  my $f = shift;
  return (join ":", ($f->seq_id, $f->primary_tag, $f->start, $f->end, $f->strand));
}

########################################################################
### This is a faster implementation to draw biased random numbers for 
### integer weights, e.g. lengths. Convert all chromosomes into 
### successive intervals, then draw a single random interger for the 
### whole range and compute overlap with the prepared interval tree
########################################################################


sub weight_to_interval{
  my %weights = @_;
  my $tree = new Set::IntervalTree;
  my $max = 0;
  while (my ($k, $v) = each %weights) {
	my $start = $max + 1;
	my $end = $max + $v;
    $tree->insert($k, $max, $max+$v);
    $max += $v - 1;
  }
  return ($tree, $max);
}

sub weighted_rand_interval {
  my ($tree, $max) = @_;
  while (1) {
    my $pos = 1 + int (rand ($max-1));
    my $res = shift @{$tree->fetch($pos, $pos+1)}; 
    warn ("invalid interval drawn $pos (max: $max)") unless $res;	
  return $res if $res;
  }
}

__END__

######################################################################
### From Perl Cookbook, Chapter 2.10. Generating Biased Random Numbers
### This code is a bit inefficient, because the whole %weights keys 
### need to be searched. Better use Set::IntervalTree
### with this implementation, the drawing of the chromosome is the 
### most time consuming operation
######################################################################

# weight_to_dist: takes a hash mapping key to weight and returns
# a hash mapping key to probability
sub _weight_to_dist {
    my %weights = @_;
    my %dist    = ();
    my $total   = 0;
    my ($key, $weight);
    local $_;

    foreach (values %weights) {
        $total += $_;
    }

    while ( ($key, $weight) = each %weights ) {
        $dist{$key} = $weight/$total;
    }

    return %dist;
}

# weighted_rand: takes a hash mapping key to probability, and
# returns the corresponding element
sub _weighted_rand {
    my %dist = @_;
    my ($key, $weight);

    while (1) {                     # to avoid floating point inaccuracies
        my $rand = rand;
        while ( ($key, $weight) = each %dist ) {
            return $key if ($rand -= $weight) < 0;
        }
    }
}
