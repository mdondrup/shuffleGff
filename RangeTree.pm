package RangeTree;

use strict;
use warnings;

use Set::IntervalTree;

sub new {
  my $class = shift;  
  my $self = {}; # container is a hash with seqnames as keys
  bless $self, $class;
}

sub addOne {
  my ($self, $seqname, $start, $end) = @_;
  $self->{$seqname} = new Set::IntervalTree unless (exists $self->{$seqname}); 
  ## Set::IntervalTree interval is half open [,) while GFF interval is closed [,]
  $self->{$seqname}->insert({}, $start, $end + 1 ); ## add 1 to the end 
  return $self;
}

sub add {
  my ($self, $seqname, $id, $start, $end, $width, $strand, $meta) = @_;
  die "start and seqname are required" unless $start && $seqname;
  die "either end or width needs to be defined" unless $end || $width;
  die "$start < $end" if  $end && $start >= $end;
  $strand ||= '*';
  $end ||= $start + $width;  
  $width = $end - $start;
  $self->{$seqname} = new Set::IntervalTree unless (exists $self->{$seqname});     
  $self->{$seqname}->insert({seqname=>$seqname, id=>$id, 
				    start=>$start, end=>$end, 
				    width=>$width, strand=>$strand, meta=>$meta}, 
				   $start, $end + 1);
  return $self;
}

sub findOverlaps {
  my ($self, $seqname, $start, $end) = @_;
  return [] unless exists $self->{$seqname};
  return $self->{$seqname}->fetch($start, $end + 1);
}

sub countOverlaps {
  my ($self, $seqname, $start, $end) = @_;
  return 0 unless exists $self->{$seqname};
  return scalar @{$self->{$seqname}->fetch($start, $end + 1)};
}

1;
