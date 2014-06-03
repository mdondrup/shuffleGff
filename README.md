shuffleGff
=============

shuffleGff -- a Perl script to shuffle genomic features leaving the gene models intact

Usage
------------------------------------------
   
   shuffleRanges.pl -chromsises chromsizes.txt -excludes xclude.gff3 -out shuffled.gff3  -minlength 1000 
   -progress -verbose input.gff3

Dependencies
------------------------------------------
   use Set::IntervalTree;
   use Bio::Tools::GFF;
   use Bio::SeqFeatureI;
   use Getopt::Long;
