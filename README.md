shuffleGff
=============

shuffleGff -- a Perl script to shuffle genomic features leaving the gene models intact

Usage
------------------------------------------
   
   shuffleRanges.pl -chromsizes chromsizes.txt -excludes xclude.gff3 -out shuffled.gff3  -minlength 1000 
   -progress -verbose input.gff3

Dependencies
------------------------------------------
   * Set::IntervalTree;
   * BioPerl;
