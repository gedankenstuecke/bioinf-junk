# SSRScanner
This small script mainly wraps around the SSR-scanner MISA and gives back a simple csv-file that includes the SSRs and the sequence in which each SSR was found

## USAGE 
Just run the SSRScanner.py with the fasta-file you'd like to get scanned as the argument from the shell.
Example: "python SSRScanner.py input.fasta"

## CONFIGURATION
MISA can be configured by editing the misa.ini 
The first lines contains the SSR-definitions. Example: 2-5 means that dinucleotide-repeats should be returned if those repeats are found at least 5 times in a row. 
The second line contains the minimum number of non-SSR-nucleotides that need to be found before a new SSR may start
 
## DEPENDENCIES
In order to use MISA/the wrapper you need some tools installed:
* Perl http://www.perl.org 
* Python 2.x http://www.python.org/
* BioPython http://biopython.org/wiki/Main_Page
