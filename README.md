# Paired end FASTQ splitting/demultiplexing with specific barcode positions

## Description
This script is used for demultiplexing paired end FASTQ files based on a list of known barcodes and known barcode coordinates. This was designed to work with the CEL-seq2 single cell protocol read structure.

This script may not work on macOS due to a limit in how many files can be open. If you encounter this problem try using iTerm2 or increasing the opened files limit.

## Installation
To install these scripts, clone this repository and link the scripts to a location in the system PATH.
```
git clone https://github.com/Shians/paired_fastq_split.git
cd paired_fastq_split
ln -s paired_fastq_split.py /usr/local/bin/ # Change /usr/local/bin as required
```
Do not delete or move scripts after linking.

## Usage
Help documentation can be found using
```
paired_fastq_split.py -h
```

which returns the call signature

```
usage: paired_fastq_split.py [-h] -file1 FILE1 -file2 FILE2 -barcodes BARCODES
                      [-bc_pos BC_POS BC_POS] -prefix PREFIX [-r2bc]
```

A basic call looks like

```
paired_fastq_split.py -file1 exp_R1.fastq.gz -file2 exp_R2.fastq.gz -barcodes barcodes.txt -prefix exp
```

by default we assume that the barcodes are on read 1 in positions 7 to 16 (length 8). The output will be a pair of files for each barcode with the filename `{prefix}_{barcode}_R(1/2).fastq` along with a pair of files for reads that did not have an exact match to any barcode.

The barcode file should be plaintext with a unique barcode on each line

```
GTAGCTCA
TAGACCAC
GGTCTATG
GTCCGAAT
...
```
