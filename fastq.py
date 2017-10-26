#!/usr/bin/env python3
# This Python file uses the following encoding: utf-8
# -*- coding: utf-8 -*-
"""
Main class is Read, takes a file connection and extracts the next FASTQ format read.
"""

def mean(num_vec):
    """Take mean of a numeric vector."""
    return float(sum(num_vec)) / len(num_vec)


class Read:
    """Representation of a single read of a FASTQ file"""
    def __init__(self, conn):
        """Takes a file connection and returns a FASTQ structured read"""
        self.anno = conn.readline().strip()
        if not self.anno:
            raise EOFError
        self.seq = conn.readline().strip()
        self.strand = conn.readline().strip()
        self.quality_str = conn.readline().strip()

    def data(self):
        """Returns the read in FASTQ format"""
        return self.anno + '\n' + self.seq + '\n' + self.strand + '\n' + self.quality_str + '\n'

    def gc_content(self):
        """Returns the GC percentage in read"""
        cnt = {'N': 0, 'G': 0, 'C': 0, 'A': 0, 'T': 0}
        for c_ind in self.seq:
            cnt[c_ind] += 1
        return 100.0 * (cnt['G'] + cnt['C']) / sum(cnt.values())

    def base_content(self, base):
        """Returns the percentage of particular base"""
        cnt = {'N': 0, 'G': 0, 'C': 0, 'A': 0, 'T': 0}
        for c_ind in self.seq:
            cnt[c_ind] += 1
        return 100.0 * cnt[base] / sum(cnt.values())

    def base_count(self, base):
        """Returns the count of particular base"""
        cnt = {'N': 0, 'G': 0, 'C': 0, 'A': 0, 'T': 0}
        for c_ind in self.seq:
            cnt[c_ind] += 1
        return cnt[base]

    def valid(self):
        """Check if a read is valid"""
        # check strand info
        if self.strand not in ['+', '-']:
            print(('Invalid strand info: ' + self.strand))
            return False
        # check length of sequence is equal to length of quality
        elif len(self.seq) != len(self.quality_str):
            print('Invalid read: sequence and quality of different length')
            return False
        else:
            return True

    def quality(self, start=None, end=None):
        """Returns the quality string either from 'start' onwards or from 'start' to 'end'"""
        if start:
            if end:
                return [ord(x) - 33 for x in self.quality_str[(start - 1):end]]
            else:
                return [ord(x) - 33 for x in self.quality_str[:start]]
        else:
            return [ord(x) - 33 for x in self.quality_str[:start]]

    def avg_quality(self, start=None, end=None):
        """Returns the average quality either from 'start' onwards or from 'start' to 'end'"""
        if start:
            if end:
                return mean([ord(x) - 33 for x in self.quality_str[(start - 1):end]])
            else:
                return mean([ord(x) - 33 for x in self.quality_str[:start]])
        else:
            return mean([ord(x) - 33 for x in self.quality_str[:start]])


class Paired:
    """Representation of paired end read"""
    def __init__(self, read1, read2):
        self.read1 = Read(read1)
        self.read2 = Read(read2)
