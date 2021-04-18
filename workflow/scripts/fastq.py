"""
fastq module to manipulate fastq files
Written by Michael Wolfe
"""

import gzip
from itertools import islice


class FastqEntry(object):
    """
    Stores information for a single fastq entry

    An example of a fastq entry is below

    @SEQ_ID
    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    +
    !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
    """
    def __init__(self, name = "", seq = "", opt = "+", qual = ""):
        self.name = name
        self.seq = seq
        self.opt = opt
        self.qual = qual

    def __str__(self):
        return "%s\n%s\n%s\n%s\n"%(self.name, self.seq, self.opt, self.qual)

class FastqFile(object):
    def __init__(self, infile = ""):
        self.infile = infile
        self.names = []
        self.data = {}

    def read_whole_file(self):
        if type(self.infile) is str and self.infile.endswith(".gz"):
            fhandle = gzip.open(self.infile, mode = "rt")
        else:
            fhandle = open(self.infile, mode = "r")
        # assumes well-formed fastq
        while True:
            entry = [line.rstrip() for line in islice(fhandle,4)]
            if entry:
                name, seq, opt, qual = entry
                self.data[name] = FastqEntry(name, seq, opt, qual)
                self.names.append(name)
            else:
                break
        fhandle.close()

    def stream_file(self):
        if type(self.infile) is str and self.infile.endswith(".gz"):
            fhandle = gzip.open(self.infile, mode = "rt")
        else:
            fhandle = open(self.infile, mode = "r")
        # this is a terrible way to do this, but I can't think of a better way
        while True:
            entry = [line.rstrip() for line in islice(fhandle,4)]
            if entry:
                name, seq, opt, qual = entry
                yield FastqEntry(name, seq, opt, qual)
            else:
                break

    def write(self, outfile):
        if type(outfile) is str and outfile.endswith(".gz"):
            fhandle = gzip.open(self.infile, mode = "wb")
        else:
            fhandle = open(self.infile, mode = "w")
        for entry in self:
            fhandle.write(str(entry))

    def __iter__(self):
        for name in self.names:
            yield self.data[name] 

    def __len__(self):
        return len(self.names)

    def names(self):
        return self.names

    def pull_entry(self, name):
        return self.data[name]

    def add_entry(self, entry):
        self.data[entry.name] = entry
        self.names.append(entry.name)

    def subset_by_names(self, names):
        out = FastqFile()
        for name in names:
            out.add_entry(self.data[entry.name])
        return out

if __name__ == "__main__":
    import sys
    infile = sys.argv[1]
    if infile == "0":
        infile = int(infile)

    outfile = sys.argv[2]
    if infile == "1":
        outfile = int(outfile)

    this_file = FastqFile(infile)
    for entry in this_file.stream_file():
        print(entry)
