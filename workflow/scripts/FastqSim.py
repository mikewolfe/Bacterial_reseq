import sys
import numpy as np
import fasta as fa
import fastq as fq
import bed_utils as bed
import gzip

class QualityDistro(object):
    def __init__(self, max_length = int(5e6)):
        self.center_per_pos = np.zeros(max_length, dtype=int)
        self.spread_per_pos = np.zeros(max_length, dtype=float)

    def generate_quality(self, length, rng):
        quals = np.zeros(length, dtype=int)
        for loc, center, spread in zip(np.arange(length), self.center_per_pos, self.spread_per_pos):
            quals[loc] = int(rng.normal(center, spread))
        return "".join([chr(val + 33) for val in quals])

    def uniform_distro(self, center, spread):
        self.center_per_pos[:] = center
        self.spread_per_pos[:] = spread

class FragmentLengthDistro(object):
    def __init__(self, center, spread):
        self.center = center
        self.spread = spread

    def generate_length(self, rng, max_size = int(5e6)):
        return min(max_size, int(rng.normal(self.center, self.spread)))

class ReadLengthDistro(object):
    def __init__(self, center, spread):
        self.center = center
        self.spread = spread

    def generate_length(self, rng, max_size = int(5e6)):
        return min(max_size, int(rng.normal(self.center, self.spread)))

class ReadSampler(object):
    def __init__(self, length_distro, fragment_distro, quality_distro, fasta):
        self.length_distro = length_distro
        self.fragment_distro = fragment_distro
        self.quality_distro = quality_distro
        self.fasta = fasta

    def generate_read(self, name, chrm, loc, strand, rng, full_length = None):
        name = "@%s"%(name)
        if full_length is not None:
            length = full_length
        else:
            length = self.length_distro.generate_length(rng)
        quality = self.quality_distro.generate_quality(length, rng)
        if strand == "-":
            # have to shift register by one to match properly
            seq = self.fasta.pull_entry(chrm).pull_seq(loc-length+1, loc+1, rc = True, circ = True)
        else:
            seq = self.fasta.pull_entry(chrm).pull_seq(loc, loc+length, rc = False, circ = True)

        return fq.FastqEntry(name = name, seq= seq, qual = quality)

    def generate_singleend(self, chrm, loc, rng):
        fragment_length = self.fragment_distro.generate_length(rng)
        fragment_strand = rng.choice(["-", "+"])
        five_loc = int(loc - fragment_length/2)
        three_loc = int(loc + fragment_length/2)
        name = "%s_%s_%s"%(chrm, loc, fragment_strand)
        if fragment_strand == "-":
            R0 = self.generate_read(name, chrm, three_loc, fragment_strand, rng, full_length = fragment_length)
        else:
            R0 = self.generate_read(name, chrm, five_loc, fragment_strand, rng, full_length = fragment_length)
        return R0

    def generate_read_pair(self, chrm, loc, rng):

        fragment_length = self.fragment_distro.generate_length(rng)
        fragment_strand = rng.choice(["-", "+"])
        five_loc = int(loc - fragment_length/2)
        three_loc = int(loc + fragment_length/2)
        name = "%s_%s_%s"%(chrm, loc, fragment_strand)
        if fragment_strand == "-":
            R1 = self.generate_read(name, chrm, three_loc, fragment_strand, rng)
            R2 = self.generate_read(name, chrm, five_loc, "+", rng)
        else:
            R1 = self.generate_read(name, chrm, five_loc, fragment_strand, rng)
            R2 = self.generate_read(name, chrm, three_loc, "-", rng)
        return (R1, R2)
        
class LocSampler(object):

    def __init__(self, FastaFile = None):
        self.weights = {}
        if FastaFile is not None:
            for entry in FastaFile:
                self.weights[entry.chrm_name()] = np.ones(len(entry), dtype=float)
        self.probs = {}
        self.determine_probs()
        self.chrm_probs = {}
        self.determine_chrm_probs()

    def determine_chrm_probs(self):
        values = []
        keys = []
        for key, value in self.weights.items():
            keys.append(key)
            values.append(np.sum(value))
        total = np.sum(values)
        values = values / total
        for key, value in zip(keys, values):
            self.chrm_probs[key] = value

    def determine_probs(self):
        for key, value in self.weights.items():
            self.probs[key] = value/np.sum(value)

    def add_enrichment_by_location(self, bedfile, footprint, max_enrichment = 16):
        values = []
        for entry in bedfile:
            values.append(entry["score"])
        # scale values between 1 and max enrichment
        values = np.array(values)
        values = (max_enrichment - 1) *\
                ((values - np.min(values))/\
                (np.max(values)-np.min(values))) + 1

        half_footprint = int(footprint/2)
        for value, entry in zip(values, bedfile):
            center = int((entry["end"] - entry["start"])/2 + entry["start"])
            these_weights = self.weights[entry["chrm"]]
            these_weights[(center-half_footprint):(center + half_footprint)] += value
            self.weights[entry["chrm"]] = these_weights
        self.determine_probs()
        self.determine_chrm_probs()

    def simulate(self, n, rng):
        locs = {}
        # which chromosomes to pull from
        possible_chrms = list(self.chrm_probs.keys())
        chrm_probs = list(self.chrm_probs.values())
        chosen_chrms = rng.choice(possible_chrms, n, p = chrm_probs)
        chrm, chrm_n = np.unique(chosen_chrms, return_counts = True)
        # which locations to pull from for each chromosome
        for this_chrm, this_chrm_n in zip(chrm, chrm_n):
            locs[this_chrm] = rng.choice(np.arange(len(self.probs[this_chrm])), 
                    this_chrm_n, p = self.probs[this_chrm])
        return locs


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Simulate fastq reads from a genome")
    parser.add_argument("infile", type=str, help="path to input fasta genome file")
    parser.add_argument("out_pre", type =str, help = "output prefix for files")
    parser.add_argument("--num_fragments", type = int, 
            help = "Number of fragments to simulate. Default = 1e6",
            default = 1000000)
    parser.add_argument("--rng_seed", type = int, help = "Seed for sampling")
    parser.add_argument("--peaks", type = str, 
            help = "Bed file specifying areas of enrichment")
    parser.add_argument("--max_enrichment", type = int, 
            help="Maximum fold enrichment above background", default = 16)
    parser.add_argument("--footprint", type = int,
            help = "Estimated size of protein footprint for binding", default = 150)
    parser.add_argument("--quality_mean", type = int,
            help = "Average quality score per base (Default = 36)", default = 36)
    parser.add_argument("--quality_std", type = int,
            help = "Std quality score per base (Default = 1)", default = 1)
    parser.add_argument("--fragment_length_mean", type = int,
            help = "Average fragment size (Default = 400)", default = 400)
    parser.add_argument("--fragment_length_std", type = int,
            help = "Std fragment size (Default = 10)", default = 10)
    parser.add_argument("--read_length_mean", type = int,
            help = "Average read length (Default = 30)", default = 30)
    parser.add_argument("--read_length_std", type = int,
            help = "Std read length (Default = 0.0001)", default = 0.0001)
    parser.add_argument("--single_end", action = "store_true",
            help = "Generate single end data rather than pairs")
    args = parser.parse_args()
    # input genome fasta
    infasta = args.infile
    # yaml file controlling the software
    out_prefix = args.out_pre
    n_fragments = args.num_fragments

    genome = fa.FastaFile()

    rng = np.random.default_rng(args.rng_seed)

    with open(infasta, mode = "r") as inf:
        genome.read_whole_file(inf)
    locations = LocSampler(genome)

    # input locations for true peaks
    if args.peaks:
        inlocs = args.peaks

        inbed = bed.BedFile()
        inbed.from_bed_file(inlocs)

        enrichment = args.max_enrichment

        locations.add_enrichment_by_location(inbed, args.footprint, max_enrichment = enrichment)
    qual_distro = QualityDistro()
    qual_distro.uniform_distro(args.quality_mean, args.quality_std)
    fragment_distro = FragmentLengthDistro(args.fragment_length_mean, args.fragment_length_std)
    read_length_distro = ReadLengthDistro(args.read_length_mean, args.read_length_std)

    reads = ReadSampler(read_length_distro, fragment_distro, qual_distro, genome)
    
    all_locs = locations.simulate(n_fragments, rng)
    if args.single_end:
        out_R0 = gzip.open(out_prefix + "_R0.fastq.gz", mode = "wb")
        for key in all_locs.keys():
            for loc in all_locs[key]:
                R0 = reads.generate_singleend(key, loc, rng) 
                out_R0.write(str(R0).encode())
        out_R0.close()
    else:
        out_R1 = gzip.open(out_prefix + "_R1.fastq.gz", mode ="wb")
        out_R2 = gzip.open(out_prefix + "_R2.fastq.gz", mode = "wb")
    
        for key in all_locs.keys():
            for loc in all_locs[key]:
                R1, R2 = reads.generate_read_pair(key, loc, rng) 
                out_R1.write(str(R1).encode())
                out_R2.write(str(R2).encode())
        out_R1.close()
        out_R2.close()
