import fasta as fa
import bed_utils as bed
import sys


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Combine fasta files into one fasta file")
    parser.add_argument('outfilepre', type=str, help="prefix for output files")
    parser.add_argument('infiles', type=str, nargs='+', help="files to combine")
    parser.add_argument('--masked_regions', type=str,
            help="A .bed file containing a list of regions to replace with Ns")
    
    args = parser.parse_args()
    genome_name = args.outfilepre
    fasta_files = args.infiles
    all_fastas = [fa.FastaFile() for fafile in fasta_files]
    # read in all fastas
    for fafile, fafilename in zip(all_fastas, fasta_files):
        with open(fafilename, mode = "r") as inf:
            fafile.read_whole_file(inf)

    final_fasta = fa.FastaFile()
    for fafile in all_fastas:
        for entry in fafile:
            final_fasta.add_entry(entry)
    lengths = {}
    for entry in final_fasta:
        chrm_name = entry.chrm_name()
        lengths[chrm_name] = len(entry)
    
    # Replace masked regions with Ns
    if args.masked_regions:
        masked_regions = bed.BedFile()
        masked_regions.from_bed_file(args.masked_regions)
        for entry in masked_regions:
            total_region_length = entry["end"] - entry["start"]
            this_chrm = final_fasta.pull_entry(entry["chrm"])
            this_chrm.mutate(entry["start"], entry["end"], "N"*total_region_length)

    # Figure out total size of each chromosome
    with open(genome_name + "_contig_sizes.tsv", mode = "w") as outf:
        for chrm, length in lengths.items():
            outf.write("%s\t%s\n"%(chrm, length))

    # figure out total mappable size of the genome
    with open(genome_name + "_mappable_size.txt", mode = "w") as outf:
        N_size = 0
        total_size = 0
        for entry in final_fasta:
            N_size += entry.seq.count("N")
            total_size += len(entry)
        outf.write("%i"%(total_size-N_size))
    # Write out final fasta
    with open(genome_name + ".fa", mode = "w") as outf:
        final_fasta.write(outf)




