import bed_utils 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Combine bed files into one bed file. Removes duplicate entries and sorts by chrm, start")
    parser.add_argument('outfilepre', type=str, help="prefix for output files")
    parser.add_argument('infiles', type=str, nargs='+', help="files to combine")
    
    args = parser.parse_args()

    bed_files = args.infiles
    all_beds = [bed_utils.BedFile() for inf in bed_files]

    for bedobj, bed_file in zip(all_beds, bed_files):
        bedobj.from_bed_file(bed_file)
    final_bed = bed_utils.BedFile()
    for bedobj in all_beds:
        for entry in bedobj:
            final_bed.add_entry(entry)
    final_bed.cleanup()
    final_bed.sort()
    final_bed.write_bed_file(args.outfilepre + ".bed")



