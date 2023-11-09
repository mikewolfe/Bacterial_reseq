from Bio import SeqIO
import sys
import argparse

#
#gi|48994873|gb|U00096.2|mod|ATCC.47076| - 1..4639676
#233 RNAs
#Location     

# example ptt file
# Pseudomonas aeruginosa PAO1 chromosome, complete genome. - 0..6264404
# 5572 proteins
# Location    Strand  Length  PID Gene    Synonym Code    COG Product
# 483..2027   +   514 -   dnaA    PA0001  -   chromosome replication initiator DnaA
# 2056..3159  +   367 -   dnaN    PA0002  -   DNA polymerase III subunit beta
# 3169..4278  +   369 -   recF    PA0003  -   DNA replication and repair protein RecF
# 4275..6695  +   806 -   gyrB    PA0004  -   DNA gyrase subunit B
# 7018..7791  -   257 -   lptA    PA0005  -   lysophosphatidic acid acyltransferase
# 7803..8339  -   178 -   -   PA0006  -   D-glycero-beta-D-manno-heptose-1,7-bisphosphate 7-phosphatase
# 8671..10377 +   568 -   -   PA0007  -   hypothetical protein
# example rnt file
# Bacteria name - 1..6264404
# 106 RNAs (number sRNAs)
# Location Strand Length PID Gene Synonym Code COG Product
# 298816..298892 - 77 110645304 - PA0263.1 - - Arg tRNA

def convert_strandval(strandval):
    if strandval == 1:
        out = "+"
    elif strandval == -1:
        out = "-"
    else:
        out = "NA"
    return out

def check_partial_start(coord):
    return "<" in str(coord)

def check_partial_end(coord):
    return ">" in str(coord)



def out_ptt_rnt(gb, ftype = "ptt", chrm = None, qual_names = ["locus_tag"]):
    # rockhopper only parses the genome name if it has >=5 fields split on
    # pipes. The 0-indexed [3] field is the name given to *_transcripts when
    # split on a "."
    # The [4] field is another name which I don't know where it goes
    if chrm:
        outchrm = str(chrm)
    else:

        outchrm = str(gb.id)
    if ftype == "rnt":
        features = ["tRNA", "rRNA", "ncRNA"]
        name =  "RNAs"
        unknown_name = "NA_RNA"
        length_func = lambda start, end: end - start
    elif ftype == "ptt":
        features = ["CDS"]
        name = "proteins"
        unknown_name = "NA_Prot"
        length_func = lambda start, end: int(((end - start)/3) - 1)
    else:
        raise ValueError("Filetype %s not recognized"%ftype)

    all_features = []
    ptt_rnt_name = set()
    for feature in gb.features:
        outstr = ""
        if feature.type in features:
        # in one-based coordinates need to adjust
            outstr += str(feature.location.start + 1) + ".." + str(feature.location.end) + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\t"
            outstr += str(length_func(feature.location.start,feature.location.end)) + "\t"
            outstr += str(feature.qualifiers.get("protein_id", ["-"])[0]) + "\t"
            # often times things need a gene name, thus will replace with locus tag or
            # if that doesn't exist. NA_number
            unique_name = check_qualifiers(feature, qual_names, unknown_name = unknown_name)
            num = 1
            id_name = unique_name
            while id_name in ptt_rnt_name:
                id_name = unique_name + "_%d"%(num)
                num += 1
            ptt_rnt_name.add(id_name)
            outstr += id_name + "\t"
            outstr += str(feature.qualifiers.get("gene", [unique_name])[0]) + "\t"
            outstr += str(feature.qualifiers.get("code", ["-"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("cog", ["-"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("product", ["-"])[0]) + "\t"
            all_features.append(outstr)
    sys.stdout.write("%s - 1..%d\n"%(outchrm, len(gb.seq)))
    sys.stdout.write("%d %s\n"%(len(all_features), name))
    sys.stdout.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
    for feature in all_features:
        sys.stdout.write(feature + "\n")


def out_tsv(gb, chrm = None, features = ["CDS"], qual_names = ["locus_tag"]):
    """
    Gives unique name, locus tag, gene name, and product for a given set
    """
    if chrm:
        outchrm = str(chrm)
    else:
        outchrm = str(gb.name)

    sys.stdout.write("chr\tstart\tend\tstrand\tfeature_type\tunique_name\tlocus_tag\tprotein_id\tgene\tproduct\n")
    tsv_name = set()
    for feature in gb.features:
        outstr = ""
        if feature.type in features:
            # skip features with partial ends
            if check_partial_start(feature.location.start) or check_partial_end(feature.location.end):
                continue
            outstr += outchrm + "\t"
            outstr += str(feature.location.start) + "\t"
            outstr += str(feature.location.end) + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\t"
            outstr += str(feature.type) + "\t"
            #figure out unique name
            unique_name = check_qualifiers(feature, qual_names)
            num = 1
            id_name = unique_name
            while id_name in tsv_name:
                id_name = unique_name + "_%d"%(num)
                num += 1
            tsv_name.add(id_name)
            outstr += id_name + "\t"
            outstr += str(feature.qualifiers.get("locus_tag", ["NA"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("protein_id", ["NA"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("gene", ["NA"])[0]) + "\t"
            outstr += str(feature.qualifiers.get("product", ["NA"])[0]) + "\n"
 #            # pull the nucleotide sequence
 #            nuc_seq = feature.extract(gb.seq)
 #            this_translation = nuc_seq.translate(table = "Bacterial", cds = True)
 #            translation = str(feature.qualifiers['translation'][0])
 #            # check that the translation matches the expected product. Have to drop the stop codon in
 #            # the nucleotide translation
 #            if not (this_translation == translation):
 #                raise ValueError("Nucleotide sequence doesn't match translation\n%s\n%s"%(this_translation, translation))
 #            outstr += str(nuc_seq) + "\t"
 #            outstr += translation + "\n"
            sys.stdout.write(outstr)


def out_bed(gb, features = ["CDS"], chrm = None, qual_names = ["locus_tag"]):
    if chrm:
        name = str(chrm)
    else:
        name = str(gb.name)

    bed_name = set()
    for feature in gb.features:
        if feature.type in features:
            # skip features with partial ends
            if check_partial_start(feature.location.start) or check_partial_end(feature.location.end):
                continue
            outstr = ""
            outstr += name + "\t"
            outstr += str(feature.location.start) + "\t"
            outstr += str(feature.location.end) + "\t"
            # often times things need a gene name, thus will replace with locus tag or
            # if that doesn't exist. NA_number
            unique_name = check_qualifiers(feature, qual_names)
            num = 1
            id_name = unique_name
            while id_name in bed_name:
                id_name = unique_name + "_%d"%(num)
                num += 1
            bed_name.add(id_name)
            outstr += id_name + "\t"
            outstr += "." + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\n"
            sys.stdout.write(outstr)

def convert_qualifiers_to_string(qualifiers_dict):
    outstr = ""
    for key in qualifiers_dict:
        outstr += "%s=%s;"%(key, ",".join(qualifiers_dict[key]))
    return outstr[:-1]


def out_gff(gb, features = ["CDS"], chrm = None):
    if chrm:
        name = str(chrm)
    else:
        name = str(gb.name)
    sys.stdout.write("#gff-version 3\n")
    for feature in gb.features:
        if feature.type in features:
            # skip features with partial ends
            if check_partial_start(feature.location.start) or check_partial_end(feature.location.end):
                continue
            outstr = ""
            outstr += name + "\t"
            outstr += "." + "\t"
            # we are not doing heiarchical relationships here so skip
            # type
            outstr += "." + "\t"
            # gffs are in 1-based format
            outstr += str(feature.location.start + 1) + "\t"
            outstr += str(feature.location.end) + "\t"
            outstr += "." + "\t"
            outstr += str(convert_strandval(feature.location.strand)) + "\t"
            outstr += "." + "\t"
            # get all qualifiers
            outstr += convert_qualifiers_to_string(feature.qualifiers) + "\n"
            sys.stdout.write(outstr)

def out_fasta(gb, features = ["CDS"], qual_names = ["locus_tag"]):
    import fasta as fa
    out_fasta = fa.FastaFile()
    fa_name = set()
    for feature in gb.features:
        if feature.type in features:
            # skip features with partial ends
            if check_partial_start(feature.location.start) or check_partial_end(feature.location.end):
                continue

            # often times things need a gene name, thus will replace with locus tag or
            # if that doesn't exist. NA_number
            unique_name = check_qualifiers(feature, qual_names)
            num = 1
            id_name = unique_name
            while id_name in fa_name:
                id_name = unique_name + "_%d"%(num)
                num += 1
            fa_name.add(id_name)
            header = ">" + id_name + \
            " " +  \
            str(feature.qualifiers.get("product", ["NA"])[0])
            seq = str(feature.extract(gb.seq))
            entry = fa.FastaEntry(header = header, seq = seq)
            out_fasta.add_entry(entry)
    out_fasta.write(sys.stdout)

def out_fna(gb, chrm = None):
    import fasta as fa
    out_fasta = fa.FastaFile()
    if chrm:
        header = ">" + str(chrm)
    else:
        header = ">" + str(gb.name)
    seq = str(gb.seq)
    entry = fa.FastaEntry(header = header, seq = seq)
    out_fasta.add_entry(entry)
    out_fasta.write(sys.stdout)

def check_qualifiers(entry, qual_names, unknown_name = "NA"):
    for qual in qual_names:
        if qual in entry.qualifiers:
            return str(entry.qualifiers.get(qual)[0])
    return unknown_name



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a genbank file for features")
    parser.add_argument("infile", type=str, help="path to input file")
    parser.add_argument("--outfmt", type=str, help='''output format. default = bed. 
                                                    Supported: [fa, fasta, bed, tsv, ptt, rnt, fna]
                                                    fa or fasta give each features nucleotide sequence as a fasta record.
                                                    ptt and rnt give format files for rockhopper input.
                                                    bed gives each feature in bed format.
                                                    fna pulls the full nucleotide sequence.
                                                    Note that tsv is a special format
                                                    that pulls several fields for a given query''',
                        default = "bed")
    parser.add_argument('--features', type=str, nargs="+", help='''
    feature types to parse for. Can specify multiple. default = "CDS"
    ''', default = ["CDS"])
    parser.add_argument('--chrm', type=str, help='''
    specify the chromosome name for the output
    ''', default = None)
    parser.add_argument('--qual_name', type=str, nargs = "+", help='''
    which qualifier fields to parse to get name of feature? If multiple are listed, will check in order until a field is filled in. Default = locus_tag
    ''', default = "locus_tag")
    args = parser.parse_args()

    fname = args.infile
    ftype = args.outfmt
    features = args.features
    chrm = args.chrm
    qual_name = args.qual_name
    with open(fname, mode = "r") as inf:
        gb = SeqIO.read(inf, "genbank")
    if ftype == "fasta" or ftype == "fa":
        out_fasta(gb, features, qual_names = qual_name)
    elif ftype == "bed":
        out_bed(gb, features, chrm = chrm, qual_names = qual_name)
    elif ftype == "ptt" or ftype == "rnt":
        out_ptt_rnt(gb, ftype, chrm = chrm, qual_names = qual_name)
    elif ftype == "gff":
        out_gff(gb, features, chrm = chrm)
    elif ftype == "fna":
        out_fna(gb, chrm = chrm)
    elif ftype == "tsv":
        out_tsv(gb, chrm = chrm, features = features, qual_names = qual_name)
