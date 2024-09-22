from Bio import AlignIO

# from Bio.SeqIO import *

file_path = "/paula/2018/acc_regions/scoring/data/knownCanonical.exonAA.fa"

# for record in parse(file_path, "fasta"):
#    print(record.id)

# https://biopython.org/wiki/Multiple_Alignment_Format
# https://biopython.org/wiki/AlignIO

alignment = AlignIO.read(open(file_path), "fasta")
print("Alignment length %i" % alignment.get_alignment_length())

AlignIO.write(alignment, "/paula/2018/acc_regions/scoring/data/knownCanonical.exonAA.maf", "maf")


# https://biopython.org/wiki/AlignIO
# File format conversion

# from Bio import AlignIO

# input_handle = open("example.phy", "rU")
# output_handle = open("example.sth", "w")
#
# alignments = AlignIO.parse(input_handle, "phylip")
# AlignIO.write(alignments, output_handle, "stockholm")

# output_handle.close()
# input_handle.close()