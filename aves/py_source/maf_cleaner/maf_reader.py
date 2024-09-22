import gzip
import shutil
from os import remove, listdir
from os.path import join, exists, splitext, basename, isfile
from Bio import AlignIO
from Bio.AlignIO import MultipleSeqAlignment


class MafReader:

    def __init__(self, maf_base_path_in, maf_base_path_out, maf_file_name=None):
        self.maf_base_path_in = maf_base_path_in
        self.maf_base_path_out = maf_base_path_out
        self.maf_file_name = maf_file_name

    # def __create_maf_index(self):
    #     dummy = 0
    #

    # def __clean_maf_species(self, species):
    #     to_remove = 'ancestral_sequences.Ancestor_'
    #     result = []
    #     for specie in species:
    #         if not specie.startswith(to_remove):
    #             dot_index = specie.index(".")
    #             specie_clean_name = specie[:dot_index]
    #             result.append(specie_clean_name)
    #     return result

    def __clean_maf_specie_name(self, specie):
        dot_index = specie.index(".")
        result = specie[:dot_index]
        return result

    # def list_maf_species(self):
    #     local_file_path = join(self.maf_base_path_in, self.maf_file_name)
    #     file_name_gz = self.maf_file_name + ".gz"
    #     gz_file_path = join(self.maf_base_path_in, file_name_gz)
    #     if not exists(local_file_path):
    #         f_in = gzip.open(gz_file_path, 'r')
    #         f_out = open(local_file_path, 'wb')
    #         shutil.copyfileobj(f_in, f_out)
    #     alignment_count = 0
    #     species = []
    #     for multiple_alignment in AlignIO.parse(local_file_path, "maf"):
    #         print("printing a new multiple alignment")
    #
    #         alignment_count = alignment_count + 1
    #         seqrec_count = 0
    #
    #         for seqrec in multiple_alignment:
    #             seqrec_count = seqrec_count + 1
    #             specie_name = seqrec.id
    #             species.append(specie_name)
    #
    #             # print("starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
    #             #      (seqrec.annotations["start"],
    #             #    seqrec.annotations["strand"],
    #             #    seqrec.annotations["srcSize"],
    #             #    seqrec.annotations["size"]))
    #
    #         print(seqrec_count)
    #     print(alignment_count)
    #
    #     species_clean = self.__clean_maf_species(species)
    #
    #     # delete duplicates
    #     result = list(set(species_clean))
    #     print(result)
    #     return result

    @staticmethod
    def complement_seq(seq_rec):
        # start en el strand contrario = srcSize - (start strand inicial + size)
        start = seq_rec.annotations["start"]
        src_size = seq_rec.annotations["srcSize"]
        size = seq_rec.annotations["size"]

        start_comp = src_size - (start + size)
        seq_comp = seq_rec.seq.reverse_complement()
        strand_comp = - seq_rec.annotations["strand"]
        result = seq_rec
        result.annotations["start"] = start_comp
        result.annotations["strand"] = strand_comp
        result.seq = seq_comp

        return result

    def read(self, maf_base_path=None, maf_file_name=None):
        if maf_base_path is None:
            maf_base_path = self.maf_base_path_in

        if maf_file_name is None:
            maf_file_name = self.maf_file_name

        local_file_path = join(maf_base_path, maf_file_name)

        if not exists(local_file_path):
            file_name, file_extension = splitext(local_file_path)
            if file_extension == ".gz":
                gz_file_path = local_file_path
            else:
                gz_file_path = file_name + file_extension + ".gz"
            f_in = gzip.open(gz_file_path, 'r')
            f_out = open(local_file_path, 'wb')
            shutil.copyfileobj(f_in, f_out)
        result = AlignIO.parse(local_file_path, "maf")
        return result

    def write(self, alignments, maf_base_path=None, maf_file_name=None):
        if maf_base_path is None:
            maf_base_path = self.maf_base_path_out
        if maf_file_name is None:
            maf_file_name = self.maf_file_name
        local_file_path = join(maf_base_path, maf_file_name)
        output_handle = open(local_file_path, "w")
        AlignIO.write(alignments, output_handle, "maf")
        output_handle.close()
        #compress maf:
        gz_file_path = local_file_path + '.gz'
        f_in = open(local_file_path, 'rb')
        f_out = gzip.open(gz_file_path, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()

    def sort(self, alignments):
        ref_seq_starts = {} # dict with key = start, value = alignment
        result = []
        for alignment in alignments:
            ref_seq_start = alignment[0].annotations["start"]
            ref_seq_starts[ref_seq_start] = alignment

        sort_starts = sorted(ref_seq_starts.keys())
        for sort_start in sort_starts:
            result.append(ref_seq_starts[sort_start])

        return result

    def clean(self, alignments, ref_specie, alignment_species):
        multiple_alignments_clean = []
        for multiple_alignment in alignments:
            seq_records = []
            _is_ref_seq_rec = True
            for seq_rec in multiple_alignment:
                specie_name = seq_rec.id
                specie_name_clean = self.__clean_maf_specie_name(specie_name)
                if specie_name_clean == ref_specie and _is_ref_seq_rec:
                    _is_ref_seq_rec = False
                    if seq_rec.annotations["strand"] == -1:
                        ref_seq_rec = MafReader.complement_seq(seq_rec)
                    else:
                        ref_seq_rec = seq_rec
                    seq_records.append(ref_seq_rec)
                elif specie_name_clean in alignment_species:
                    seq_records.append(seq_rec)

            align = MultipleSeqAlignment(seq_records)
            multiple_alignments_clean.append(align)

            #    species.append(specie_name)

                # print("starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
                #      (seqrec.annotations["start"],
                #    seqrec.annotations["strand"],
                #    seqrec.annotations["srcSize"],
                #    seqrec.annotations["size"]))

        # species_clean = self.__clean_maf_species(species)

        result = self.sort(multiple_alignments_clean)
        return result


# class MafJoiner:
#     def __init__(self):
#         dummy = 0
#
#     def out_file_name(self, in_file_name):
#         chrom = -1
#         file_name_parts = in_file_name.split(".")
#         if len(file_name_parts) > 2:
#             chrom_part = file_name_parts[2]
#             parts = chrom_part.split("_")
#             if len(parts) > 0:
#                 chrom = parts[0]
#         maf_joiner_config = MafJoinerConfig()
#         result = maf_joiner_config.out_file_name_head() + str(chrom) + maf_joiner_config.out_file_name_tail()
#         return(result)
#
#     # mover este metodo a otra clase porque usa un conjutno de alineamientos no solo uno
#     def join(self, maf_files):
#
#         result = None
#
#         if len(maf_files) > 0:
#             maf_joiner_config = MafJoinerConfig()
#
#             maf_base_path_in = maf_joiner_config.local_input_directory()
#             maf_base_path_out = maf_joiner_config.local_output_directory()
#             maf_reader = MafReader(maf_base_path_in, maf_base_path_out)
#
#             alignments_all = []
#             for maf_file in maf_files:
#                 file_name, file_extension = splitext(maf_file)
#                 if file_extension != ".gz":
#                     file_name = file_name + "." + file_extension
#                 alignments = maf_reader.read(None, file_name)
#                 for multiple_alignment in alignments:
#                     alignments_all.append(multiple_alignment)
#             result = maf_reader.sort(alignments_all)
#
#             maf_file_name = self.out_file_name(maf_files[0])
#             maf_reader.write(result, maf_base_path_out, maf_file_name)
#
#             # zip file:
#             maf_file_path_in = join(maf_base_path_out, maf_file_name)
#             maf_file_path_out = maf_file_path_in + ".gz"
#             f_in = open(maf_file_path_in, 'rb')
#             f_out = gzip.open(maf_file_path_out, 'wb')
#             shutil.copyfileobj(f_in, f_out)
#             f_out.close()
#             f_in.close()
#             remove(maf_file_path_in)
#
#         return maf_file_path_out


# class MafReaderChr:
#
#     def __init__(self, chrom, maf_base_path, maf_file_name_head, maf_file_name_tail):
#         local_file_name = maf_file_name_head + str(chrom) + maf_file_name_tail
#         local_file_path = join(maf_base_path, local_file_name)
#         self.align = AlignIO.parse(local_file_path, "maf")
#
#     def alignment(self):
#         return self.align

