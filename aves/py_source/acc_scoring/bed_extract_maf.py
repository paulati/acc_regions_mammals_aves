import gzip
import shutil
from pathlib import Path
from Bio import AlignIO
from Bio.AlignIO import MafIO
import pandas as pd
from config import *
from os.path import join


class BedExtractMAF:

    def __init__(self, bed_file_path, multiple_alignment_file_index_path,
                 multiple_alignment_file_path, index_target_seqname):

        self.__bed_elements = BedExtractMAF.__read_acc_bed(bed_file_path)

        # load maf index, create it if necessary:
        self.__maf_alignment_index = MafIO.MafIndex(multiple_alignment_file_index_path, multiple_alignment_file_path,
                                                    index_target_seqname)

    @staticmethod
    def __ungz(gz_file_path, extract_to_dir, extract_to_file_name):
        extract_to_file_path = join(extract_to_dir, extract_to_file_name)
        f_in = gzip.open(gz_file_path, 'r')
        f_out = open(extract_to_file_path, 'wb')
        shutil.copyfileobj(f_in, f_out)

    @staticmethod
    def __uncompress_maf(data_base_path, multiple_alignment_file_name):
        multiple_alignment_file_path = join(data_base_path, multiple_alignment_file_name)
        multiple_alignment_file_path_gz = multiple_alignment_file_path + ".gz"

        maf_file = Path(multiple_alignment_file_path)
        if not maf_file.exists():
            BedExtractMAF.__ungz(multiple_alignment_file_path_gz, data_base_path, multiple_alignment_file_name)

    @staticmethod
    def __read_acc_bed(bed_file_path):
        result = pd.read_csv(bed_file_path, sep='\t', header=None)
        result.columns = ['chr', 'start', 'end']
        return result

    def __extract_acc_elements_maf(self, acc_elements_start, acc_elements_end, acc_elem_maf_file_path_base):

        # tengo que escribir un archivo por cada elemento porque no hay forma de cortarlos de long 25-50
        # con los metodos disponibles
        for i in range(len(acc_elements_start)):

            # TODO: ver si es (start - 1) o (end + 1)
            elem_start = acc_elements_start[i] - 1
            elem_end = acc_elements_end[i]

            acc_elements_alignment = self.__maf_alignment_index.get_spliced([elem_start], [elem_end], strand=1)

            acc_maf_file_path = acc_elem_maf_file_path_base + "_" + str(
                elem_start) + "_" + str(elem_end) + ".maf"
            print(acc_maf_file_path)
            AlignIO.write(acc_elements_alignment, acc_maf_file_path, "maf")

    def process(self, acc_elem_maf_file_path_base):

        acc_elements_start = self.__bed_elements.start.values
        acc_elements_end = self.__bed_elements.end.values

        self.__extract_acc_elements_maf(acc_elements_start, acc_elements_end, acc_elem_maf_file_path_base)


def main():

    configuration = BedExtractMafConfig()

    chromosome = "22"

    file_name_head = configuration.file_name_head() + chromosome
    elem_file_name_base = file_name_head + configuration.file_name_tail()

    acc_elem_bed_file_name = elem_file_name_base + configuration.bed_file_name_tail()
    acc_bed_file_path = join(configuration.bed_base_path(), acc_elem_bed_file_name)

    multiple_alignment_file_name = configuration.maf_file_name_head() + chromosome + configuration.maf_file_name_tail()
    multiple_alignment_file_index_name = configuration.index_maf_file_name_head() + multiple_alignment_file_name
    multiple_alignment_file_index_path = join(configuration.index_maf_base_path(), multiple_alignment_file_index_name)
    multiple_alignment_file_path = join(configuration.maf_base_path(), multiple_alignment_file_name)

    index_target_seqname = configuration.index_target_seqname() + chromosome

    acc_elem_maf_file_path_base = join(configuration.out_maf_base_path(), file_name_head, elem_file_name_base)

    bed_extract_maf = BedExtractMAF(acc_bed_file_path, multiple_alignment_file_index_path, multiple_alignment_file_path,
                                    index_target_seqname)
    bed_extract_maf.process(acc_elem_maf_file_path_base)


# main()
