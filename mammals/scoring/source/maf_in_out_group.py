from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from os import listdir
from os.path import isfile, join
from config import *


class MafInOutGroup:

    @staticmethod
    def __clean_specie(data):
        result = []
        for line in data:
            tmp_1 = str.strip(line)
            if len(tmp_1) > 0:
                parts = str.split(tmp_1)
                for part in parts:
                    tmp_2 = str.strip(part)
                    if len(tmp_2) > 0:
                        result.append(tmp_2)
        return result

    @staticmethod
    def __process_maf(alignment, species):

        sorted_array = [None] * len(species)  # sets up an empty array that good sequences can be added to

        for sequence in alignment:

            sequence_specie = sequence.id.split(".")[0]
            print(sequence_specie)

            if sequence_specie in species:
                # do nothing
                print(sequence.id + " agregar")
                index = species.index(sequence_specie)
                print(index)
                # TODO: ver que estoy poniendo todas las letras de la secuencia en mayusculas
                sequence.seq = sequence.seq.upper()
                print(sequence.seq)
                sorted_array[index] = sequence
            # else:
                # do nothing

        # veo que ninguno de los elementos sea None y si hay alguno lo quito,
        # el paso previo me asegura que los tengo ordenados

        alignment_sorted_array = []
        for element in sorted_array:
            if element is not None:
                alignment_sorted_array.append(element)

        result = MultipleSeqAlignment(alignment_sorted_array)

        return result

    @staticmethod
    def load_species(sorted_species_file_path):
        file = open(sorted_species_file_path, "r")
        file_text = file.readlines()
        species = MafInOutGroup.__clean_specie(file_text)
        file.close()
        return species

    @staticmethod
    def filter_maf(maf_base_path, maf_file_name, species):

        maf_file_path = join(maf_base_path, maf_file_name)

        alignment = AlignIO.read(maf_file_path, "maf")  # my input alignment

        result = MafInOutGroup.__process_maf(alignment, species)

        return result

    @staticmethod
    def process_chromosome(chromosome, sorted_species_file_path, in_maf_base_path_head, out_maf_base_path_head):

        # read species from config file:
        # file = open(sorted_species_file_path, "r")
        # file_text = file.readlines()
        # species = MafInOutGroup.__clean_specie(file_text)
        # file.close()
        # species = ['ornAna1', 'allMis1', 'cheMyd1', 'chrPic2', 'pelSin1', 'apaSpi1', 'anoCar2', 'xenTro7', 'latCha1']

        species = MafInOutGroup.load_species(sorted_species_file_path)

        out_maf_base_path = out_maf_base_path_head + str(chromosome)

        in_maf_base_path = in_maf_base_path_head + str(chromosome)

        only_files = [f for f in listdir(in_maf_base_path) if isfile(join(in_maf_base_path, f))]

        for maf_file in only_files:

            result = MafInOutGroup.filter_maf(in_maf_base_path, maf_file, species)

            # maf_file_path = join(in_maf_base_path, maf_file)
            # alignment = AlignIO.read(maf_file_path, "maf")   # my input alignment
            # result = MafInOutGroup.__process_maf(alignment, species)

            out_maf_file_path = join(out_maf_base_path, maf_file)

            AlignIO.write(result, out_maf_file_path, "maf")  # writes a maf alignment containing only passing seqs


def main(chromosome):

    configuration = MafInOutGroupConfig()

    sorted_species_file_path = configuration.outgroup_sorted_species_file_path()
    in_maf_base_path_head = configuration.maf_base_path()
    out_maf_base_path_head = configuration.out_outgroup_base_path_head()
    MafInOutGroup.process_chromosome(chromosome,  sorted_species_file_path, in_maf_base_path_head,
                                     out_maf_base_path_head)

    sorted_species_file_path = configuration.ingroup_sorted_species_file_path()
    in_maf_base_path_head = configuration.maf_base_path()
    out_maf_base_path_head = configuration.out_ingroup_base_path_head()
    MafInOutGroup.process_chromosome(chromosome,  sorted_species_file_path, in_maf_base_path_head,
                                     out_maf_base_path_head)


# chrom = "22"
# main(chrom)
