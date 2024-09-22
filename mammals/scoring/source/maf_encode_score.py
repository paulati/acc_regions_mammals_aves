from Bio import AlignIO
from os import listdir
import numpy as np
from maf_consensus import MafConsensus
from config import MafConsensusConfig, MafEncodeScoreConfig
from os.path import isfile, join
from scipy.spatial import distance


class MafEncodeScore:

    def __init__(self, ingroup_maf_base_path, chromosome, shift_specie):
        only_files = [f for f in listdir(ingroup_maf_base_path) if isfile(join(ingroup_maf_base_path, f))]
        print("cantidad de archivos maf a procesar:")
        print(str(len(only_files)))
        self.maf_files = {}
        for maf_file_name in only_files:
            self.maf_files[maf_file_name] = join(ingroup_maf_base_path, maf_file_name)
        self.__load_consensus(chromosome)
        self.shift_specie = shift_specie

    def __load_consensus(self, chromosome):
        maf_consensus = MafConsensus(chromosome)

        maf_consensus_config = MafConsensusConfig()

        # load_consensus_dict(self, consensus_base_path, consensus_file_name_head, consensus_file_name_tail)
        outgroup_consensus_base_path = maf_consensus_config.outgroup_consensus_base_path()
        outgroup_consensus_file_name_head = maf_consensus_config.outgroup_consensus_file_name_head()
        outgroup_consensus_file_name_tail = maf_consensus_config.outgroup_consensus_file_name_tail()
        self.outgroup_consensus_dict = maf_consensus.load_consensus_dict(outgroup_consensus_base_path,
                                                                         outgroup_consensus_file_name_head,
                                                                         outgroup_consensus_file_name_tail)

        ingroup_consensus_base_path = maf_consensus_config.ingroup_consensus_base_path()
        ingroup_consensus_file_name_head = maf_consensus_config.ingroup_consensus_file_name_head()
        ingroup_consensus_file_name_tail = maf_consensus_config.ingroup_consensus_file_name_tail()
        self.ingroup_consensus_dict = maf_consensus.load_consensus_dict(ingroup_consensus_base_path,
                                                                        ingroup_consensus_file_name_head,
                                                                        ingroup_consensus_file_name_tail)

    def process(self, out_file_path):

        result = {}

        file = open(out_file_path, "w")

        print("cantidad de archivos a procesar:")
        print(len(self.maf_files.keys()))

        for maf_file_name in self.maf_files.keys():

            maf_file_path = self.maf_files[maf_file_name]

            maf_score_data = self.__process_maf_file(maf_file_name)

            # TODO: cambiar este return a struct o clase
            # (cols_sc, block_sc, shifts_c, hamming_in, hamming_out, hamming_ratio) = self.__process_maf_file(
            #     maf_file_name)
            #
            result[maf_file_path] = maf_score_data

        for k in result.keys():
            maf_score_data = result[k]

            score_columns = str(maf_score_data.columns_score).replace("\n", "\t")

            line = k + "\t" + \
                score_columns + "\t" + \
                str(maf_score_data.block_score) + "\t" + \
                str(maf_score_data.elem_size) + "\t" + \
                str(maf_score_data.gap_count) + "\t" + \
                str(maf_score_data.shift_count) + "\t" + \
                str(maf_score_data.hamming_ingroup_incl_gaps) + "\t" + \
                str(maf_score_data.hamming_outgroup_incl_gaps) + "\t" + \
                str(maf_score_data.hamming_in_out_incl_gaps) + "\t" + \
                str(maf_score_data.hamming_ingroup_excl_gaps) + "\t" + \
                str(maf_score_data.hamming_outgroup_excl_gaps) + "\t" + \
                str(maf_score_data.hamming_in_out_excl_gaps) + "\n"

            # print(line)
            file.write(line)

        file.close()

    # def __count_shifts(self, alignment, maf_file_name):

    @staticmethod
    def get_specie_sequence(specie_name, alignment):

        result = ""

        for sequence in alignment:

            sequence_specie = sequence.id.split(".")[0]
            # print(sequence_specie)

            if sequence_specie == specie_name:
                result = sequence.seq

        return result

    @staticmethod
    def __alignment_bit_matrix(alignment_matrix):
        result = np.zeros(alignment_matrix.shape, dtype=int, order='C')
        row_count = alignment_matrix.shape[0]
        col_count = alignment_matrix.shape[1]
        for col_index in range(col_count):
            # la ultima fila no la miro y reviso hasta la de indice 0
            for row_index in range(row_count-2, -1, -1):
                # esto NO debe ser case sensitive:
                if(str.lower(str(alignment_matrix[row_index][col_index])) ==
                        str.lower(str(alignment_matrix[row_index+1][col_index]))):
                    result[row_index][col_index] = 1
                else:
                    result[row_index][col_index] = 0
        # print(result)
        return result

    @staticmethod
    def __alignment_matrix(alignment):

        rows_count = len(alignment)
        if rows_count > 0:
            cols_count = len(alignment[0].seq)
        else:
            cols_count = 0

        result = np.chararray((rows_count, cols_count))
        result[:] = '.'

        row_index = 0
        for seq_rec in alignment:
            result[row_index] = seq_rec.seq
            row_index = row_index + 1

        return result

    def __process_maf_file(self, maf_file_name):

        maf_file_path = self.maf_files[maf_file_name]

        alignment = AlignIO.read(maf_file_path, "maf")  # my input alignment
        align_matrix = MafEncodeScore.__alignment_matrix(alignment)
        align_bit_mat = MafEncodeScore.__alignment_bit_matrix(align_matrix)

        if maf_file_name in self.outgroup_consensus_dict.keys():
            outgroup_consensus_seq = self.outgroup_consensus_dict[maf_file_name]
        else:
            outgroup_consensus_seq = ""

        if maf_file_name in self.ingroup_consensus_dict.keys():
            ingroup_consensus_seq = self.ingroup_consensus_dict[maf_file_name]
        else:
            ingroup_consensus_seq = ""

        shift_seq = MafEncodeScore.get_specie_sequence(self.shift_specie, alignment)

        result = MafEncodeScore.MafScoreData(maf_file_path, outgroup_consensus_seq, ingroup_consensus_seq,
                                             shift_seq, align_bit_mat)

        return result

    class MafScoreData:

        # def __init__(self, maf_encode_score):
        def __init__(self, id, outgroup_consensus_seq, ingroup_consensus_seq, shift_seq, align_bit_mat):

            self._align_bit_mat = align_bit_mat
            self._outgroup_consensus_seq = outgroup_consensus_seq
            self._ingroup_consensus_seq = ingroup_consensus_seq
            self._shift_seq = shift_seq

            self.id = id
            self.elem_size = len(self._shift_seq)
            self.columns_score = self.__score_columns()
            self.block_score = self.__score_block()

            # se completan cuando cuento gaps:
            self._outgroup_consensus_seq_gaps = None
            self._ingroup_consensus_seq_gaps = None

            # se completan cuando cuento shifts:
            self._outgroup_shifts = None
            self._ingroup_shifts = None

            if len(outgroup_consensus_seq) > 0 and len(ingroup_consensus_seq) > 0 and len(shift_seq) > 0:

                # self.maf_encode_score = maf_encode_score
                self.shift_count = self.__count_shifts()
                self.gap_count = self.__count_gaps()

                self.hamming_ingroup_incl_gaps = self.__hamming(self._ingroup_consensus_seq, self._shift_seq)
                self.hamming_outgroup_incl_gaps = self.__hamming(self._outgroup_consensus_seq, self._shift_seq)

                # excluyo solo los gaps que tambien son shifts
                # si excluyo todos los gaps en general se reduce mucho el tam del array
                # y la distancia de hamming resulta mayor que la que tengo incluyendo gaps (yo queria lo contratio)
                # (ver si es util calcular la medida excuyendo TODOS los gaps)
                ingroup_gap_and_shift = np.logical_and(self._ingroup_consensus_seq_gaps, self._ingroup_shifts)
                ingroup_gap_and_shift_index = ingroup_gap_and_shift.astype(int)
                shift_seq_excl_gaps_ingroup = MafEncodeScore.MafScoreData.__substring_by_indexes(
                    self._shift_seq, ingroup_gap_and_shift_index)

                ingroup_consensus_seq_excl_gaps = MafEncodeScore.MafScoreData.__substring_by_indexes(
                    self._ingroup_consensus_seq, ingroup_gap_and_shift_index)

                self.hamming_ingroup_excl_gaps = self.__hamming(ingroup_consensus_seq_excl_gaps, shift_seq_excl_gaps_ingroup)

                outgroup_gap_and_shift = np.logical_and(self._outgroup_consensus_seq_gaps, self._outgroup_shifts)
                outgroup_gap_and_shift_index = outgroup_gap_and_shift.astype(int)
                shift_seq_excl_gaps_outgroup = MafEncodeScore.MafScoreData.__substring_by_indexes(
                    self._shift_seq, outgroup_gap_and_shift_index)

                outgroup_consensus_seq_excl_gaps = MafEncodeScore.MafScoreData.__substring_by_indexes(
                    self._outgroup_consensus_seq, outgroup_gap_and_shift_index)

                self.hamming_outgroup_excl_gaps = self.__hamming(outgroup_consensus_seq_excl_gaps, shift_seq_excl_gaps_outgroup)

            else:
                self.shift_count = 0.0
                self.gap_count = 0.0
                self.hamming_ingroup_incl_gaps = 0.0
                self.hamming_outgroup_incl_gaps = 0.0
                self.hamming_ingroup_excl_gaps = 0.0
                self.hamming_outgroup_excl_gaps = 0.0

            # distancia a la recta (0,0),(1,1):
            self.hamming_in_out_incl_gaps = MafEncodeScore.MafScoreData.__distance_to_eq(
                self.hamming_ingroup_incl_gaps, self.hamming_outgroup_incl_gaps)
            self.hamming_in_out_excl_gaps = MafEncodeScore.MafScoreData.__distance_to_eq(
                self.hamming_ingroup_excl_gaps, self.hamming_outgroup_excl_gaps)

        @staticmethod
        def __substring_by_indexes(seq, gaps_indexes):
            gaps_count = sum(gaps_indexes)
            result_len = len(gaps_indexes) - gaps_count
            result = [None] * int(result_len)
            result_index = 0
            for i in range(0, len(gaps_indexes)):
                if gaps_indexes[i] == 0:
                    result[result_index] = seq[i]
                    result_index += 1
            return result

        # def __hamming(self, consensus_seq, include_gaps):
        def __hamming(self, consensus_seq, shift_seq):
            consensus_vec = list(consensus_seq)
            hamming = distance.hamming(consensus_vec, shift_seq)
            return hamming

        def __get_gaps(self, consensus_sequence):
            shift_sequence = self._shift_seq
            result = None
            if len(shift_sequence) == len(consensus_sequence):
                result = np.zeros(len(shift_sequence))
                # print(shift_sequence)
                # print(consensus_sequence)
                for i in range(len(shift_sequence)):
                    base_consensus_sequence = consensus_sequence[i]
                    if base_consensus_sequence == "-":
                        result[i] = 1
                    else:
                        result[i] = 0
            else:
                print("error, las secuencias deben ser de la misma longitud")
            return result

        def __count_gaps(self):

            # load instance variables needed in hamming without gaps
            self._outgroup_consensus_seq_gaps = self.__get_gaps(self._outgroup_consensus_seq)
            self._ingroup_consensus_seq_gaps = self.__get_gaps(self._ingroup_consensus_seq)

            outgroup_gaps = self._outgroup_consensus_seq_gaps

            ingroup_gaps = self._ingroup_consensus_seq_gaps

            result_gaps = np.zeros(len(outgroup_gaps))

            for i in range(len(outgroup_gaps)):
                elem_consensus_outgroup_gap = outgroup_gaps[i]
                elem_consensus_ingroup_gap = ingroup_gaps[i]
                result_gaps[i] = elem_consensus_outgroup_gap == 1 or elem_consensus_ingroup_gap == 1

            gaps_count = np.sum(result_gaps)

            return gaps_count

        def __get_shifts(self, consensus_sequence):
            shift_sequence = self._shift_seq
            result = None
            if len(shift_sequence) == len(consensus_sequence):
                result = np.zeros(len(shift_sequence))
                # print(shift_sequence)
                # print(consensus_sequence)
                for i in range(len(shift_sequence)):
                    base_shift_sequence = shift_sequence[i]
                    base_consensus_sequence = consensus_sequence[i]
                    if base_consensus_sequence == "X" or base_consensus_sequence == "-":
                        result[i] = 0
                    else:
                        if base_consensus_sequence == base_shift_sequence:
                            result[i] = 0
                        else:
                            result[i] = 1
            else:
                print("error, las secuencias deben ser de la misma longitud")
            return result

        def __count_shifts(self):

            self._outgroup_shifts = self.__get_shifts(self._outgroup_consensus_seq)
            outgroup_shifts = self._outgroup_shifts
            self._ingroup_shifts = self.__get_shifts(self._ingroup_consensus_seq)
            ingroup_shifts = self._ingroup_shifts

            result_shifts = np.zeros(len(outgroup_shifts))

            for i in range(len(outgroup_shifts)):
                elem_consensus_outgroup_shift = outgroup_shifts[i]
                elem_consensus_ingroup_shift = ingroup_shifts[i]
                result_shifts[i] = elem_consensus_outgroup_shift == 1 and elem_consensus_ingroup_shift == 0

            shifts_count = np.sum(result_shifts)

            return shifts_count

        def __score_columns(self):

            column_score = np.zeros(self._align_bit_mat.shape[1])

            for column_index in range(self._align_bit_mat.shape[1]):

                # print(align_bit_mat[:, column_index])

                # para cada columna, veo cuales son las filas que tienen un 0
                column = self._align_bit_mat[:, column_index]
                zero_indexes = np.where(column == 0)[0]
                # print(zero_indexes)

                # score column:

                col_score = 0
                from_index = 0
                for zero_index in zero_indexes:
                    to_index = zero_index + 1
                    # print("from: " + str(from_index))
                    # print("to: " + str(to_index))
                    block = column[from_index:to_index]  # incluye al limite inferior, NO incluye al limite superior
                    block_sum = np.sum(block)
                    block_len = len(block)

                    # print(block)
                    # print("len: " + str(block_len))
                    # print("sum: " + str(block_sum))
                    from_index = to_index
                    col_score = col_score + float(block_sum) / float(block_len)
                    # print("col_score: " + str(col_score))

                column_score[column_index] = col_score

            # print(column_score)

            return column_score

        def __score_block(self):
            score_block = np.sum(self.columns_score) / len(self.columns_score)
            return score_block

        # interesa el signo de la distancia, no solo el valor absoluto
        @staticmethod
        def __distance(p1, p2, p3):
        # distancia de p3 a la recta formada por p1 y p2

            x = p3[0]
            y = p3[1]
            if x > y:
                sign = -1
            else:
                sign = 1

            denominator = np.linalg.norm(p2 - p1)
            numerator = np.linalg.norm(np.cross(p2 - p1, p1 - p3))
            result = sign * numerator / denominator
            #d = norm(np.cross(p2 - p1, p1 - p3)) / norm(p2 - p1)
            return result

        @staticmethod
        def __distance_to_eq(coord_x, coord_y):
            p1 = np.array([0, 0], dtype=np.float32)
            p2 = np.array([1, 1], dtype=np.float32)
            #p3 = np.array([hamming_ingroup, hamming_outgroup], dtype=np.float32)
            p3 = np.array([coord_x, coord_y], dtype=np.float32)
            result = MafEncodeScore.MafScoreData.__distance(p1, p2, p3)
            return result


def main():

    chromosome = 22

    maf_encode_score_config = MafEncodeScoreConfig()

    ingroup_maf_base_path = maf_encode_score_config.ingroup_maf_base_path() + str(chromosome)

    shift_specie = maf_encode_score_config.shift_specie()

    maf_encode_score = MafEncodeScore(ingroup_maf_base_path, chromosome, shift_specie)

    out_file_path = maf_encode_score_config.out_file_path_head() + str(chromosome) + maf_encode_score_config.\
        out_file_path_tail()

    maf_encode_score.process(out_file_path)


# main()
