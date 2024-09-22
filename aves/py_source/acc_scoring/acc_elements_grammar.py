from pyparsing import *


class ElemScoreGrammar:

    def __init__(self):

        acc_element_file_path_head = Literal("/home/rstudio/disco_tmp/data/acc_maf_ingroup/")

        self.acc_element_file_path = Combine(acc_element_file_path_head + OneOrMore(Word(printables))
                                             + WordEnd(".maf"))

        decimal = Combine(Word(nums) + ZeroOrMore("." + Word(nums))) ^ Combine(Word(nums) + Literal("."))

        self.acc_element_bp_score = Combine(Literal("[") + OneOrMore(decimal + ZeroOrMore(White())) + Literal("]"))

        acc_element_score_num = decimal
        acc_element_score_max = Literal("MAX")
        acc_element_score = acc_element_score_num ^ acc_element_score_max
        acc_element_score_sgn = Combine(Optional("-") + acc_element_score)

        # acc_element_score_bases_matrix_row = Combine(OneOrMore(acc_element_score + matrix_separator))
        # acc_element_score_bases_matrix = Combine("[" + OneOrMore(acc_element_score_bases_matrix_row) + "]")

        self.acc_element_score_average = acc_element_score

        self.acc_element_size = Word(nums)

        self.acc_element_gap_count = decimal

        self.acc_element_shift_count = decimal

        self.acc_element_hamming_ingroup_incl_gaps = acc_element_score

        self.acc_element_hamming_outgroup_incl_gaps = acc_element_score

        self.acc_element_hamming_in_out_incl_gaps = acc_element_score_sgn

        self.acc_element_hamming_ingroup_excl_gaps = acc_element_score

        self.acc_element_hamming_outgroup_excl_gaps = acc_element_score

        self.acc_element_hamming_in_out_excl_gaps = acc_element_score_sgn

        self.elem_data = self.acc_element_file_path.setResultsName("file_path") + \
                        self.acc_element_bp_score.setResultsName("bp_score") + \
                        self.acc_element_score_average.setResultsName("score_average") + \
                        self.acc_element_size.setResultsName("size") + \
                        self.acc_element_gap_count.setResultsName("gap_count") + \
                        self.acc_element_shift_count.setResultsName("shift_count") + \
                        self.acc_element_hamming_ingroup_incl_gaps.setResultsName("hamming_ingroup_incl_gaps") + \
                        self.acc_element_hamming_outgroup_incl_gaps.setResultsName("hamming_outgroup_incl_gaps")  + \
                        self.acc_element_hamming_in_out_incl_gaps.setResultsName("hamming_in_out_incl_gaps") + \
                        self.acc_element_hamming_ingroup_excl_gaps.setResultsName("hamming_ingroup_excl_gaps") + \
                        self.acc_element_hamming_outgroup_excl_gaps.setResultsName("hamming_outgroup_excl_gaps") + \
                        self.acc_element_hamming_in_out_excl_gaps.setResultsName("hamming_in_out_excl_gaps")

        # self.elem_data = self.acc_element_file_path + self.acc_element_bp_score

        # self.elem_last_line_tail = OneOrMore(self.acc_element_score) + Literal("]") + \
        #     self.acc_element_score_average.setResultsName("element_score_average") + \
        #     self.acc_element_shift_count.setResultsName("shift_count") + \
        #     self.acc_element_hamming_ingroup.setResultsName("hamming_ingroup") + \
        #     self.acc_element_hamming_outgroup.setResultsName("hamming_outgroup") + \
        #     self.acc_element_hamming_relative.setResultsName("hamming_relative")

        # self.elem_data = self.acc_element_file_path
                         # + \
        # self.acc_element_bp_score.setResultsName("element_bp_score") + \
        # self.acc_element_score_average.setResultsName("element_score_average") + \
        #               self.acc_element_shift_count.setResultsName("shift_count") + \
        #               self.acc_element_hamming_ingroup.setResultsName("hamming_ingroup") + \
        #               self.acc_element_hamming_outgroup.setResultsName("hamming_outgroup") + \
        #               self.acc_element_hamming_relative.setResultsName("hamming_relative")


class ElemScore:

    # TODO: checkear los tipos de datos de los atributos

    #def __init__(self, elem_id, elem_last_line):
    def __init__(self, line):

        grammar = ElemScoreGrammar()
        data = grammar.elem_data.parseString(line)
        self.elem_id = data["file_path"]
        self.bp_score = data["bp_score"]
        self.score_average = float(data["score_average"])
        self.size = float(data["size"])
        self.gap_count = float(data["gap_count"])
        self.shift_count = float(data["shift_count"])
        self.hamming_ingroup_incl_gaps = float(data["hamming_ingroup_incl_gaps"])
        self.hamming_outgroup_incl_gaps = float(data["hamming_outgroup_incl_gaps"])
        self.hamming_in_out_incl_gaps = float(data["hamming_in_out_incl_gaps"])
        self.hamming_ingroup_excl_gaps = float(data["hamming_ingroup_excl_gaps"])
        self.hamming_outgroup_excl_gaps = float(data["hamming_outgroup_excl_gaps"])
        self.hamming_in_out_excl_gaps = float(data["hamming_in_out_excl_gaps"])


class ElemScoreParser:

    def __init__(self, file_path):

        f = open(file_path, 'r')

        lines = f.readlines()

        self.score_elements = []

        for line in lines:
            print(line)
            elem_score = ElemScore(line)
            if elem_score.shift_count > 0:
                self.score_elements.append(elem_score)

        f.close()




# filename = "/paula/2018/acc_regions/scoring/data/results/chr13_results_25.txt"
# filename = "/paula/2018/acc_regions/scoring/data/results/201901/chr22_results_25.txt"
# elem_score_parser = ElemScoreParser(filename)
# result = elem_score_parser.score_elements
# print(result[0])

# line = "/home/rstudio/disco_tmp/data/acc_maf_ingroup/25/chr13/chr13_acc_25_105359296_105359321.maf	[0.94117647 0.94117647 0.94117647 0.94117647 2.3        1.42857143	 2.77380952 1.42307692 0.94117647 0.94117647 2.19047619 0.94117647	 0.94117647 0.94117647 0.94117647 0.94117647 0.94117647 1.5952381	 0.94117647 0.94117647 0.94117647 0.94117647 0.94117647 0.94117647	 0.94117647]	1.1837410040939453	25	0.0	0.0	0.08	0.04	-0.028284271099323854	0.08	0.04	-0.028284271099323854"
# grammar = ElemScoreGrammar()
# data = grammar.elem_data.parseString(line)
# dummy = 1