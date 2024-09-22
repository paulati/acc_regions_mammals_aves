from os import listdir
from os.path import isfile, join
import numpy as np

from acc_elements_grammar import *


# interesa el signo de la distancia, no solo el valor absoluto
# def distance(p1, p2, p3):
# # distancia de p3 a la recta formada por p1 y p2
#
#     x = p3[0]
#     y = p3[1]
#     if x > y:
#         sign = -1
#     else:
#         sign = 1
#
#     denominator = np.linalg.norm(p2 - p1)
#     numerator = np.linalg.norm(np.cross(p2 - p1, p1 - p3))
#     result = sign * numerator / denominator
#     #d = norm(np.cross(p2 - p1, p1 - p3)) / norm(p2 - p1)
#     return result


def process_file(elements_base_path, elements_file_name):

    file = join(elements_base_path, elements_file_name)

    elem_score_parser = ElemScoreParser(file)

    elements = elem_score_parser.score_elements

    output_file_path = join(elements_base_path, "filter", "filter_" + elements_file_name)

    f = open(output_file_path, 'w')

    for element in elements:

        # hamming_ingroup = element.hamming_ingroup
        #
        # hamming_outgroup = element.hamming_outgroup
        #
        # p1 = np.array([0, 0], dtype=np.float32)
        # p2 = np.array([1, 1], dtype=np.float32)
        # p3 = np.array([hamming_ingroup, hamming_outgroup], dtype=np.float32)
        #
        # d = distance(p1, p2, p3)
        #
        parts_slash = str(element.elem_id).split("/")
        length = parts_slash[6]
        chromosome = parts_slash[7]
        maf_file = parts_slash[8]
        parts_point = maf_file.split(".")
        parts_underscore = parts_point[0].split("_")
        start = parts_underscore[3]
        end = parts_underscore[4]

        if (element.hamming_in_out_excl_gaps != 0) and (element.shift_count != 0):
            line = str(element.elem_id) + "\t" + \
                length + "\t" + \
                str(element.size) + "\t" + \
                str(element.gap_count) + "\t" + \
                str(element.shift_count) + "\t" + \
                str(element.hamming_ingroup_incl_gaps) + "\t" + \
                str(element.hamming_outgroup_incl_gaps) + "\t" + \
                str(element.hamming_in_out_incl_gaps) + "\t" + \
                str(element.hamming_ingroup_excl_gaps) + "\t" + \
                str(element.hamming_outgroup_excl_gaps) + "\t" + \
                str(element.hamming_in_out_excl_gaps) + "\t" + \
                str(chromosome) + ":" + str(start) + "-" + str(end) + "\t" + \
                str(chromosome) + "\t" + \
                str(start) + "\t" + \
                str(end) + "\n"

            # los ultimos tres campos los uso para armar un bed file
            f.write(line)

        print(element.elem_id)
        print(element.shift_count)
        print(element.hamming_in_out_excl_gaps)

        # print(hamming_ingroup)
        # print(hamming_outgroup)


def join_filtered_files(filtered_elements_base_path):

    only_files = [f for f in listdir(filtered_elements_base_path) if isfile(join(filtered_elements_base_path, f))]

    f_out_path = join(filtered_elements_base_path, "join", "join_filtered_elements.csv")
    f_out = open(f_out_path, "w")

    for file in only_files:
        file_path = join(filtered_elements_base_path, file)
        f_in = open(file_path, "r")
        lines_in = f_in.readlines()
        f_out.writelines(lines_in)


def main(elements_base_path):

    # only_files = [f for f in listdir(elements_base_path) if isfile(join(elements_base_path, f))]

    # for file in only_files:
    #     print(file)
    #     process_file(elements_base_path, file)

    filtered_elements_base_path = join(elements_base_path, "filter")
    join_filtered_files(filtered_elements_base_path)

    #file = "/paula/2018/acc_regions/scoring/data/results/chr13_results_25.txt"
    #file_name = "chr13_results_25.txt"
    #print(file)




elements_base_path = "/paula/2018/acc_regions/scoring/data/results/201901"
main(elements_base_path)


#
# filename = "/paula/2018/acc_regions/scoring/data/results/chr13_results_25.txt"
# elem_score_parser = ElemScoreParser(filename)
# result = elem_score_parser.score_elements
#
# print(result[0].elem_id)
# print(result[0].hamming_outgroup)
# print(result[0].shift_count)
# print(result[0].elem_score_average)
# print(result[0].hamming_ingroup)
# print(result[0].hamming_outgroup)
# print(result[0].hamming_relative)
#

# test distance:
#
# hamming_ingroup = 0
# hamming_outgroup = 0.5
# p1 = np.array([0, 0], dtype=np.float32)
# p2 = np.array([1, 1], dtype=np.float32)
# p3 = np.array([hamming_ingroup, hamming_outgroup], dtype=np.float32)
#
# test_1 = distance(p1, p2, p3)
# print(test_1)
#
# hamming_ingroup = 0
# hamming_outgroup = 0.8
# p1 = np.array([0, 0], dtype=np.float32)
# p2 = np.array([1, 1], dtype=np.float32)
# p3 = np.array([hamming_ingroup, hamming_outgroup], dtype=np.float32)
#
# test_2 = distance(p1, p2, p3)
# print(test_2)
#
# hamming_ingroup = 0.5
# hamming_outgroup = 0
# p1 = np.array([0, 0], dtype=np.float32)
# p2 = np.array([1, 1], dtype=np.float32)
# p3 = np.array([hamming_ingroup, hamming_outgroup], dtype=np.float32)
#
#
# test_3 = distance(p1, p2, p3)
# print(test_3)
#
#
# hamming_ingroup = 0.8
# hamming_outgroup = 0
# p1 = np.array([0, 0], dtype=np.float32)
# p2 = np.array([1, 1], dtype=np.float32)
# p3 = np.array([hamming_ingroup, hamming_outgroup], dtype=np.float32)
#
# test_4 = distance(p1, p2, p3)
# print(test_4)
#




















