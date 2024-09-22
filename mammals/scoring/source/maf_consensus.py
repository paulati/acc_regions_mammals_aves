from Bio.Align import AlignInfo
from os import listdir
from os.path import isfile, join
from maf_in_out_group import MafInOutGroup
from config import MafConsensusConfig


class MafConsensus:

    def __init__(self, chromosome):
        self.chromosome = str(chromosome)

    def load_consensus_dict(self, consensus_base_path, consensus_file_name_head, consensus_file_name_tail):

        # outgroup_consensus_base_path: "/paula/2018/acc_regions/scoring/data/acc_maf_outgroup_consensus"
        # outgroup_consensus_file_name_head: "chr"
        # outgroup_consensus_file_name_tail: "_consensus.txt"

        # consensus_base_path = "/paula/2018/acc_regions/scoring/data/acc_maf_outgroup_consensus"
        # consensus_file_name = "chr" + self.chromosome + "_consensus.txt"
        # out_maf_base_path = join(consensus_base_path, consensus_file_name)

        consensus_file_name = consensus_file_name_head + self.chromosome + consensus_file_name_tail
        consensus_path = join(consensus_base_path, consensus_file_name)
        # print(consensus_path)

        consensus = {}
        f = open(consensus_path)
        for line in f:
            tmp = line.split()
            if len(tmp) == 2:
                (key, val) = tmp
                consensus[key] = val

        return consensus

    def calculate_consensus(self, maf_base_path_head, out_base_path, consensus_file_name_head, consensus_file_name_tail,
                            threshold, species_consensus_file_path):

        maf_base_path = maf_base_path_head + self.chromosome

        only_files = [f for f in listdir(maf_base_path) if isfile(join(maf_base_path, f))]

        file_lines = []

        consensus_species = MafInOutGroup.load_species(species_consensus_file_path)

        # print(consensus_species)

        for maf_file in only_files:

            # maf_file_path = join(maf_base_path, maf_file)
            # alignment = AlignIO.read(maf_file_path, "maf")

            alignment = MafInOutGroup.filter_maf(maf_base_path, maf_file, consensus_species)

            align_info = AlignInfo.SummaryInfo(alignment)

            consensus = align_info.gap_consensus(threshold=threshold, ambiguous='X', consensus_alpha=None,
                                                 require_multiple=0)

            line = maf_file + "\t" + str(consensus)

            file_lines.append(line)

        consensus_file_name = consensus_file_name_head + self.chromosome + consensus_file_name_tail
        consensus_file_path = join(out_base_path,  consensus_file_name)

        # print(consensus_file_path)

        # print(consensus_file_path)
        out_file = open(consensus_file_path, 'w')
        for line in file_lines:
            out_file.write(line + '\n')
        out_file.close()


def main(chromosome):

    # outgroup consensus
    configuration = MafConsensusConfig()
    maf_base_path = configuration.outgroup_maf_base_path()
    out_base_path = configuration.outgroup_consensus_base_path()
    threshold = configuration.outgroup_consensus_threshold()
    outgroup_consensus_file_name_head = configuration.outgroup_consensus_file_name_head()
    outgroup_consensus_file_name_tail = configuration.outgroup_consensus_file_name_tail()
    outgroup_species_consensus_file_path = configuration.outgroup_sorted_species_consensus_file_path()

    maf_consensus = MafConsensus(chromosome)
    maf_consensus.calculate_consensus(maf_base_path, out_base_path, outgroup_consensus_file_name_head,
                                      outgroup_consensus_file_name_tail, threshold,
                                      outgroup_species_consensus_file_path)

    # ingroup consensus
    configuration = MafConsensusConfig()
    maf_base_path = configuration.ingroup_maf_base_path()
    out_base_path = configuration.ingroup_consensus_base_path()
    threshold = configuration.ingroup_consensus_threshold()
    ingroup_consensus_file_name_head = configuration.ingroup_consensus_file_name_head()
    ingroup_consensus_file_name_tail = configuration.ingroup_consensus_file_name_tail()
    ingroup_species_consensus_file_path = configuration.ingroup_sorted_species_consensus_file_path()

    maf_consensus = MafConsensus(chromosome)
    maf_consensus.calculate_consensus(maf_base_path, out_base_path, ingroup_consensus_file_name_head,
                                      ingroup_consensus_file_name_tail, threshold, ingroup_species_consensus_file_path)


# chrom = "22"
# main(chrom)
