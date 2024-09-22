from maf_aws_s3 import MafAwsS3
from bed_extract_maf import BedExtractMAF
from maf_in_out_group import MafInOutGroup
from maf_consensus import MafConsensus
from maf_encode_score import MafEncodeScore
from config import BedExtractMafConfig, MafInOutGroupConfig, MafConsensusConfig, MafEncodeScoreConfig
from os.path import join


def process_chromosome(chromosome):

    maf_aws_s3 = MafAwsS3()
    maf_aws_s3.prepare(chromosome)

    # bed_extract_maf:

    configuration = BedExtractMafConfig()

    file_name_head = configuration.file_name_head() + str(chromosome)
    elem_file_name_base = file_name_head + configuration.file_name_tail()

    acc_elem_bed_file_name = elem_file_name_base + configuration.bed_file_name_tail()
    acc_bed_file_path = join(configuration.bed_base_path(), acc_elem_bed_file_name)

    multiple_alignment_file_name = configuration.maf_file_name_head() + str(chromosome) + configuration.\
        maf_file_name_tail()
    multiple_alignment_file_index_name = configuration.index_maf_file_name_head() + multiple_alignment_file_name
    multiple_alignment_file_index_path = join(configuration.index_maf_base_path(), multiple_alignment_file_index_name)
    multiple_alignment_file_path = join(configuration.maf_base_path(), multiple_alignment_file_name)

    index_target_seqname = configuration.index_target_seqname() + str(chromosome)

    acc_elem_maf_file_path_base = join(configuration.out_maf_base_path(), file_name_head, elem_file_name_base)

    bed_extract_maf = BedExtractMAF(acc_bed_file_path, multiple_alignment_file_index_path, multiple_alignment_file_path,
                                    index_target_seqname)
    bed_extract_maf.process(acc_elem_maf_file_path_base)

    # maf_in_out_group:

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

    # maf_consensus:

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
                                      outgroup_consensus_file_name_tail, threshold, outgroup_species_consensus_file_path
                                      )
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

    # maf_encode_score:

    maf_encode_score_config = MafEncodeScoreConfig()
    ingroup_maf_base_path = maf_encode_score_config.ingroup_maf_base_path() + str(chromosome)
    shift_specie = maf_encode_score_config.shift_specie()
    maf_encode_score = MafEncodeScore(ingroup_maf_base_path, chromosome, shift_specie)
    out_file_path = maf_encode_score_config.out_file_path_head() + str(chromosome) + maf_encode_score_config.\
        out_file_path_tail()
    maf_encode_score.process(out_file_path)


for chrom in range(1, 10):
    process_chromosome(chrom)

