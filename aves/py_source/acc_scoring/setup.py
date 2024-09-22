import os
from os.path import join

# create folders structure

# TODO: move to config
base_path = "/home/rstudio/disco_tmp"

data1_directory = join(base_path, "data")
if not os.path.exists(data1_directory):
    os.makedirs(data1_directory)

data_directory = join(data1_directory, "scoring")
if not os.path.exists(data_directory):
    os.makedirs(data_directory)


# ###################

maf_directory = join(data_directory, "maf")
if not os.path.exists(maf_directory):
    os.makedirs(maf_directory)

maf_index_directory = join(data_directory, "maf_indexes")
if not os.path.exists(maf_index_directory):
    os.makedirs(maf_index_directory)

acc_maf_directory = join(data_directory, "acc_maf")
if not os.path.exists(acc_maf_directory):
    os.makedirs(acc_maf_directory)

acc_maf_outgroup_directory = join(data_directory, "acc_maf_outgroup")
if not os.path.exists(acc_maf_outgroup_directory):
    os.makedirs(acc_maf_outgroup_directory)

acc_maf_ingroup_directory = join(data_directory, "acc_maf_ingroup")
if not os.path.exists(acc_maf_ingroup_directory):
    os.makedirs(acc_maf_ingroup_directory)

# ###################

acc_maf_25_directory = join(acc_maf_directory, "25")
if not os.path.exists(acc_maf_25_directory):
    os.makedirs(acc_maf_25_directory)

for i in range(1, 23):
    chr_directory = join(acc_maf_25_directory, "chr" + str(i))
    if not os.path.exists(chr_directory):
        os.makedirs(chr_directory)

acc_maf_50_directory = join(acc_maf_directory, "50")
if not os.path.exists(acc_maf_50_directory):
    os.makedirs(acc_maf_50_directory)

for i in range(1, 23):
    chr_directory = join(acc_maf_50_directory, "chr" + str(i))
    if not os.path.exists(chr_directory):
        os.makedirs(chr_directory)

# ###################

acc_maf_outgroup_25_directory = join(acc_maf_outgroup_directory, "25")
if not os.path.exists(acc_maf_outgroup_25_directory):
    os.makedirs(acc_maf_outgroup_25_directory)

for i in range(1, 23):
    chr_directory = join(acc_maf_outgroup_25_directory, "chr" + str(i))
    if not os.path.exists(chr_directory):
        os.makedirs(chr_directory)


acc_maf_outgroup_50_directory = join(acc_maf_outgroup_directory, "50")
if not os.path.exists(acc_maf_outgroup_50_directory):
    os.makedirs(acc_maf_outgroup_50_directory)

for i in range(1, 23):
    chr_directory = join(acc_maf_outgroup_50_directory, "chr" + str(i))
    if not os.path.exists(chr_directory):
        os.makedirs(chr_directory)

# ###################

acc_maf_ingroup_25_directory = join(acc_maf_ingroup_directory, "25")
if not os.path.exists(acc_maf_ingroup_25_directory):
    os.makedirs(acc_maf_ingroup_25_directory)

for i in range(1, 23):
    chr_directory = join(acc_maf_ingroup_25_directory, "chr" + str(i))
    if not os.path.exists(chr_directory):
        os.makedirs(chr_directory)


acc_maf_ingroup_50_directory = join(acc_maf_ingroup_directory, "50")
if not os.path.exists(acc_maf_ingroup_50_directory):
    os.makedirs(acc_maf_ingroup_50_directory)

for i in range(1, 23):
    chr_directory = join(acc_maf_ingroup_50_directory, "chr" + str(i))
    if not os.path.exists(chr_directory):
        os.makedirs(chr_directory)

# ###################

acc_maf_outgroup_consensus_directory = join(data_directory, "acc_maf_outgroup_consensus")
if not os.path.exists(acc_maf_outgroup_consensus_directory):
    os.makedirs(acc_maf_outgroup_consensus_directory)


acc_maf_ingroup_consensus_directory = join(data_directory, "acc_maf_ingroup_consensus")
if not os.path.exists(acc_maf_ingroup_consensus_directory):
    os.makedirs(acc_maf_ingroup_consensus_directory)
