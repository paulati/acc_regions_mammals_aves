from os.path import join, exists
from aws_s3 import AwsS3
from config import BedExtractMafConfig, AwsS3Config
import shutil
import gzip
from Bio.AlignIO import MafIO


class MafAwsS3(AwsS3):

    def __init__(self):
        self.aws_s3_config = AwsS3Config()
        bucket_name = self.aws_s3_config.bucket_name()
        AwsS3.__init__(self, bucket_name)
        self.bed_extract_maf_config = BedExtractMafConfig()

    def prepare(self, chromosome):
        # download maf file
        self.__download_maf_file(chromosome)

        # unzip maf file
        self.__ungz_maf(chromosome)

        # if index does not exists in aws, create and upload index
        self.__create_maf_index(chromosome)

        # if index does not exists local  download index file
        self.__download_maf_index(chromosome)

    def __download_maf_file(self, chromosome):

        maf_base_path = self.bed_extract_maf_config.maf_base_path()
        maf_file_name_head = self.bed_extract_maf_config.maf_file_name_head()
        maf_file_name_tail = self.bed_extract_maf_config.maf_file_name_tail()
        local_file_name_gz = maf_file_name_head + str(chromosome) + maf_file_name_tail + ".gz"
        local_file_key = join(maf_base_path, local_file_name_gz)

        if not exists(local_file_key):
            remote_base_folder_path = self.aws_s3_config.remote_maf_base_path()
            file_name_gz = self.aws_s3_config.remote_maf_file_name_head() + str(chromosome) + \
                           self.aws_s3_config.remote_maf_file_name_tail()
            remote_file_key = join(remote_base_folder_path, file_name_gz)
            self.download_file(remote_file_key, local_file_key)

    def __ungz_maf(self, chromosome):

        maf_base_path = self.bed_extract_maf_config.maf_base_path()
        maf_file_name_head = self.bed_extract_maf_config.maf_file_name_head()
        maf_file_name_tail = self.bed_extract_maf_config.maf_file_name_tail()
        local_file_name = maf_file_name_head + str(chromosome) + maf_file_name_tail
        local_file_path = join(maf_base_path, local_file_name)

        file_name_gz = self.aws_s3_config.remote_maf_file_name_head() + str(chromosome) + self.aws_s3_config.\
            remote_maf_file_name_tail()
        gz_file_path = join(maf_base_path, file_name_gz)

        if not exists(local_file_path):
            f_in = gzip.open(gz_file_path, 'r')
            f_out = open(local_file_path, 'wb')
            shutil.copyfileobj(f_in, f_out)

    def __gz_maf_index(self, chromosome):
        multiple_alignment_file_name = self.bed_extract_maf_config.maf_file_name_head() + str(chromosome) + \
                                       self.bed_extract_maf_config.maf_file_name_tail()
        multiple_alignment_file_index_name = self.bed_extract_maf_config.index_maf_file_name_head() + \
                                             multiple_alignment_file_name
        multiple_alignment_file_index_path = join(self.bed_extract_maf_config.index_maf_base_path(),
                                                  multiple_alignment_file_index_name)
        multiple_alignment_file_index_path_gz = multiple_alignment_file_index_path + ".gz"
        f_in = open(multiple_alignment_file_index_path, 'rb')
        f_out = gzip.open(multiple_alignment_file_index_path_gz, 'wb')
        shutil.copyfileobj(f_in, f_out)

    def __ungz_maf_index(self, chromosome):

        multiple_alignment_file_name = self.bed_extract_maf_config.maf_file_name_head() + str(chromosome) + \
                                       self.bed_extract_maf_config.maf_file_name_tail()
        multiple_alignment_file_index_name = self.bed_extract_maf_config.index_maf_file_name_head() + \
                                             multiple_alignment_file_name
        multiple_alignment_file_index_path = join(self.bed_extract_maf_config.index_maf_base_path(),
                                                  multiple_alignment_file_index_name)
        multiple_alignment_file_index_path_gz = multiple_alignment_file_index_path + ".gz"

        if not exists(multiple_alignment_file_index_path):
            f_in = gzip.open(multiple_alignment_file_index_path_gz, 'r')
            f_out = open(multiple_alignment_file_index_path, 'wb')
            shutil.copyfileobj(f_in, f_out)



    def __create_maf_index(self, chromosome):
        remote_base_folder_path = self.aws_s3_config.remote_index_maf_base_path()
        file_name_gz = self.bed_extract_maf_config.index_maf_file_name_head() + \
                       self.aws_s3_config.remote_maf_file_name_head() + str(chromosome) + \
                       self.aws_s3_config.remote_maf_file_name_tail()
        remote_file_key = join(remote_base_folder_path, file_name_gz)
        file_size = self.file_exists(remote_file_key)
        if (file_size is None) or (file_size == 0):

            multiple_alignment_file_name = self.bed_extract_maf_config.maf_file_name_head() + str(chromosome) + \
                                           self.bed_extract_maf_config.maf_file_name_tail()
            multiple_alignment_file_index_name = self.bed_extract_maf_config.index_maf_file_name_head() + \
                                                 multiple_alignment_file_name
            multiple_alignment_file_index_path = join(self.bed_extract_maf_config.index_maf_base_path(),
                                                      multiple_alignment_file_index_name)

            multiple_alignment_file_index_path_gz = multiple_alignment_file_index_path + ".gz"

            multiple_alignment_file_path = join(self.bed_extract_maf_config.maf_base_path(),
                                                multiple_alignment_file_name)

            index_target_seqname = self.bed_extract_maf_config.index_target_seqname() + str(chromosome)

            index = MafIO.MafIndex(multiple_alignment_file_index_path, multiple_alignment_file_path, index_target_seqname)

            multiple_alignment_file_index_key_remote = join(self.aws_s3_config.remote_index_maf_base_path(),
                                                            multiple_alignment_file_index_name)

            multiple_alignment_file_index_key_remote_gz = multiple_alignment_file_index_key_remote + ".gz"

            self.__gz_maf_index(chromosome)

            self.upload_file(multiple_alignment_file_index_key_remote_gz, multiple_alignment_file_index_path_gz)

    def __download_maf_index(self, chromosome):

        multiple_alignment_file_name = self.bed_extract_maf_config.maf_file_name_head() + str(chromosome) + \
                                       self.bed_extract_maf_config.maf_file_name_tail()

        multiple_alignment_file_index_name = self.bed_extract_maf_config.index_maf_file_name_head() + \
                                             multiple_alignment_file_name

        multiple_alignment_file_index_name_gz = multiple_alignment_file_index_name + ".gz"

        multiple_alignment_file_index_key_remote = join(self.aws_s3_config.remote_index_maf_base_path(),
                                                        multiple_alignment_file_index_name_gz)

        multiple_alignment_file_index_path = join(self.bed_extract_maf_config.index_maf_base_path(),
                                                  multiple_alignment_file_index_name)

        multiple_alignment_file_index_path_gz = join(self.bed_extract_maf_config.index_maf_base_path(),
                                                     multiple_alignment_file_index_name_gz)

        if not exists(multiple_alignment_file_index_path):
            if not exists(multiple_alignment_file_index_path_gz):
                self.download_file(multiple_alignment_file_index_key_remote, multiple_alignment_file_index_path_gz)
            self.__ungz_maf_index(chromosome)


def main():
    chromosome = 22
    maf_aws_s3 = MafAwsS3()
    maf_aws_s3.prepare(chromosome)


# main()