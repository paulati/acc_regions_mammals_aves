import wget
from os.path import join, isfile
from source.config.config import MafDownloadConfig
from source.common.aws_s3 import *


# def download_chr_files(chr, file_part_count):
def download_chr_files(chr):

    config = MafDownloadConfig()

    output_directory = config.local_output_directory()

    filename = config.file_name_head() + str(chr) + config.file_name_tail()
    print(filename)
    url_base = config.url_base()
    url = url_base + filename

    out_file_path = join(url_base, filename)

    if isfile(out_file_path):
        result = out_file_path
    else:
        result = wget.download(url, out=output_directory)
    return result


def upload_chr_files(chr):
    config = AwsS3Config()
    bucket_name = config.bucket_name()
    aws_s3 = AwsS3(bucket_name)

    filename = config.remote_file_name_head() + str(chr) + config.remote_file_name_tail()
    remote_file_key = join(config.remote_maf_base_path(), filename)

    maf_config = MafDownloadConfig()
    output_directory = maf_config .local_output_directory()

    local_file_key = join(output_directory, filename)

    aws_s3.upload_file(remote_file_key, local_file_key)


####################################################
#config = MafDownloadConfig()
#files_per_chr = config.files_count_per_chr()

# for chrom in range(33, 34):
#
#     download_chr_files(chrom)
#
#     # upload raw maf to s3
#     upload_chr_files(chrom)


