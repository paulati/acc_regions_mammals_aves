from source.config.config import *
from source.common.aws_s3 import *
from os.path import join, isfile

class S3MafDownloader:

    def __init__(self):
        self.aws_s3_config = AwsS3Config()
        self.maf_downloader_config = MafDownloadConfig()

    def chrom_remote_files_name(self, chrom):
        remote_file_name = self.aws_s3_config.remote_file_name_head() + str(chrom) + \
                        self.aws_s3_config.remote_file_name_tail()
        result = [remote_file_name]
        return result

    def download(self, chrom):
        remote_file_names = self.chrom_remote_files_name(chrom)
        bucket_name = self.aws_s3_config.bucket_name()
        aws_s3 = AwsS3(bucket_name)
        for remote_file_name in remote_file_names:
            remote_file_key = join(self.aws_s3_config.remote_maf_base_path(), remote_file_name)
            local_file_key = join(self.maf_downloader_config.local_output_directory(), remote_file_name)

            if isfile(local_file_key):
                pass
            else:
                aws_s3.download_file(remote_file_key, local_file_key)

