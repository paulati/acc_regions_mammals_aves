import boto3
import botocore
from boto3.s3.transfer import *
from source.config.config import AwsS3Config


class AwsS3:

    def __init__(self, bucket_name):
        aws_s3_config = AwsS3Config()
        self.s3_client = boto3.client('s3', aws_access_key_id=aws_s3_config.account_key(),
                                      aws_secret_access_key=aws_s3_config.account_secret())
        self.bucket_name = bucket_name

    def download_file(self, remote_file_key, local_file_key):
        try:
            self.s3_client.download_file(self.bucket_name, remote_file_key, local_file_key)
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] == "404":
                print("The object does not exist.")
            else:
                raise

    def upload_file(self, remote_file_key, local_file_key):

        config = TransferConfig(multipart_threshold=1024 * 25, max_concurrency=10, multipart_chunksize=1024 * 25,
                                use_threads=True)
        try:
            self.s3_client.upload_file(local_file_key, self.bucket_name, remote_file_key, Config=config)
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] == "404":
                print("bla bla")
            else:
                raise

    def file_exists(self, key):
        try:
            obj = self.s3_client.head_object(Bucket=self.bucket_name, Key=key)
            return obj['ContentLength']
        except ClientError as exc:
            if exc.response['Error']['Code'] != '404':
                raise
