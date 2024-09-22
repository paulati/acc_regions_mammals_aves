
# remote_base_folder_path = "download/alignments/multiz100way"
        # file_name_gz = "chr" + str(chromosome) + ".maf.gz"
        # remote_file_key = join(remote_base_folder_path, file_name_gz)
        # print(remote_file_key)
        #
        # local_base_folder_path = "/paula"
        # local_file_name_gz = "chr" + str(chromosome) + ".maf.gz"
        # local_file_key = join(local_base_folder_path, local_file_name_gz)

        # print(bucket.name)




# account.key <- ""
# account.secret <- ""
#
# chromosome = 22
# bucket_name = "acc-regions-2018"
#
#
# remote_base_folder_path = "download/alignments/multiz100way"
#
# #bucket_prefix = remote_base_folder_path
# #objs = bucket.objects.filter(Prefix=bucket_prefix)
#
# # for obj in objs:
# #     print(obj)
#
#     # path, filename = split(obj.key)
#     # # boto3 s3 download_file will throw exception if folder not exists
#     # try:
#     #     #os.makedirs(path)
#     #     print(path)
#     #     print(filename)
#     # except FileExistsError:
#     #     pass
#     #mybucket.download_file(obj.key, obj.key)
#
#
# # buckets = s3.buckets.all()
# # for bucket in buckets:
# #    print (bucket.name)
#

#

# bucket_prefix="/some/prefix/here"
# objs = mybucket.objects.filter(
#     Prefix = bucket_prefix)
#
# for obj in objs:
#     path, filename = os.path.split(obj.key)
#     # boto3 s3 download_file will throw exception if folder not exists
#     try:
#         os.makedirs(path)
#     except FileExistsError:
#         pass
#     mybucket.download_file(obj.key, obj.key)



#
#
#
#
# remote_base_folder_path = "/download/alignments/multiz100way/"
#
# local_base_folder_path = "/home/rstudio/2018/phastCons/data/tmp"
#
# file_name_gz = "chr11.maf.gz"
#
# local_file_name_gz = "chr11.maf.gz"
#
# # setwd(local.base.folder.path)
#
# # object.name < - paste(remote.base.folder.path, file.name.gz, sep="")
#
# access_key = account_key
# secret_key = account_secret
#
# conn = boto.connect_s3(
#         aws_access_key_id=access_key,
#         aws_secret_access_key=secret_key,
#         host='s3.console.aws.amazon.com/s3/home?region=us-east-1#',
#         is_secure=False,               # uncomment if you are not using ssl
#         calling_format=boto.s3.connection.OrdinaryCallingFormat(),
#         )
#
# for bucket in conn.get_all_buckets():
#     text = "{name}\t{created}".format(name = bucket.name, created = bucket.creation_date)
#     print(text)
#
# dummy = 1

#
# setwd(local.base.folder.path)
#
# aws.s3::save_object(object=object.name,
#                     key=account.key,
#                     secret=account.secret,
#                     bucket=bucket.name,
#                     file=local.file.name.gz)
#
# gc()
#
# raw.object < - aws.s3::get_object(object=object.name,
#                                   key=account.key,
#                                   secret=account.secret,
#                                   bucket=bucket.name)

#
#
# def upload_file():
#
#
#     chromosome = 22
#
#     s3_client = boto3.client('s3')
#
#     bucket_name = "acc-regions-2018"
#
#     remote_base_folder_path = "preparation/indexes_multiz100way"
#     file_name_gz = "chr" + str(chromosome) + ".maf.gz"
#     remote_file_key = join(remote_base_folder_path, file_name_gz)
#
#     print(remote_file_key)
#
#     local_base_folder_path = "/paula"
#     local_file_name_gz = "chr" + str(chromosome) + ".maf.gz"
#     local_file_key = join(local_base_folder_path, local_file_name_gz)
#
#     # FILE_PATH = '/path/to/file/'
#     # KEY_PATH = "/path/to/s3key/"
#
#     config = TransferConfig(multipart_threshold=1024 * 25, max_concurrency=10, multipart_chunksize=1024 * 25,
#                             use_threads=True)
#     # file = FILE_PATH + filename
#     file = local_file_key
#
#     # key = KEY_PATH + filename
#     key = remote_file_key
#
#     s3_client.upload_file(file, bucket_name, key,
#                           # ExtraArgs={'ACL': 'public-read', 'ContentType': 'video/mp4'},
#                           Config=config
#                           # Callback=ProgressPercentage(file)
#     )
#
#
#
# upload_file()



# test:
#
# s3_client = boto3.client('s3', aws_access_key_id='',
#                                       aws_secret_access_key='')
# chromosome = 22
#
# bucket_name = "acc-regions-2018"
#
# remote_base_folder_path = "preparation/indexes_multiz100way"
# file_name_gz = "chr" + str(chromosome) + ".maf.gz"
# remote_file_key = join(remote_base_folder_path, file_name_gz)
#
# print(remote_file_key)
#
# local_base_folder_path = "/paula"
# local_file_name_gz = "chr" + str(chromosome) + ".maf.gz"
# local_file_key = join(local_base_folder_path, local_file_name_gz)
#
# # FILE_PATH = '/path/to/file/'
# # KEY_PATH = "/path/to/s3key/"
#
# config = TransferConfig(multipart_threshold=1024 * 25, max_concurrency=10, multipart_chunksize=1024 * 25,
#                         use_threads=True)
# # file = FILE_PATH + filename
# file = local_file_key
#
# # key = KEY_PATH + filename
# key = remote_file_key
#
# s3_client.upload_file(file, bucket_name, key,
#                       # ExtraArgs={'ACL': 'public-read', 'ContentType': 'video/mp4'},
#                       Config=config
#                       # Callback=ProgressPercentage(file)
# )

