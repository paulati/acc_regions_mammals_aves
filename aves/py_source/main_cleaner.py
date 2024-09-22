from source.maf_cleaner.s3_maf_downloader import *
from source.maf_cleaner.maf_reader import *
from os.path import splitext, basename


def main(chrom, branch, deleteFiles):

    #bajar los alineamientos de s3 a local
    maf_downloader = S3MafDownloader()
    maf_downloader.download(chrom)
    maf_file_names = maf_downloader.chrom_remote_files_name(chrom)

    # read + clean maf file + translate strand - in ref specie + sort alignment respect refseq start
    if branch == 'aves':
        maf_reader_config = AvesMafReaderConfig()
    elif branch == 'sarcopterygii':
        maf_reader_config = SarcopterygiiMafReaderConfig()

    ref_specie = maf_reader_config.maf_ref_specie()
    maf_species = maf_reader_config.maf_alignment_species()

    maf_base_path_in = maf_reader_config.local_input_directory()
    maf_base_path_out = maf_reader_config.local_output_directory()

    for local_file_name in maf_file_names:
        file_name, file_extension = splitext(local_file_name)
        if file_extension != ".gz":
            file_name = file_name + "." + file_extension
        print(file_name)
        maf_reader = MafReader(maf_base_path_in, maf_base_path_out, file_name)
        alignments = maf_reader.read()
        alignments_clean = maf_reader.clean(alignments, ref_specie, maf_species)
        maf_reader.write(alignments_clean)


        #upload to s3
        maf_file_path_out = join(maf_base_path_out, file_name)
        maf_file_path_out += '.gz'
        if branch == 'aves':
            aws_s3_config = AvesAwsS3Config()
        elif branch == 'sarcopterygii':
            aws_s3_config = SarcopterygiiAwsS3Config()

        bucket_name = aws_s3_config.bucket_name()
        aws_s3 = AwsS3(bucket_name)
        remote_file_name = basename(maf_file_path_out)
        print(remote_file_name)
        remote_file_key = join(aws_s3_config.remote_chr_maf_base_path(), remote_file_name)
        aws_s3.upload_file(remote_file_key, maf_file_path_out)

        # delete local files
        # delete all files in this folders:
        if deleteFiles:
            maf_base_path_in = maf_reader_config.local_input_directory()
            file_list = [f for f in listdir(maf_base_path_in) if f.endswith(".maf") or f.endswith(".gz")]
            for f in file_list :
                remove(join(maf_base_path_in, f))
            maf_base_path_out = maf_reader_config.local_output_directory()
            file_list = [f for f in listdir(maf_base_path_out) if f.endswith(".maf") or f.endswith(".gz")]
            for f in file_list:
                remove(join(maf_base_path_out, f))


chrom = 30

branch = 'sarcopterygii'
deleteFiles = False
main(chrom, branch, deleteFiles)

branch = 'aves'
deleteFiles = True
main(chrom, branch, deleteFiles)

