import yaml


# TODO: ver si tengo que hacer esta clase abstracta
class Config:

    def __init__(self, section, config_file_name=None):
        if config_file_name:
            stream = open('./config/' + config_file_name, 'r')
        else:
            stream = open('./config/config_multiz77way.yaml', 'r')    # 'document.yaml' contains a single YAML document.
        self.config = yaml.load(stream, Loader=yaml.FullLoader)
        self.__env_config__()
        self.__env_config_section__(section)

    def __env_config__(self):
        tmp_config = self.config
        environment = tmp_config['environment']
        self.env_config = tmp_config[environment]
        self.env_config_base = tmp_config['common']

    def __env_config_section__(self, section):
        if section in self.env_config:
            self.env_config_section = self.env_config[section]
        else:
            self.env_config_section = None
        if section in self.env_config_base:
            self.env_config_section_base = self.env_config_base[section]
        else:
            self.env_config_section_base = None

    def ___get_value__(self, key):

        if self.env_config_section is not None and self.env_config_section_base is not None:
            if key in self.env_config_section:
                result = self.env_config_section[key]
            elif key in self.env_config_section_base:
                result = self.env_config_section_base[key]
            else:
                result = None
        elif self.env_config_section is not None:
            if key in self.env_config_section:
                result = self.env_config_section[key]
            else:
                result = None
        elif self.env_config_section_base is not None:
            if key in self.env_config_section_base:
                result = self.env_config_section_base[key]
            else:
                result = None
        else:
            result = None

        return result

    # def file_name_head(self):
    #     return self.env_config_section['file_name_head']
    #
    # def file_name_tail(self):
    #     return self.env_config_section['file_name_tail']


class MafReaderConfig(Config):

    def __init__(self, config_file_name=None):
        section = 'maf_reader'
        Config.__init__(self, section, config_file_name)

    def maf_ref_specie(self):
        result = self.___get_value__('maf_ref_specie')
        return result

    def maf_alignment_species(self):
        result = self.___get_value__('maf_aln_species')
        return result

    def local_input_directory(self):
        result = self.___get_value__('local_input_directory')
        return result

    def local_output_directory(self):
        result = self.___get_value__('local_output_directory')
        return result

    def maf_file_name_head(self):
        result = self.___get_value__('maf_file_name_head')
        return result

    def maf_file_name_tail(self):
        result = self.___get_value__('maf_file_name_tail')
        return result

    # def maf_base_path(self):
    #     return self.env_config_section['maf_base_path']
    #
    # def maf_file_name_head(self):
    #     return self.env_config_section['maf_file_name_head']
    #
    # def maf_file_name_tail(self):
    #     return self.env_config_section['maf_file_name_tail']
    #
    # def bed_base_path(self):
    #     return self.env_config_section['bed_base_path']
    #
    # def bed_file_name_tail(self):
    #     return self.env_config_section['bed_file_name_tail']
    #
    # def index_maf_file_name_head(self):
    #     return self.env_config_section['index_maf_file_name_head']
    #
    # def index_maf_base_path(self):
    #     return self.env_config_section['index_maf_base_path']
    #
    # def index_target_seqname(self):
    #     return self.env_config_section['index_target_seqname']
    #
    # def out_maf_base_path(self):
    #     return self.env_config_section['out_maf_base_path']


class SarcopterygiiMafReaderConfig(MafReaderConfig):
    def __init__(self):
        section = 'maf_reader'
        config_file_name = 'config_multiz77way_sarcopterygii.yaml'
        Config.__init__(self, section, config_file_name)


class AvesMafReaderConfig(MafReaderConfig):
    def __init__(self):
        section = 'maf_reader'
        config_file_name = 'config_multiz77way_aves.yaml'
        Config.__init__(self, section, config_file_name)


class MafDownloadConfig(Config):

    def __init__(self):
        section = 'maf_downloader'
        Config.__init__(self, section)

    def url_base(self):
        result = self.___get_value__('url_base')
        return result

    def local_output_directory(self):
        result = self.___get_value__('local_output_directory')
        return result

    def file_name_head(self):
        result = self.___get_value__('file_name_head')
        return result

    def file_name_tail(self):
        result = self.___get_value__('file_name_tail')
        return result

    # def files_count_per_chr(self):
    #     result = self.___get_value__('files_per_chr')
    #     return result


class AwsS3Config(Config):

    def __init__(self, config_file_name=None):
        section = 'aws_s3'
        Config.__init__(self, section, config_file_name)

    def account_key(self):
        result = self.___get_value__('key')
        return result

    def account_secret(self):
        result = self.___get_value__('secret')
        return result

    def bucket_name(self):
        result = self.___get_value__('bucket_name')
        return result

    def remote_maf_base_path(self):
        result = self.___get_value__('remote_maf_base_path')
        return result

    def remote_file_name_head(self):
        result = self.___get_value__('remote_file_name_head')
        return result

    def remote_file_name_tail(self):
        result = self.___get_value__('remote_file_name_tail')
        return result

    def remote_chr_maf_base_path(self):
        result = self.___get_value__('remote_chr_maf_base_path')
        return result


class AvesAwsS3Config(AwsS3Config):

    def __init__(self):
        config_file_name = 'config_multiz77way_aves.yaml'
        AwsS3Config.__init__(self, config_file_name)


class SarcopterygiiAwsS3Config(AwsS3Config):

    def __init__(self):
        config_file_name = 'config_multiz77way_sarcopterygii.yaml'
        AwsS3Config.__init__(self, config_file_name)


class MafJoinerConfig(Config):
    def __init__(self):
        section = 'maf_joiner'
        Config.__init__(self, section)

    def local_input_directory(self):
        result = self.___get_value__('local_input_directory')
        return result

    def local_output_directory(self):
        result = self.___get_value__('local_output_directory')
        return result

    def out_file_name_head(self):
        result = self.___get_value__('out_file_name_head')
        return result

    def out_file_name_tail(self):
        result = self.___get_value__('out_file_name_tail')
        return result
