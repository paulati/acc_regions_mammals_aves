import yaml


# TODO: ver si tengo que hacer esta clase abstracta
class Config:

    def __init__(self, section):
        stream = open('./config/config.yaml', 'r')    # 'document.yaml' contains a single YAML document.
        self.config = yaml.load(stream)
        self.__env_config__()
        self.__env_config_section__(section)

    def __env_config__(self):
        tmp_config = self.config
        environment = tmp_config['environment']
        self.env_config = tmp_config[environment]

    def __env_config_section__(self, section):
        self.env_config_section = self.env_config[section]

    def file_name_head(self):
        return self.env_config_section['file_name_head']

    def file_name_tail(self):
        return self.env_config_section['file_name_tail']


class BedExtractMafConfig(Config):

    def __init__(self):
        section = 'bed_extract_maf'
        Config.__init__(self, section)

    def maf_base_path(self):
        return self.env_config_section['maf_base_path']

    def maf_file_name_head(self):
        return self.env_config_section['maf_file_name_head']

    def maf_file_name_tail(self):
        return self.env_config_section['maf_file_name_tail']

    def bed_base_path(self):
        return self.env_config_section['bed_base_path']

    def bed_file_name_tail(self):
        return self.env_config_section['bed_file_name_tail']

    def index_maf_file_name_head(self):
        return self.env_config_section['index_maf_file_name_head']

    def index_maf_base_path(self):
        return self.env_config_section['index_maf_base_path']

    def index_target_seqname(self):
        return self.env_config_section['index_target_seqname']

    def out_maf_base_path(self):
        return self.env_config_section['out_maf_base_path']


class MafInOutGroupConfig(Config):

    def __init__(self):
        section = 'maf_in_out_group'
        Config.__init__(self, section)

    def outgroup_sorted_species_file_path(self):
        return self.env_config_section['outgroup_sorted_species_file_path']

    def out_outgroup_base_path_head(self):
        return self.env_config_section['out_outgroup_base_path_head']

    def ingroup_sorted_species_file_path(self):
        return self.env_config_section['ingroup_sorted_species_file_path']

    def out_ingroup_base_path_head(self):
        return self.env_config_section['out_ingroup_base_path_head']

    def maf_base_path(self):
        return self.env_config_section['maf_base_path']

# TODO: ver si esto debe ser singleton


class MafConsensusConfig(Config):

    def __init__(self):
        section = 'maf_consensus'
        Config.__init__(self, section)

    def outgroup_consensus_base_path(self):
        return self.env_config_section['outgroup_consensus_base_path']

    def outgroup_consensus_file_name_head(self):
        return self.env_config_section['outgroup_consensus_file_name_head']

    def outgroup_consensus_file_name_tail(self):
        return self.env_config_section['outgroup_consensus_file_name_tail']

    def outgroup_consensus_threshold(self):
        return self.env_config_section['outgroup_consensus_threshold']

    def outgroup_maf_base_path(self):
        return self.env_config_section['outgroup_maf_base_path']

    def outgroup_sorted_species_consensus_file_path(self):
        return self.env_config_section['outgroup_sorted_species_consensus_file_path']

    def ingroup_consensus_base_path(self):
        return self.env_config_section['ingroup_consensus_base_path']

    def ingroup_consensus_file_name_head(self):
        return self.env_config_section['ingroup_consensus_file_name_head']

    def ingroup_consensus_file_name_tail(self):
        return self.env_config_section['ingroup_consensus_file_name_tail']

    def ingroup_consensus_threshold(self):
        return self.env_config_section['ingroup_consensus_threshold']

    def ingroup_maf_base_path(self):
        return self.env_config_section['ingroup_maf_base_path']

    def ingroup_sorted_species_consensus_file_path(self):
        return self.env_config_section['ingroup_sorted_species_consensus_file_path']


class MafEncodeScoreConfig(Config):

    def __init__(self):
        section = 'maf_score'
        Config.__init__(self, section)

    def ingroup_maf_base_path(self):
        return self.env_config_section['ingroup_maf_base_path']

    def shift_specie(self):
        return self.env_config_section['shift_specie']

    def out_file_path_head(self):
        return self.env_config_section['out_file_path_head']

    def out_file_path_tail(self):
        return self.env_config_section['out_file_path_tail']


class AwsS3Config(Config):
    def __init__(self):
        section = 'aws_s3'
        Config.__init__(self, section)

    def account_key(self):
        return self.env_config_section['key']

    def account_secret(self):
        return self.env_config_section['secret']

    def bucket_name(self):
        return self.env_config_section['bucket_name']

    def remote_maf_base_path(self):
        return self.env_config_section['remote_maf_base_path']

    def remote_maf_file_name_head(self):
        return self.env_config_section['remote_maf_file_name_head']

    def remote_maf_file_name_tail(self):
        return self.env_config_section['remote_maf_file_name_tail']

    def remote_index_maf_base_path(self):
        return self.env_config_section['remote_index_maf_base_path']

# config = Config('bed_extract_maf')
# print(config.bed_base_path())
