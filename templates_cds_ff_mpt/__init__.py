# -*- coding: utf-8 -*-

import os
import pkg_resources

import yaml

_yaml_file = pkg_resources.resource_filename(__name__, os.path.join('data', 'tech_params.yaml'))

with open(_yaml_file, 'r') as f:
    config = yaml.load(f)
