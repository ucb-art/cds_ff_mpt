# -*- coding: utf-8 -*-

import os
import pkg_resources

from bag.io import read_yaml

_yaml_file = pkg_resources.resource_filename(__name__, os.path.join('data', 'tech_params.yaml'))

config = read_yaml(_yaml_file)
