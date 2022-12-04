from typing import *

import os
import yaml
from yamlinclude import YamlIncludeConstructor


class dotdict(dict):
    """
    dot.notation access to dictionary attributes
    TODO: keep somehow reference to the original YAML in order to report better
    KeyError origin.
    """
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            return self.__getattribute__(item)

    @classmethod
    def create(cls, cfg : Any):
        """
        - recursively replace all dicts by the dotdict.
        """
        if isinstance(cfg, dict):
            items = ( (k, cls.create(v)) for k,v in cfg.items())
            return dotdict(items)
        elif isinstance(cfg, list):
            return [cls.create(i) for i in cfg]
        elif isinstance(cfg, tuple):
            return tuple([cls.create(i) for i in cfg])
        else:
            return cfg


def load_config(path):
    """
    Load configuration from given file replace, dictionaries by dotdict
    uses pyyaml-tags namely for:
    include tag:
        geometry: <% include(path="config_geometry.yaml")>
    """
    YamlIncludeConstructor.add_to_loader_class(loader_class=yaml.FullLoader, base_dir=os.path.dirname(path))
    with open(path) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
    return dotdict.create(cfg)