from typing import *
import yaml

class dotdict(dict):
    """
    dot.notation access to dictionary attributes
    """
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getattr__(self, item):
        return self[item]

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
    """
    with open(path) as f:
        cfg = yaml.safe_load(f)
        return dotdict.create(cfg)