from dataclasses import dataclass
from typing import *

import os
import yaml
import re
from glob import iglob

from yamlinclude import YamlIncludeConstructor
from yamlinclude.constructor import WILDCARDS_REGEX, get_reader_class_by_name



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

Key = Union[str, int]
Path = Tuple[Key]

@dataclass
class PathIter:
    path: Path
    # full address path
    i: int = 0
    # actual level of the path; initial -1 is before first call to `idx` or `key`.

    def is_leaf(self):
        return self.i == len(self.path)

    def idx(self):
        try:
            return int(self.path[self.i]), PathIter(self.path, self.i + 1)
        except ValueError:
            raise IndexError(f"Variant substitution: IndexError at address: '{self.address()}'.")

    def key(self):
        key = self.path[self.i]
        if len(key) > 0 and not key[0].isdigit():
            return key, PathIter(self.path, self.i + 1)
        else:
            raise KeyError(f"Variant substitution: KeyError at address: '{self.address()}'.")

    def address(self):
        sub_path = self.path[:self.i + 1]
        return '/'.join([str(v) for v in sub_path])


def _item_update(key:Key, val:dotdict, sub_path:Key, sub:dotdict):
    sub_key, path = sub_path
    if key == sub_key:
        if path.empty():
            # Recursion termination
            return sub
        else:
            return deep_update(val, path, sub)
    else:
        return val

def deep_update(cfg: dotdict, iter:PathIter, substitute:dotdict):
    if iter.is_leaf():
        return substitute
    new_cfg = cfg.copy()
    if isinstance(cfg, list):
        key, sub_path = iter.idx()
    elif isinstance(cfg, (dict, dotdict)):
        key, sub_path = iter.key()
    else:
        raise TypeError(f"Variant substitution: Unknown type {type(cfg)}")
    new_cfg[key] = deep_update(cfg[key], sub_path, substitute)
    return new_cfg



def apply_variant(cfg:dotdict, variant:Dict[str, dotdict]) -> dotdict:
    """
    In the `variant` dict the keys are interpreted as the address
    in the YAML file. The address is a list of strings and ints separated by '/'
    and representing an item of the YAML file.
    For every `(address, value)` item of the `variant` dict the referenced item
    in `cfg` is replaced by `value`.

    Implemented by recursion with copy of changed collections.
    May be slow for too many variant items and substitution of the large collection.
    :param cfg:
    :param variant: dictionary path -> dotdict
    :return:
    """
    new_cfg = cfg
    for path_str, val in variant.items():
        path = path_str.split('/')
        assert path
        new_cfg = deep_update(new_cfg, PathIter(path), val)
    return new_cfg

class YamlInclude(YamlIncludeConstructor):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.included_files = []

    def load(
            self,
            loader,
            pathname: str,
            recursive: bool = False,
            encoding: str = '',
            reader: str = ''
    ):  # pylint:disable=too-many-arguments
        if not encoding:
            encoding = self._encoding or self.DEFAULT_ENCODING
        if self._base_dir:
            pathname = os.path.join(self._base_dir, pathname)
        reader_clz = None
        if reader:
            reader_clz = get_reader_class_by_name(reader)
        if re.match(WILDCARDS_REGEX, pathname):
            result = []
            iterable = iglob(pathname, recursive=recursive)
            for path in filter(os.path.isfile, iterable):
                self.included_files.append(path)
                if reader_clz:
                    result.append(reader_clz(path, encoding=encoding, loader_class=type(loader))())
                else:
                    result.append(self._read_file(path, loader, encoding))
            return result
        self.included_files.append(pathname)
        if reader_clz:
            return reader_clz(pathname, encoding=encoding, loader_class=type(loader))()
        return self._read_file(pathname, loader, encoding)


def load_config(path, collect_files=False):
    """
    Load configuration from given file replace, dictionaries by dotdict
    uses pyyaml-tags namely for:
    include tag:
        geometry: <% include(path="config_geometry.yaml")>
    """
    instance = YamlInclude.add_to_loader_class(loader_class=yaml.FullLoader, base_dir=os.path.dirname(path))
    cfg_dir = os.path.dirname(path)
    with open(path) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
    dd = dotdict.create(cfg)
    if collect_files:
        referenced = instance.included_files
        referenced.append(path)
        referenced.extend(collect_referenced_files(dd, ['.', cfg_dir]))
        dd['_file_refs'] = referenced
    return dd


def path_search(filename, path):
    for dir in path:
        if not isinstance(filename, str):
            continue
        full_name = os.path.join(dir, filename)
        if os.path.isfile(full_name):
            return [os.path.abspath(full_name)]
    return []

FilePath = NewType('FilePath', str)
def collect_referenced_files(cfg:dotdict, search_path:List[str]) -> List[FilePath]:
    referenced = []
    if isinstance(cfg, (dict, dotdict)):
        referenced = [collect_referenced_files(v, search_path) for v in cfg.values()]
    elif isinstance(cfg, (list, tuple)):
        referenced = [collect_referenced_files(v, search_path) for v in cfg]
    else:
        return path_search(cfg, search_path)
    # flatten
    return [i for l in referenced for i in l]
