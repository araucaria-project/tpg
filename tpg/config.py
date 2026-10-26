"""Configuration loading shared by tpg and sfs.

Kept free of heavy imports (astropy, pyaraucaria, ctc), so `sfs` can load the config without
importing the whole plan generator.
"""

import os

import yaml

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))


def load_config(log=None):
    """Packaged config.yaml merged with the local override file (if it exists).

    @param log: optional callable taking a message, called for every file loaded
    @return: merged config dict
    """
    config_file = os.path.join(PACKAGE_DIR, "config.yaml")
    with open(config_file, 'r') as cfg_file:
        cfg = yaml.safe_load(cfg_file)
    if log:
        log(f'config loaded from {config_file}')

    local_config_file = os.path.expanduser(cfg["local_config"])
    if os.path.exists(local_config_file):
        with open(local_config_file) as f:
            local_cfg = yaml.safe_load(f) or {}
        deep_update(cfg, local_cfg)
        if log:
            log(f'local config loaded from {local_config_file}')

    return cfg


def resolve_path(path):
    """Path as given if it exists, else relative to the package dir, else to the repo root.

    Lets the config hold both absolute production paths and repo-relative ones
    (e.g. `master_files/...`). Returns the path unchanged if none of them exists.
    """
    path = os.path.expanduser(path)
    for base in ("", PACKAGE_DIR, os.path.dirname(PACKAGE_DIR)):
        candidate = os.path.join(base, path) if base else path
        if os.path.isfile(candidate):
            return candidate
    return path


def deep_update(base: dict, extra: dict) -> dict:
    # merge extra into base in place, key by key at every level - a local config
    # overrides single fields of a telescope block instead of replacing the whole block
    for k, v in extra.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        else:
            base[k] = v
    return base
