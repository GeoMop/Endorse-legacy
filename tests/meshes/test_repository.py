# Test src/repository_mesh.py
import yaml

def main():
    conf_file = "../test_data/config.yaml"
    with open(conf_file) as f:
        conf = yaml.safe_load(f)