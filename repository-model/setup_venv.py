import venv
import os
import pip
import importlib

env_dir='env'

print("(SETUP ENV): Setup virtual environment for the project.")
print("(SETUP ENV): Creating python environment 'env'.")
venv.create(env_dir)

#source ./load_modules.sh

activate_script = os.path.abspath(os.path.join(env_dir, "bin", "activate_this.py"))
with open(activate_script) as fpy:
    exec(f.read())

#print("(SETUP ENV): 'env' created python and pip in use:")
#python --version
#which python
#which pip

print("(SETUP ENV): Installing necessary packages.")
packages=['pyyaml', 'attrs', 'numpy', 'matplotlib', 'bgem']
modules={}
for p in packages:
    pip.main(["install", "--prefix", env_dir, p])
    modules[p] = importlib.import_module(p)

