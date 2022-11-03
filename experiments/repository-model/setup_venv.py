import pip
import venv
import os
import pip
import importlib
import shutil

env_dir='env'
print("(SETUP ENV): Setup virtual environment for the project.")
print("(SETUP ENV): Creating python environment 'env'.")
builder = venv.EnvBuilder(system_site_packages=False, clear=False, symlinks=False, upgrade=False, with_pip=True)
builder.create(env_dir)


#source ./load_modules.sh
# TODO: find a way how to do it directly from python without using shell

env_python = os.path.join(env_dir, "bin", "python")
#shutil.copyfile("activate_this.py", activate_script)
#with open(activate_script) as fpy:
    #exec(fpy.read())

#print("(SETUP ENV): 'env' created python and pip in use:")
#python --version
#which python
#which pip

print("(SETUP ENV): Installing necessary packages.")
import sys, subprocess
packages=['pyyaml', 'attrs', 'numpy', 'matplotlib', 'bgem']
subprocess.check_call([env_python, '-m', 'pip', 'install'] + packages)

    

