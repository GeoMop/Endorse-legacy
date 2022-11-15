from endorse import common
import os
import time
import pytest

script_dir = script_dir = os.path.dirname(os.path.realpath(__file__))



@common.memoize
def file_func(f:common.File, out_file:str) -> common.File:
    print(f"file_func(f={f}, out_file={out_file})")
    with open(f.path, "r") as ff:
        content = ff.read()

    f_name = out_file
    with common.File.open(f_name, "w") as ff:      # need to use File.open that would check that the file does'nt exist
        ff.write(content)
        print(f"Appended: {out_file}", file=ff)

    return common.File.from_handle(ff)  # Have to create the File object explicitely after handle is closed.


def test_file_memoization():
    cache = common.EndorseCache.instance()
    cache.expire_all()

    input_file = "sandbox/memoize_file.txt"
    output_file = "sandbox/output_file.txt"
    with open(input_file, "w") as ff:
        ff.write(f"First line.")
    try:
        os.remove(output_file)
    except OSError:
        pass
    ####
    f = common.File(input_file)

    f1 = file_func(f, output_file)    # Create output_file.txt
    f2 = file_func(f, output_file)    # Re create output_file.txt, skipped.
    with pytest.raises(FileExistsError):
        f3 = file_func(f1, output_file)   # Trying to overwrite the created file. Should raise.



def func(a:int):
    print("Compute func.")
    time.sleep(2)
    return a * a

def test_memoization_reporting():
    common.EndorseCache.instance(report_calls=True)
