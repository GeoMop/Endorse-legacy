import os
import sys
import shutil
import argparse
import subprocess
import tarfile

from argparse import RawTextHelpFormatter


def mprint(*margs, **mkwargs):
    print(*margs, file=sys.stdout, flush=True, **mkwargs)


def oscommand(command_string):
    mprint(command_string)
    mprint(os.popen(command_string).read())


def create_ssh_agent():
    mprint("creating ssh agent...")
    p = subprocess.Popen('ssh-agent -s',
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True, universal_newlines=True)
    outinfo, errinfo = p.communicate('ssh-agent cmd\n')
    # print(outinfo)

    lines = outinfo.split('\n')
    for line in lines:
        # trim leading and trailing whitespace
        line = line.strip()
        # ignore blank/empty lines
        if not line:
            continue
        # break off the part before the semicolon
        left, right = line.split(';', 1)
        if '=' in left:
            # get variable and value, put into os.environ
            varname, varvalue = left.split('=', 1)
            mprint("setting variable from ssh-agent:", varname, "=", varvalue)
            os.environ[varname] = varvalue


def dirs_tar(scratch_dir, tar_filepath, dirs_to_scratch):
    """
    :param scratch_dir: path to scratch directory
    :param tar_filepath: path to the tar file that is currently being created
    :param dirs_to_scratch: list of dirs that are being compressed and moved to the scratch
    """
    tar_names = []
    with tarfile.open(tar_filepath, mode='w:xz') as archive:
        for path in dirs_to_scratch:
            if os.path.isfile(path) and tarfile.is_tarfile(path):
                filename = os.path.basename(path)
                print("shutil copy path ", path)
                print("scratch dir ", scratch_dir)
                shutil.copy(path, scratch_dir)
                tar_names.append(filename)
            archive.add(path, recursive=True)
    print("tar names ", tar_names)
    return tar_names


def mpiexec_value(arg_val):
    if "mpiexec" not in arg_val and arg_val != "None":
        raise argparse.ArgumentTypeError("Invalid value. Set path to mpiexec or None for not using mpiexec")
    return arg_val


if __name__ == "__main__":

    mprint("================== singularity_exec.py START ==================")
    script_dir = os.getcwd()

    parser = argparse.ArgumentParser(
        description='Auxiliary executor for parallel programs running inside (Singularity) container under PBS.',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d', '--debug', action='store_true',
                        help='use testing files and print the final command')
    parser.add_argument('-i', '--image', type=str, required=True,
                        help='Singularity SIF image or Docker image (will be converted to SIF)')
    parser.add_argument('-B', '--bind', type=str, metavar="PATH,...", default="", required=False,
                        help='comma separated list of paths to be bind to Singularity container')
    parser.add_argument('-m', '--mpiexec', type=mpiexec_value, metavar="PATH", default="", required=False,
                        help="path (inside the container) to mpiexec to be run, default is 'mpiexec'")
    parser.add_argument('-s', '--scratch_copy', type=str, metavar="PATH", default="", required=False,
                        help='''
                        directory path, its content will be copied to SCRATCHDIR;
                        ''')
    parser.add_argument('-ds', '--dirs-to-scratch', nargs='+', default=[])
    # if file path, each user defined path inside the file will be copied to SCRATCHDIR
    parser.add_argument('prog', nargs=argparse.REMAINDER,
                        help='''
                        mpiexec arguments and the executable, follow mpiexec doc:
                        "mpiexec args executable pgmargs [ : args executable pgmargs ... ]"

                        still can use MPMD (Multiple Program Multiple Data applications):
                        -n 4 program1 : -n 3 program2 : -n 2 program3 ...
                        ''')

    # create the parser for the "prog" command
    # parser_prog = parser.add_subparsers().add_parser('prog', help='program to be run and all its arguments')
    # parser_prog.add_argument('args', nargs="+", help="all arguments passed to 'prog'")

    # parser.print_help()
    # parser.print_usage()
    args = parser.parse_args()

    print("args ", args)
   

    # get debug variable
    debug = args.debug

    # get program and its arguments
    prog_args = args.prog[1:]

    # get program and its arguments, set absolute path
    if os.path.isfile(args.image):
        image = os.path.abspath(args.image)
    elif args.image.startswith('docker://'):
        image = args.image
    else:
        raise Exception("Invalid image: not a file nor docker hub link: " + args.image)

    mprint("Hostname: ", os.popen('hostname').read())
    # mprint("os.environ", os.environ)

    ###################################################################################################################
    # Process node file and setup ssh access to given nodes.
    ###################################################################################################################

    # get nodefile, copy it to local dir so that it can be passed into container mpiexec later
    if debug:
        node_file = "testing_hostfile"
    else:
        mprint("getting host file...")
        orig_node_file = os.environ['PBS_NODEFILE']
        node_file = os.path.join(script_dir, os.path.basename(orig_node_file))
        shutil.copy(orig_node_file, node_file)
        # mprint(os.popen("ls -l").read())

    # Get ssh keys to nodes and append it to $HOME/.ssh/known_hosts
    ssh_known_hosts_to_append = []
    if debug:
        # ssh_known_hosts_file = 'testing_known_hosts'
        ssh_known_hosts_file = 'xxx/.ssh/testing_known_hosts'
    else:
        assert 'HOME' in os.environ
        ssh_known_hosts_file = os.path.join(os.environ['HOME'], '.ssh/known_hosts')

    mprint("host file name:", ssh_known_hosts_file)

    ssh_known_hosts = []
    if os.path.exists(ssh_known_hosts_file):
        with open(ssh_known_hosts_file, 'r') as fp:
            ssh_known_hosts = fp.readlines()
    else:
        mprint("creating host file...")
        dirname = os.path.dirname(ssh_known_hosts_file)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

    mprint("reading host file...")
    with open(node_file) as fp:
        node_names_read = fp.read().splitlines()
        # remove duplicates
        node_names = list(dict.fromkeys(node_names_read))

    mprint("connecting nodes...")
    for node in node_names:
        # touch all the nodes, so that they are accessible also through container
        os.popen('ssh ' + node + ' exit')
        # add the nodes to known_hosts so the fingerprint verification is skipped later
        # in shell just append # >> ~ /.ssh / known_hosts
        # or sort by 3.column in shell: 'sort -k3 -u ~/.ssh/known_hosts' and rewrite
        ssh_keys = os.popen('ssh-keyscan -H ' + node).readlines()
        ssh_keys = list((line for line in ssh_keys if not line.startswith('#')))
        for sk in ssh_keys:
            splits = sk.split(" ")
            if not splits[2] in ssh_known_hosts:
                ssh_known_hosts_to_append.append(sk)

    mprint("finishing host file...")
    with open(ssh_known_hosts_file, 'a') as fp:
        fp.writelines(ssh_known_hosts_to_append)

    # mprint(os.environ)
    create_agent = 'SSH_AUTH_SOCK' not in os.environ
    if not create_agent:
        create_agent = os.environ['SSH_AUTH_SOCK'] == ''

    if create_agent:
        create_ssh_agent()
    assert 'SSH_AUTH_SOCK' in os.environ
    assert os.environ['SSH_AUTH_SOCK'] != ""

    ###################################################################################################################
    # Create Singularity container commands.
    ###################################################################################################################

    mprint("assembling final command...")

    scratch_dir_path = None
    if 'SCRATCHDIR' in os.environ:
        scratch_dir_path = os.environ['SCRATCHDIR']
        mprint("Using SCRATCHDIR:", scratch_dir_path)

        mprint("copying to SCRATCHDIR on all nodes...")
        username = os.environ['USER']
        # get source files
        source = None
        if os.path.isdir(args.scratch_copy):
            # source = args.scratch_copy + "/."
            # paths = [os.path.join(args.scratch_copy,fp) for fp in os.listdir(args.scratch_copy)]
            # source = ' '.join(paths)
            source = args.scratch_copy
        else:
            raise Exception("--scratch_copy argument is not a valid directory: " + args.scratch_copy)
            # with open(args.scratch_copy) as fp:
            #   paths = fp.read().splitlines()
            #   source = ' '.join(paths)

        if source is None or source is []:
            mprint(args.scratch_copy, "is empty")

        # create tar
        source_tar_filename = 'scratch.tar'
        source_tar_filepath = os.path.join(script_dir, source_tar_filename)
        command = ' '.join(['cd', source, '&&', 'tar -cvf', source_tar_filepath, '.', '&& cd', script_dir])
        oscommand(command)

        dirs_to_scratch = args.dirs_to_scratch
        tar_names = []
        if len(dirs_to_scratch) > 0:
            source_dirs_tar_filename = 'dirs_scratch.tar.xz'
            source_dirs_tar_filepath = os.path.join(script_dir, source_dirs_tar_filename)
            tar_names = dirs_tar(scratch_dir_path, source_dirs_tar_filepath, dirs_to_scratch)
            tar_names.append(source_dirs_tar_filename)

        for node in node_names:
            destination_name = username + "@" + node
            destination_path = destination_name + ':' + scratch_dir_path
            print("scratch dir path ", scratch_dir_path)
            print("destination path ", destination_path)
            command = ' '.join(['scp', source_tar_filepath, destination_path])
            oscommand(command)

            # command = ' '.join(['ssh', destination_name, 'cd', scratch_dir_path, '&&', 'tar --strip-components 1 -xf', source_tar_filepath, '-C /'])

            command_list = ['ssh', destination_name, '"cd', scratch_dir_path, '&&', 'tar -xf', source_tar_filename,
                                '&&', 'rm ', source_tar_filename]

            if len(tar_names) > 0:
                for tar_name in tar_names:
                    command = command_list.extend(['&&', 'tar -xf', tar_name, '&&', 'rm ', tar_name])

            command_list.append('"')

            command = ' '.join(command_list)
            print("command ", command)

            oscommand(command)

        # remove the scratch tar
        oscommand(' '.join(['rm', source_tar_filename]))

    # A] process bindings, exclude ssh agent in launcher bindings
    bindings = "-B " + os.environ['SSH_AUTH_SOCK']
    # possibly add current dir to container bindings
    # bindings = bindings + "," + script_dir + ":" + script_dir
    bindings_in_launcher = ""
    if args.bind != "":
        bindings = bindings + "," + args.bind
        bindings_in_launcher = "-B " + args.bind

    if scratch_dir_path:
        bindings = bindings + "," + scratch_dir_path
        if args.bind == "":
            bindings_in_launcher = "-B " + scratch_dir_path
        else:
            bindings_in_launcher = bindings_in_launcher + "," + scratch_dir_path

    sing_command = ' '.join(['singularity', 'exec', '--nv', bindings, image])
    sing_command_in_launcher = ' '.join(['singularity', 'exec', '--nv', bindings_in_launcher, image])

    mprint('sing_command:', sing_command)
    mprint('sing_command_in_ssh:', sing_command_in_launcher)

    # B] prepare node launcher script
    mprint("creating launcher script...")
    launcher_path = os.path.join(script_dir, "launcher.sh")
    launcher_lines = [
        '#!/bin/bash',
        '\n',
        'echo $(hostname) >> launcher.log',
        'echo $(pwd) >> launcher.log',
        'echo $@ >> launcher.log',
        'echo "singularity container: $SINGULARITY_NAME" >> launcher.log',
        '\n',
        'ssh $1 $2 ' + sing_command_in_launcher + ' ${@:3}',
    ]
    with open(launcher_path, 'w') as f:
        f.write('\n'.join(launcher_lines))
    oscommand('chmod +x ' + launcher_path)

    if args.mpiexec != "None":
        # C] set mpiexec path inside the container
        # if container path to mpiexec is provided, use it
        # otherwise try to use the default
        mpiexec_path = "mpiexec"
        if args.mpiexec != "":
            mpiexec_path = args.mpiexec

        # test_mpiexec = os.popen(sing_command + ' which ' + 'mpiexec').read()
        # # test_mpiexec = os.popen('singularity exec docker://flow123d/geomop:master_8d5574fc2 which flow123d').read()
        # mprint("test_mpiexec: ", test_mpiexec)
        # if mpiexec_path == "":
        #     raise Exception("mpiexec path '" + mpiexec_path + "' not found in container!")

        # D] join mpiexec arguments
        mpiexec_args = " ".join([mpiexec_path, '-f', node_file, '-launcher-exec', launcher_path])

        # F] join all the arguments into final singularity container command
        final_command = " ".join([sing_command, mpiexec_args, *prog_args])
    else:
        print("prog args ", prog_args)
        #p_args_split = prog_args[0].split()
        #p_args_split[2] = os.path.join(scratch_dir_path, p_args_split[2][1:])
        #p_args_split[3] = os.path.join(scratch_dir_path, p_args_split[3][1:])
        #print("p_args_split[2] exists ", os.path.exists(p_args_split[2]))
        #print("p_args_split[3] exists ", os.path.exists(p_args_split[3]))


        #prog_args = [" ".join(p_args_split)]
        print("prog args scratch", prog_args)


        final_command = " ".join([sing_command, *prog_args])


    ###################################################################################################################
    # Final call.
    ###################################################################################################################
    if scratch_dir_path:
        mprint("Entering SCRATCHDIR:", scratch_dir_path)
        os.chdir(scratch_dir_path)

    mprint("current directory:", os.getcwd())
    # mprint(os.popen("ls -l").read())
    mprint("final command:")
    mprint("=================== singularity_exec.py END ===================")
    if not debug:
        mprint("================== Program output START ==================")
        os.system(final_command)
        mprint("=================== Program output END ===================")
