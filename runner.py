import subprocess
import os
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def collate(file, chunk_size=10):
    status_ = MPI.Status()

    ranks = {}
    for i in range(comm.Get_size() - 2):
        _ = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status_)
        ranks[i + 2] = []
    print("Found everyone and said hi.")
    assigner = 0

    print("Assigning")
    with open(file, 'r') as f:
        size = comm.Get_size() - 2
        # init phase
        for pos, line in (enumerate(f)):
            spl = line.split(' ')
            if len(spl) != 2:
                continue
            smile, name = spl[0].strip(), spl[1].strip()
            ranks[assigner + 2].append((pos, smile, name))
            assigner += 1
            assigner = assigner % (comm.Get_size() - 2)
            if pos != 0 and pos % (size * chunk_size) == 0:
                print(pos)
                break

        for k, v in ranks.items():
            print("Sending init data")
            comm.send(v, dest=k, tag=11)

        while True:
            try:
                _ = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status_)
                print("got request for more. ")
                source = status_.Get_source()
                new_data = []
                for new_pos, line in (enumerate(f)):
                    spl = line.split(' ')
                    if len(spl) != 2:
                        continue
                    smile, name = spl[0].strip(), spl[1].strip()
                    new_data.append((pos + new_pos, smile, name))
                    if new_pos != 0 and new_pos % chunk_size == 0:
                        print(new_pos, pos + new_pos)
                        pos = pos + new_pos
                        comm.send(new_data, dest=source, tag=11)
                        break
            except:
                pass


def setup_server(name):
    opentype = 'w'
    while True:
        try:
            with open(name, opentype, buffering=1) as f:
                f.write("pos,name,smiles,Dock,Dock_U,dbase,target\n")
                while True:
                    data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
                    for line in data:
                        f.write(line)
        except:
            opentype = 'a'


def worker(path_root, dbase_name, target_name, buffer_size, receptor_file=None):
    buffer = []

    while True:
        comm.send('need data', dest=1, tag=11)
        data = comm.recv(tag=11)

        with open("tmp/" + str(rank) + ".csv", 'w') as f:
            for pos in data:
                pos, smiles, name = pos
                path = path_root + str(pos) + "/"
                f.write(",".join([str(smiles), receptor_file, path, dbase_name, target_name, str(pos), name]) + '\n')

        try:
            try:
                call_string = ['python', 'theta_dock.py', "tmp/" + str(rank) + ".csv"]
                byteOutput = subprocess.check_output(call_string, shell=False)
                byteOutput = byteOutput.decode('UTF-8').rstrip() + "\n"
            except subprocess.CalledProcessError as e:
                print("Error in ls -a:\n", e.output)
                continue
            except subprocess.TimeoutExpired as e:
                print("Error in ls -a:\n", e.output)
                continue
            except Exception as e:
                print("Error rank", rank, e)
                continue

            res = byteOutput
            # score, res = interface_functions.RunDocking_(smiles, receptor_file, path, dbase_name, target_name,
            #                                              pos=pos, write=True,
            #                                              receptor_file=receptor_file, name=name,
            #                                              docking_only=docking_only)
            if res is not None:
                buffer.append(res)
            if len(buffer) > buffer_size:
                comm.send(buffer, dest=0, tag=11)
                buffer = []

        except KeyboardInterrupt as e:
            print("Error rank", rank, e)
            exit()
        except subprocess.CalledProcessError as e:
            print("Error rank", rank, e)
        except IndexError as e:
            print("Error rank", rank, e)
        except RuntimeError as e:
            print("Error rank", rank, e)
        except Exception as e:
            print("Error rank", rank, e)


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles', type=str, required=True)
    parser.add_argument('--path', type=str, required=True)
    parser.add_argument('--target_name', type=str, required=True)
    parser.add_argument('--dbase_name', type=str, required=True)
    parser.add_argument('--receptor_file', type=str, required=True)
    parser.add_argument('--chunk_size', type=int, required=False, default=10)
    parser.add_argument('--buffer_size', type=int, required=False, default=2)

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    path_root = args.path

    if rank == 0:
        if not os.path.exists(path_root):
            os.mkdir(path_root)
        setup_server(args.target_name + "_" + args.dbase_name + "_output.csv")
    elif rank == 1:
        collate(args.smiles, args.chunk_size)
    else:
        worker(path_root + "/rank", args.dbase_name, args.target_name, args.buffer_size, receptor_file=args.receptor_file)
