from mpi4py import MPI
from impress_md import interface_functions
import subprocess

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def collate(file, chunk=15):
    status_ = MPI.Status()

    print("loadintg data")
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
            if pos != 0 and pos % (size * chunk) == 0:
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
                    if new_pos != 0 and new_pos % chunk == 0:
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


def worker(path_root, dbase_name, target_name, docking_only=False, receptor_file=None):
    from rdkit import Chem
    from rdkit.Chem import Lipinski, Descriptors
    struct = "input/"
    docker, recept = interface_functions.get_receptr(receptor_file=receptor_file)
    mols_docked = 0

    buffer = []

    while True:
        comm.send('hi', dest=1, tag=11)
        data = comm.recv(tag=11)
        for pos in data:
            pos, smiles, name = pos

            # try:
            #     mol = Chem.MolFromSmiles(smiles)
            #     if mol is None:
            #         print("passing on mol", smiles)
            #         continue
            #     else:
            #         if mol.GetNumAtoms() < 12 or Descriptors.MolWt(mol) > 500:
            #             print("passing on mol", smiles)
            #             continue
            # except:
            #     print("passing on mol", smiles)
            #     continue

            path = path_root + str(pos) + "/"
            try:
                score, res = interface_functions.RunDocking_(smiles, struct, path, dbase_name, target_name, pos=pos,
                                                             dock_obj=docker, write=True, recept=recept,
                                                             receptor_file=receptor_file, name=name,
                                                             docking_only=docking_only)
                mols_docked += 1

                if docking_only:
                    if res is not None:
                        buffer.append(res)
                    if len(buffer) > 5:
                        comm.send(buffer, dest=0, tag=11)
                        buffer = []

            except KeyboardInterrupt:
                exit()
            except subprocess.CalledProcessError as e:
                print("Error rank", rank, e)
            except IndexError as e:
                print("Error rank", rank, e)
            except RuntimeError as e:
                print("Error rank", rank, e)
            except:
                pass


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--dock_only', action='store_true')
    parser.add_argument('--smiles', type=str, required=True)
    parser.add_argument('--path', type=str, required=True)
    parser.add_argument('--target_name', type=str, required=True)
    parser.add_argument('--dbase_name', type=str, required=True)
    parser.add_argument('--receptor_file', type=str, required=True)
    return parser.parse_args()


if __name__ == '__main__':
    import os

    args = get_args()

    path_root = args.path

    if rank == 0:
        if not os.path.exists(path_root):
            os.mkdir(path_root)
        setup_server(args.path + "_" + args.target_name + "_" + args.dbase_name + "_output.csv")
    elif rank == 1:
        collate(args.smiles)
    else:
        worker(path_root + "/rank", args.dbase_name, args.target_name, docking_only=args.dock_only,
               receptor_file=args.receptor_file)
