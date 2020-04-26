# Workflow 0
import argparse

import numpy as np
from mpi4py import MPI
from openeye import oechem
import concurrent.futures
from impress_md import interface_functions
import queue
import threading
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
world_size = comm.Get_size()
import time
import copy
WORKTAG, DIETAG = 11, 13

# class timeout:
#     def __init__(self, seconds=1, error_message='Timeout'):
#         self.seconds = seconds
#         self.error_message = error_message
#
#     def handle_timeout(self, signum, frame):
#         raise TimeoutError(self.error_message)
#
#     def __enter__(self):
#         signal.signal(signal.SIGALRM, self.handle_timeout)
#         signal.alarm(self.seconds)
#
#     def __exit__(self, type, value, traceback):
#         signal.alarm(0)


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help='input csv for smiles', required=True, type=str)
    parser.add_argument("-o", help='output file for data', required=True, type=str)
    parser.add_argument("-r", help='receptor file', required=True, type=str)
    parser.add_argument('-v', help='verbose (1 print every 1000, 2 print every thing)', type=int, choices=[0,1,2], default=1)
    parser.add_argument("-n", type=int, default=1)
    parser.add_argument("-l", type=str, default=None, required=False)
    parser.add_argument("-w", type=int, default=3, required=False)
    parser.add_argument('-c', type=int, default=9, required=False)
    parser.add_argument('--queue_lim', type=int, default=120, required=False)
    return parser.parse_args()


def get_root_protein_name(file_name):
    return file_name.split("/")[-1].split(".")[0]


def get_smiles_col(col_names):
    return int(np.where(['smile' in s.lower() for s in col_names])[0][0])


def get_ligand_name_col(col_names):
    return int(np.where(['id' in s.lower() or 'title' in s.lower() or "name" in s.lower() for s in col_names])[0][0])

def loader(fname, queue, done, lim):
    mol = oechem.OEMol()

    ifs = oechem.oemolistream(fname)
    ifs.SetConfTest(oechem.OEAbsCanonicalConfTest())
    iterator = iter(enumerate(ifs.GetOEMols()))

    while True:
        data = []

        try:
            for _ in range(CHUNK):
                pos, mol = next(iterator)
                smiles = oechem.OEMol(mol)
                ligand_name = smiles.GetTitle()
                data.append((pos, smiles, ligand_name))
            queue.put(data)

            fors =0
            while queue.qsize() > lim:
                time.sleep(0.5)
                fors += 1
                if fors > 5 * 60 * 2:
                    break
        except StopIteration:
            if len(data) > 0:
                queue.put(data)
            break
    done.set()
    ifs.close()

def master():
    lim = args.queue_lim
    q = queue.Queue()
    done = threading.Event()
    done.clear()

    t = threading.Thread(target=loader, kwargs={'fname' : input_smiles_file, 'queue' : q, 'done' : done, 'lim' : lim})
    t.start()

    while not done.is_set() or not q.empty():
            #wait until asked
            status = MPI.Status()
            comm.recv(source=MPI.ANY_SOURCE, tag=WORKTAG, status=status)
            print('recevoed')
            if not q.empty():
                data = q.get()
            else:
                break
            pos = data[0][0]

            ## SEND
            rank_from = status.Get_source()
            comm.send(data, dest=rank_from, tag=WORKTAG)

            if args.v >= 1 and pos % 1000 == 0:
                print("sent", pos, "jobs")
    for i in range(1, world_size):
        comm.send([], dest=i, tag=DIETAG)
    t.join()

    comm.Barrier()

def run_dock_timeout(**kwargs):
    try:
        return interface_functions.RunDocking_conf(**kwargs)
    except Exception as e:
        print("unkown error",e)
    return None, None, None


def slave():
    ofs = oechem.oemolostream(output_poses)

    dockers = [interface_functions.get_receptor(target_file, use_hybrid=use_hybrid,high_resolution=high_resolution)[0] for _ in range(CHUNK)]

    reruns = {}

    while True:
        with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as executor:
            comm.send([], dest=0, tag=WORKTAG)
            status = MPI.Status()
            poss = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)

            if status.Get_tag() == DIETAG:
                break
            # Start the load operations and mark each future with its URL
            future_to_url = {executor.submit(run_dock_timeout, ligand=oechem.OEMol(smiles),
                                                                     dock_obj=dockers[i],
                                                                     pos=copy.copy(pos),
                                                                     name=copy.copy(ligand_name),
                                                                     target_name=copy.copy(pdb_name),
                                                                     force_flipper=copy.copy(force_flipper)) : False for i, (pos, smiles, ligand_name) in enumerate(poss)}
            future_to_url.update(reruns)
            try:
                for ct, future in enumerate(concurrent.futures.as_completed(future_to_url, timeout=90 * len(future_to_url) / WORKERS)):
                    # pos, smiles, ligand = future_to_url[future]
                    try:
                        _, res, ligand = future.result()
                        if args.v == 2:
                            print("RANK {}:".format(rank), res, end='')
                            # pass
                        if ofs and ligand is not None:
                            oechem.OEWriteMolecule(ofs, ligand)
                        future_to_url[future] = True

                    except Exception as exc:
                        pass
                    if len(future_to_url) - ct < (WORKERS):
                        raise(concurrent.futures.TimeoutError())

            except concurrent.futures.TimeoutError:
                reruns = {}
                for k, v in future_to_url.items():
                    if not v:
                        reruns[k] = v

    if ofs is not None:
        ofs.close()
    comm.Barrier()

if __name__ == '__main__':
    args = getargs()
    WORKERS = args.w
    CHUNK = args.c
    input_smiles_file = args.i
    target_file = args.r  # twenty of these
    basename, file_ending = ".".join(args.o.split(".")[:-1]), args.o.split(".")[-1]
    output_poses = basename + "_" + str(rank) + "." + file_ending
    pdb_name = get_root_protein_name(target_file)
    ## setting don't change

    use_hybrid = True
    force_flipper = False
    high_resolution = False

    if rank == 0:
        master()
    else:
        slave()
    exit()
