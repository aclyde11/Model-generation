# Workflow 0
import argparse

import numpy as np
from mpi4py import MPI
from openeye import oechem

from impress_md import interface_functions

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
world_size = comm.Get_size()
import time
import signal

WORKTAG, DIETAG = 11, 13


class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help='input csv for smiles', required=True, type=str)
    parser.add_argument("-o", help='output file for data', required=True, type=str)
    parser.add_argument("-r", help='receptor file', required=True, type=str)
    parser.add_argument('-v', help='verbose (1 print every 1000, 2 print every thing)', type=int, choices=[0,1,2], default=1)
    parser.add_argument("-n", type=int, default=1)
    parser.add_argument("-l", type=str, default=None, required=False)

    return parser.parse_args()


def get_root_protein_name(file_name):
    return file_name.split("/")[-1].split(".")[0]


def get_smiles_col(col_names):
    return int(np.where(['smile' in s.lower() for s in col_names])[0][0])


def get_ligand_name_col(col_names):
    return int(np.where(['id' in s.lower() or 'title' in s.lower() or "name" in s.lower() for s in col_names])[0][0])


def master():
    mol = oechem.OEMol()

    ifs = oechem.oemolistream(input_smiles_file)
    ifs.SetConfTest(oechem.OEAbsCanonicalConfTest())

    rstart = time.time()
    for pos, mol in enumerate(ifs.GetOEMols()):
        rend = time.time()
        smiles = oechem.OEMol(mol)
        ligand_name = smiles.GetTitle()
        status = MPI.Status()
        waitstart = time.time()
        comm.recv(source=MPI.ANY_SOURCE, tag=WORKTAG, status=status)
        waitend = time.time()
        sstart = time.time()
        rank_from = status.Get_source()
        data = (pos, smiles, ligand_name)
        comm.send(data, dest=rank_from, tag=WORKTAG)
        send = time.time()
        if args.v == 1 and pos % 1000 == 0:
            print("sent", pos, "jobs")
        if pos % 10 == 0:
            print('master rtime', rend - rstart, 'stime', send - sstart, 'waitime', waitend - waitstart)
        rstart = time.time()
    for i in range(1, world_size):
        comm.send([], dest=i, tag=DIETAG)
    comm.Barrier()


def slave():
    docker, receptor = interface_functions.get_receptor(target_file, use_hybrid=use_hybrid,
                                                        high_resolution=high_resolution)

    comm.send([], dest=0, tag=11)
    poss = comm.recv(source=0, tag=WORKTAG)

    while True:
        pos, smiles, ligand_name = poss
        try:
            with timeout(seconds=120):
                dstart = time.time()

                score, res, ligand = interface_functions.RunDocking_conf(smiles,
                                                                         dock_obj=docker,
                                                                         pos=pos,
                                                                         name=ligand_name,
                                                                         target_name=pdb_name,
                                                                         force_flipper=force_flipper)
                dend = time.time()

                if args.v == 2:
                    print("RANK {}:".format(rank), res, end='')
                if ofs and ligand is not None:
                    wstart = time.time()
                    oechem.OEWriteMolecule(ofs, ligand)
                    wend = time.time()
                if rank in [3, 43, 218, 32]:
                    print("rank {} dtime".format(rank), dend - dstart, "wtime", wend - wstart)
        except TimeoutError:
            print("TIMEOUT", smiles, ligand_name)
            continue

        wstart = time.time()
        comm.send([], dest=0, tag=WORKTAG)
        status = MPI.Status()
        poss = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        wend = time.time()
        print("rank {} dtime".format(rank), "waittime", wend - wstart)

        if status.Get_tag() == DIETAG:
            break

    if ofs is not None:
        ofs.close()
    comm.Barrier()

if __name__ == '__main__':
    args = getargs()
    input_smiles_file = args.i
    target_file = args.r  # twenty of these
    basename, file_ending = ".".join(args.o.split(".")[:-1]), args.o.split(".")[-1]
    output_poses = basename + "_" + str(rank) + "." + file_ending
    pdb_name = get_root_protein_name(target_file)
    ## setting don't change

    use_hybrid = True
    force_flipper = False
    high_resolution = False
    ofs = oechem.oemolostream(output_poses)

    if rank == 0:
        master()
    else:
        slave()
    exit()
