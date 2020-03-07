import sys

from impress_md import interface_functions

if __name__ == '__main__':
    smiles, receptor_file, path, dbase_name, target_name, pos, name = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[
        4], sys.argv[5], sys.argv[6], sys.argv[7]
    interface_functions.RunDocking_(smiles, receptor_file, path, dbase_name, target_name,
                                    pos=pos, write=True,
                                    receptor_file=receptor_file, name=name,
                                    docking_only=True)
