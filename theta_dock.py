import sys
import pandas as pd
from impress_md import interface_functions

if __name__ == '__main__':
    file_name = sys.argv[1]
    ress = ""
    with open(file_name, 'r') as f:
        for line in f:
            smiles, receptor_file, path, dbase_name, target_name, pos, name = line.split(',')
            name = name.strip()
            _, res = interface_functions.RunDocking_(smiles, receptor_file, path, dbase_name, target_name,
                                            pos=pos, write=True,
                                            receptor_file=receptor_file, name=name,
                                            docking_only=True)
            ress += res
    print(ress)
