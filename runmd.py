import sys
from impress_md import interface_functions

path = sys.argv[1]
interface_functions.RunMMGBSA_(path, path, gpu=True)
