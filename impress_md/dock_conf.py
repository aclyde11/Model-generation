from openeye import oechem, oedocking

def PrepareReceptor(pdb,padding=4,outpath=""):
    """
    Prepares a receptor from a pdb with a crystalized ligand
    Padding controls the docking region.
    If outpath is given, PrepareReceptor will write an openeye binary (oeb) of the receptor structure. This will be faster than rebuilding the receptor every time.
    """
    print("STOP CALLING THIS FUNCTION")
    exit()
    com = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    if ifs.open(pdb):
        oechem.OEReadPDBFile(ifs, com)
        ifs.close()

    """
    Sorry, this requires some explanation. Openeye wasn't recognizing the previously docked ligand, so I tried to find other ways.
    The next blocks of code take our system and split it based on its connected components, for which its REQUIRED that our protein
      only has a single chain. It assumes that the last component is the ligand. It then creates the ligand (lig) and protein (prot)
      as separate molecules. Next, it finds the minimum and maximum 3D coordinates of the current ligand and produces a box around
      it with the specified padding. Finally it uses this box to create a 'receptor' object into which ligands can be docked.
    Only the receptor is returned.
    Openeye's docking shouldn't be this involved, but I couldn't get it to run the typical 'hybrid' docking without error.
    """
    oechem.OEDetermineConnectivity(com)
    nparts, connect = oechem.OEDetermineComponents(com)
    if(nparts != 2):
        print("ERR in dock_conf::prepareReceptor. PDB doesn't have 2 connected components")
        exit()
        ## TODO: What is a good way to catch errors?
    # Get apo
    pred = oechem.OEPartPredAtom(connect)
    pred.SelectPart(nparts)
    lig = oechem.OEGraphMol()
    oechem.OESubsetMol(lig, com, pred)
    print(lig)
    
    # Get protein
    pred = oechem.OEPartPredAtom(connect)
    pred.SelectPart(1)
    prot = oechem.OEGraphMol()
    oechem.OESubsetMol(prot, com, pred)
    
    # Get box dimensions by iterating over ligand
    x_min = y_min = z_min = float('inf')
    x_max = y_max = z_max = -float('inf')
    crd = lig.GetCoords()
    print("CRD", crd)
    for atm in crd:
        x,y,z = crd[atm]
        if x < x_min:
            x_min = x
        if y < y_min:
            y_min = y
        if z < z_min:
            z_min = z
        if x > x_max:
            x_max = x
        if y > y_max:
            y_max = y
        if z > z_max:
            z_max = z
    x_min -= padding
    y_min -= padding
    z_min -= padding
    x_max += padding
    y_max += padding
    z_max += padding
    print(x_min,y_min,z_max, y_max)
    # Now prepare the receptor
    receptor = oechem.OEGraphMol()
    box = oedocking.OEBox()
    box.Setup(x_max, y_max, z_max, x_min, y_min, z_min)
    oedocking.OEMakeReceptor(receptor, prot, box)
    
    if not outpath == "":
        oedocking.OEWriteReceptorFile(receptor,f'{outpath}/receptor.oeb')
    return receptor

def PrepareReceptorFromBinary(filename):
    """
    If a receptor.oeb has been created, use this version.
    If not, run 'PrepareReceptor' and a receptor.oeb will be created.
    """
    receptor = oechem.OEGraphMol()
    oedocking.OEReadReceptorFile(receptor,filename)
    return receptor

def DockConf(pdb_file, mol, MAX_POSES = 1, dock=None):
    if dock is None:
        receptor = oechem.OEGraphMol()
        oedocking.OEReadReceptorFile(receptor, pdb_file)
        dock = oedocking.OEDock(oedocking.OEScoreType_Chemgauss4, oedocking.OESearchResolution_High)
        dock.Initialize(receptor)
    else:
        receptor = None
    lig = DockConf_(dock, mol, MAX_POSES=MAX_POSES)
    return lig, receptor

def DockConf_(dock, mol, lig, MAX_POSES = 1, receptor_filename="Not Available"):
    err = dock.DockMultiConformerMolecule(lig,mol,MAX_POSES)
    if (err != oedocking.OEDockingReturnCode_Success):
        print("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(err))
    sdtag = dock.GetName()
    oedocking.OESetSDScore(lig, dock, sdtag)
    dock.AnnotatePose(lig)
    oechem.OESetSDData(lig, "receptor", receptor_filename)
    return lig

def WriteStructures(receptor, lig, apo_path, lig_path, com_path=None, oe=False):
    ofs = oechem.oemolostream()
    success = True
    if ofs.open(apo_path):
        oechem.OEWriteMolecule(ofs,receptor)
        ofs.close()
    else:
        success = False
    # TODO: If MAX_POSES != 1, we should select the top pose to save
    if ofs.open(lig_path):
        oechem.OEWriteMolecule(ofs,lig)
        ofs.close()
    else:
        success = False

    if com_path is not None and ofs.open(com_path):
        oechem.OEWriteMolecule(ofs, receptor)
        oechem.OEWriteMolecule(ofs, lig)
        ofs.close()


    if oe:
        import re
        ofs = oechem.oemolostream()
        success = True
        if ofs.open(re.sub(".pdb", ".oeb", apo_path)):
            oechem.OEWriteMolecule(ofs, receptor)
            ofs.close()
        else:
            success = False
        # TODO: If MAX_POSES != 1, we should select the top pose to save
        if ofs.open(re.sub(".pdb", ".oeb", lig_path)):
            oechem.OEWriteMolecule(ofs, lig)
            ofs.close()
        else:
            success = False

        if com_path is not None and ofs.open(re.sub(".pdb", ".oeb", com_path)):
            oechem.OEWriteMolecule(ofs, receptor)
            oechem.OEWriteMolecule(ofs, lig)
            ofs.close()

    return success

### Returns an array of length MAX_POSES from above. This is the range of scores
def LigandScores(dock, lig):
    return [ dock.ScoreLigand(conf) for conf in lig.GetConfs() ]

def BestDockScore(dock, lig):
    return LigandScores(dock,lig)[0]

def ScoreRange(dock,lig):
    tmp = LigandScores(dock,lig)
    return tmp[0],tmp[-1]
