from Bio.PDB import *
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import math

def distSquared(coords1, coords2):
    return (coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 + (coords1[2]-coords2[2])**2

def getAtomsFromModel(modelResidues, refResidues):
    firstIter = True
    for residue in refResidues.get_residues():
        for atom in residue.get_atoms():
            coords = modelResidues[residue.id[1]][atom.get_id()].get_coord()
            if firstIter:
                modelCoords = np.array(coords)
                firstIter = False
            else:
                modelCoords = np.vstack([modelCoords,coords])
    return(modelCoords)

def transformCoords(refCoords, modelCoords):
    sup = SVDSuperimposer()
    sup.set(refCoords,modelCoords)
    sup.run()
    return sup.get_transformed()

def calculateRMSD(refCoords, modelCoords):
    rmsd = 0
    for i in range(len(refCoords)):
        rmsd += distSquared(refCoords[i], modelCoords[i])
    rmsd /= len(refCoords)
    rmsd = math.sqrt(rmsd)
    return(rmsd)

if __name__ == "__main__":
    parser = PDBParser()
    modelFile = parser.get_structure("Models", "R1107TS081.pdb")
    refFile = parser.get_structure("Reference", "R1107_reference.pdb")
    refResidues = refFile[0]["0"]
    refCoords = getAtomsFromModel(refResidues, refResidues)
    f = open("results/rmsdForModels.txt", "w")
    for model in modelFile.get_models():
        modelResidues = model["0"]
        modelCoords = getAtomsFromModel(modelResidues, refResidues)
        transformedCoords = transformCoords(refCoords, modelCoords)
        rmsd = calculateRMSD(refCoords, transformedCoords)
        f.write("Model " + str(model.id + 1) + "\n" + str(rmsd) + "\n")
