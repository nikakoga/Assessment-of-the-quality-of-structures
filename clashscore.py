import argparse
from Bio.PDB import *
import numpy as np
import math

radii = {
    "O": 1.52,
    "P": 1.8,
    "N": 1.55,
    "H": 1.2,
    "C": 1.7
}

def dist(coords1, coords2):
    return math.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 + (coords1[2]-coords2[2])**2)

def atomCount(chain):
    atoms = 0
    for atom in chain.get_atoms():
        atoms += 1
    return atoms

def numberOfClashes2Res(res1, res2,threshold):
    clashes = 0
    for atom1 in res1.get_atoms():
        for atom2 in res2.get_atoms():
            if not ((atom1.name == "O3'" and atom2.name == "P") or (atom2.name == "O3'" and atom1.name == "P")):
                if (atom1-atom2)<(radii[atom1.element] + radii[atom2.element] - threshold):
                    clashes+=1
    return clashes

def clashScore(PDBFile, threshold, chain_id, model_id):
    clashes = 0
    parser = PDBParser()
    try:
        structure = parser.get_structure("structure", PDBFile)
        if model_id == None:
            for i in structure.get_models():
                model = i
                break
        else:
            model = structure[model_id]
        if chain_id == None:
            for i in model:
                chain = i
                break
        else:
            chain = model[chain_id]
        for i, residue in enumerate(chain):
            for j in range(i+1, len(chain)):
                clashes += numberOfClashes2Res(chain[i+1], chain[j+1],threshold)
        return 1000*(clashes/atomCount(chain))
    except FileNotFoundError:
        print("Error! File not found. Please make sure that filename and path to file are correct.")
        return -1
    except KeyError:
        print("Error! Model number or chain code are incorrect.")
        return -1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",help="Input file in PDB format", required=True)
    parser.add_argument("-t","--threshold",help="Overlap threshold for clash in Angstrom. Default value: 0.4", type = float, default=0.4)
    parser.add_argument("-c","--chain",help="Chain code from PDB file. If not provided, program takes first chain from model.")
    parser.add_argument("-m","--model",help="Model number from PDB file. If not provided, program takes first model from file.")
    parser.add_argument("-o","--output",help="Output .txt file in which the clashscore will be saved. If not provided, clashscore will not be saved to your computer.")
    args = parser.parse_args()
    score = clashScore(args.input, args.threshold, args.chain, args.model)
    if args.output:
        if not args.output.endswith(".txt"):
            args.output += ".txt"
        file = open(args.output, "w")
        file.write("Clashscore: " + str(round(score,2)))
    print("Clashscore: " + str(round(score,2)))