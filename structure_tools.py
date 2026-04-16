#!/usr/bin/env python3
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa, PPBuilder

class ProteinAnalyser:
    def __init__(self, pdb_file, chain_id='A'):
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('protein', pdb_file)
        self.chain = self.structure[0][chain_id]

    def get_ca_atoms(self):
        return [res['CA'] for res in self.chain if is_aa(res) and 'CA' in res]

    def compute_rmsd(self, other_pdb, other_chain='A'):
        other = self.parser.get_structure('other', other_pdb)
        fixed = self.get_ca_atoms()
        moving = [res['CA'] for res in other[0][other_chain] if is_aa(res) and 'CA' in res]
        n = min(len(fixed), len(moving))
        sup = Superimposer()
        sup.set_atoms(fixed[:n], moving[:n])
        return sup.rms

    def contact_map_plot(self, output, threshold=8.0):
        ca = self.get_ca_atoms()
        n = len(ca)
        cmap = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                if ca[i] - ca[j] < threshold:
                    cmap[i,j] = cmap[j,i] = 1
        fig, ax = plt.subplots(figsize=(8,8))
        ax.imshow(cmap, cmap='Blues', origin='lower')
        ax.set_title(f'Contact Map ({threshold}A)')
        plt.savefig(output, dpi=150); plt.close()

    def ramachandran_plot(self, output):
        ppb = PPBuilder()
        angles = [(np.degrees(phi), np.degrees(psi)) for pp in ppb.build_peptides(self.chain) for phi, psi in pp.get_phi_psi_list() if phi and psi]
        if not angles: return
        phis, psis = zip(*angles)
        fig, ax = plt.subplots(figsize=(8,8))
        ax.scatter(phis, psis, s=10, alpha=0.6)
        ax.set_xlim(-180,180); ax.set_ylim(-180,180)
        ax.set_title('Ramachandran Plot')
        plt.savefig(output, dpi=150); plt.close()

    def bfactor_plot(self, output):
        residues = [res for res in self.chain if is_aa(res)]
        nums = [res.id[1] for res in residues]
        bfacs = [np.mean([a.get_bfactor() for a in res]) for res in residues]
        fig, ax = plt.subplots(figsize=(12,4))
        ax.plot(nums, bfacs); ax.fill_between(nums, bfacs, alpha=0.3)
        ax.set_title('B-factor Profile')
        plt.savefig(output, dpi=150); plt.close()
