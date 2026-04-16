# Protein Structure Analysis Toolkit

> Bioinformatics isn't just about sequences. Proteins fold into three-dimensional structures, and those structures determine function. This toolkit gives you the essential structural analyses from PDB files.

## What Each Analysis Tells You

**RMSD**: How similar are two structures? < 1 angstrom = nearly identical. 1-3 = similar fold. > 3 = significant changes. Use this to compare wild-type vs. mutant, or AlphaFold prediction vs. experimental structure.

**Contact Map**: Which residues are physically close in 3D? Two residues far apart in sequence might be neighbours in the fold. Reveals domain architecture, binding interfaces, and validates predictions.

**Ramachandran Plot**: Are backbone angles physically reasonable? Residues in disallowed regions suggest structural errors. First quality check for any new structure.

**B-factors**: How much does each residue move? High = flexible (often loops involved in binding). Low = rigid core.

## Usage
```python
from structure_tools import ProteinAnalyser
pa = ProteinAnalyser('protein.pdb')
pa.ramachandran_plot('rama.png')
pa.bfactor_plot('bfactors.png')
```
