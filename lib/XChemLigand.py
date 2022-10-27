import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


def generate_steroisomers(input_smiles):
    """
    Function to generate all stereoisomers from a SMILES String. Only works on latest version of RDKit.
    """
    mol = Chem.MolFromSmiles(input_smiles)
    out_mols = AllChem.EnumerateStereoisomers(
        mol, AllChem.StereoEnumerationOptions(onlyUnassigned=False, tryEmbedding=True)
    )
    return [Chem.MolToSmiles(x, isomericSmiles=True) for x in out_mols]
