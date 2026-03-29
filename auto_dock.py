import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy




def prepare_mol(mol):
    mol = Chem.AddHs(mol)

    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    elif not mol.GetConformer().Is3D():
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    return mol


def get_name(mol, idx):
    if mol.HasProp("name"):
        name = mol.GetProp("name")
        if name:
            return name
    return f"ligand_{idx}"


def mol_to_pdbqt(mol, out_path, flexible_rings=False):
    preparator = MoleculePreparation()

    if flexible_rings:
        preparator.min_ring_size = 6

    molsetup = preparator.prepare(mol)[0]

    result = PDBQTWriterLegacy.write_string(molsetup)
    pdbqt_string = result[0]

    with open(out_path, "w") as f:
        f.write(pdbqt_string)


def process_sdf(sdf_path, out_dir, flexible_rings=False):
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)

    for i, mol in enumerate(supplier):
        if mol is None:
            continue

        mol = prepare_mol(mol)
        name = get_name(mol, i)

        out_path = os.path.join(out_dir, f"{name}.pdbqt")
        mol_to_pdbqt(mol, out_path, flexible_rings)


def process_smi(smi_path, out_dir, flexible_rings=False):
    with open(smi_path) as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if not parts:
                continue

            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"ligand_{i}"

            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            mol = prepare_mol(mol)

            out_path = os.path.join(out_dir, f"{name}.pdbqt")
            mol_to_pdbqt(mol, out_path, flexible_rings)


def main(input_file, output_dir, flexible_rings=False):
    os.makedirs(output_dir, exist_ok=True)

    if input_file.endswith(".sdf"):
        process_sdf(input_file, output_dir, flexible_rings)

    elif input_file.endswith(".smi") or input_file.endswith(".smiles"):
        process_smi(input_file, output_dir, flexible_rings)

    else:
        raise ValueError("Unsupported file format")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("-o", "--output", default="pdbqt_out")
    parser.add_argument("--flexible-rings-6plus", action="store_true")

    args = parser.parse_args()

    main(
        args.input_file,
        args.output,
        flexible_rings=args.flexible_rings_6plus,
    )