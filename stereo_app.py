import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from stmol import showmol
import py3Dmol
import numpy as np


# ==============================
# Allene axial chirality detector
# ==============================
def detect_allene_axes(mol):
    axes = []
    for bond1 in mol.GetBonds():
        if bond1.GetBondType() != Chem.BondType.DOUBLE:
            continue

        a1 = bond1.GetBeginAtom()
        a2 = bond1.GetEndAtom()

        for bond2 in a2.GetBonds():
            if bond2.GetIdx() == bond1.GetIdx():
                continue
            if bond2.GetBondType() != Chem.BondType.DOUBLE:
                continue

            a3 = bond2.GetOtherAtom(a2)

            if a1.GetSymbol()=="C" and a2.GetSymbol()=="C" and a3.GetSymbol()=="C":

                left_subs = [n.GetIdx() for n in a1.GetNeighbors() if n.GetIdx()!=a2.GetIdx()]
                right_subs = [n.GetIdx() for n in a3.GetNeighbors() if n.GetIdx()!=a2.GetIdx()]

                if len(left_subs)==2 and len(right_subs)==2:
                    if left_subs[0]!=left_subs[1] and right_subs[0]!=right_subs[1]:
                        axes.append((a1.GetIdx(), a2.GetIdx(), a3.GetIdx()))

    return axes


# ==============================
# Ra / Sa assignment
# ==============================
def assign_allene_ra_sa(mol, axis):
    left_idx, center_idx, right_idx = axis
    left = mol.GetAtomWithIdx(left_idx)
    right = mol.GetAtomWithIdx(right_idx)

    left_subs = [n for n in left.GetNeighbors() if n.GetIdx()!=center_idx]
    right_subs = [n for n in right.GetNeighbors() if n.GetIdx()!=center_idx]

    if len(left_subs)<2 or len(right_subs)<2:
        return "Undetermined"

    left_sorted = sorted(left_subs, key=lambda a: a.GetAtomicNum(), reverse=True)
    right_sorted = sorted(right_subs, key=lambda a: a.GetAtomicNum(), reverse=True)

    molH = Chem.AddHs(mol)
    AllChem.EmbedMolecule(molH, AllChem.ETKDG())
    conf = molH.GetConformer()

    def vec(a,b):
        pa = np.array(conf.GetAtomPosition(a))
        pb = np.array(conf.GetAtomPosition(b))
        return pb-pa

    axis_vec = vec(left_idx, right_idx)
    left_vec = vec(left_idx, left_sorted[0].GetIdx())
    right_vec = vec(right_idx, right_sorted[0].GetIdx())

    cross = np.cross(left_vec, axis_vec)
    dot = np.dot(cross, right_vec)

    if dot > 0:
        return "Ra"
    else:
        return "Sa"


# ==============================
# Page config
# ==============================
st.set_page_config(page_title="Advanced Chemical Isomer Analysis", layout="wide")

st.markdown("""
<style>
    .stApp { background-color: white; color: black; }
</style>
<h2 style='color: #800000;'>Chemical Isomer Analysis System 2.0</h2>
""", unsafe_allow_html=True)


# ==============================
# 3D viewer
# ==============================
def render_3d(mol, title):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mblock = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=400, height=300)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    view.zoomTo()
    st.write(f"**{title}**")
    showmol(view, height=300, width=400)


# ==============================
# UI
# ==============================
compound_name = st.text_input("Enter Structure Name:", "")

if st.button("Analyze & Visualize Isomers"):

    if not compound_name:
        st.warning("Please enter a compound name first.")

    else:
        try:
            results = pcp.get_compounds(compound_name, 'name')

            if not results:
                st.error(f"No compound found for: {compound_name}")

            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)

                # ======================
                # Allene analysis
                # ======================
                st.subheader("Allene Axial Chirality Analysis")

                axes = detect_allene_axes(mol)

                if axes:
                    st.success(f"Detected {len(axes)} chiral allene axis/axes")

                    for ax in axes:
                        config = assign_allene_ra_sa(mol, ax)
                        st.write(f"Allene axis {ax} → configuration: {config}")

                else:
                    st.info("No chiral allene axis detected")

                # ======================
                # Normal stereoisomer pipeline
                # ======================
                mol_no_stereo = Chem.Mol(mol)

                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)

                for atom in mol_no_stereo.GetAtoms():
                    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

                isomers = list(EnumerateStereoisomers(mol_no_stereo))

                # 2D
                st.subheader("2D Structures")

                labels = []
                for i, iso in enumerate(isomers):
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    centers = Chem.FindMolChiralCenters(iso)
                    label = f"Isomer {i+1}: {centers}" if centers else f"Isomer {i+1}"
                    labels.append(label)

                if isomers:
                    img = Draw.MolsToGridImage(isomers, molsPerRow=3, subImgSize=(300, 300), legends=labels)
                    st.image(img, use_container_width=True)

                # 3D
                st.subheader("3D Structures")
                cols = st.columns(3)

                for i, iso in enumerate(isomers):
                    with cols[i % 3]:
                        render_3d(iso, labels[i])

        except Exception as e:
            st.error(f"Error: {e}")
