import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from stmol import showmol
import py3Dmol
import numpy as np


# ==============================
# Allene axial chirality detector (robust)
# ==============================
def detect_allene_axes(mol):
    axes = []
    patt = Chem.MolFromSmarts("C=C=C")
    matches = mol.GetSubstructMatches(patt)

    for match in matches:
        left_idx, center_idx, right_idx = match

        left = mol.GetAtomWithIdx(left_idx)
        right = mol.GetAtomWithIdx(right_idx)

        left_subs = [n for n in left.GetNeighbors() if n.GetIdx()!=center_idx]
        right_subs = [n for n in right.GetNeighbors() if n.GetIdx()!=center_idx]

        if len(left_subs)==2 and len(right_subs)==2:
            if left_subs[0].GetAtomicNum()!=left_subs[1].GetAtomicNum() and \
               right_subs[0].GetAtomicNum()!=right_subs[1].GetAtomicNum():
                axes.append((left_idx, center_idx, right_idx))

    return axes


# ==============================
# Ra / Sa assignment
# ==============================
def assign_allene_ra_sa(mol, axis):

    molH = Chem.AddHs(mol)
    AllChem.EmbedMolecule(molH, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(molH)

    left_idx, center_idx, right_idx = axis

    left = molH.GetAtomWithIdx(left_idx)
    right = molH.GetAtomWithIdx(right_idx)

    left_subs = [n for n in left.GetNeighbors() if n.GetIdx()!=center_idx]
    right_subs = [n for n in right.GetNeighbors() if n.GetIdx()!=center_idx]

    left_sorted = sorted(left_subs, key=lambda a: a.GetAtomicNum(), reverse=True)
    right_sorted = sorted(right_subs, key=lambda a: a.GetAtomicNum(), reverse=True)

    conf = molH.GetConformer()

    def vec(a,b):
        pa = np.array(conf.GetAtomPosition(a))
        pb = np.array(conf.GetAtomPosition(b))
        return pb-pa

    axis_vec = vec(left_idx, right_idx)
    left_vec = vec(left_idx, left_sorted[0].GetIdx())
    right_vec = vec(right_idx, right_sorted[0].GetIdx())

    cross = np.cross(axis_vec, left_vec)
    dot = np.dot(cross, right_vec)

    if dot > 0:
        return "Ra"
    else:
        return "Sa"


# ==============================
# 3D viewer (stable)
# ==============================
def render_3d(mol, title):

    mol = Chem.AddHs(mol)

    if mol.GetNumConformers()==0:
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if result != 0:
            st.warning(f"3D embedding failed for {title}")
            return
        AllChem.UFFOptimizeMolecule(mol)

    mblock = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=400, height=300)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {}})
    view.zoomTo()

    st.write(f"**{title}**")
    showmol(view, height=300, width=400)


# ==============================
# UI
# ==============================
st.set_page_config(page_title="Allene Stereochemistry Analyzer", layout="wide")

st.markdown("<h2 style='color:#800000'>Allene Stereochemistry Analyzer</h2>", unsafe_allow_html=True)

compound_name = st.text_input("Enter Structure Name:", "")

if st.button("Analyze & Visualize Isomers"):

    if not compound_name:
        st.warning("Please enter a compound name")

    else:
        try:
            results = pcp.get_compounds(compound_name, 'name')

            if not results:
                st.error("Compound not found")

            else:
                smiles = results[0].smiles
                mol = Chem.MolFromSmiles(smiles)

                # ======================
                # Allene analysis
                # ======================
                st.subheader("Allene Axial Chirality Analysis")

                axes = detect_allene_axes(mol)

                if axes:
                    st.success(f"{len(axes)} chiral allene axis detected")

                    for ax in axes:
                        config = assign_allene_ra_sa(mol, ax)
                        st.write(f"Axis {ax} → {config}")

                else:
                    st.info("No chiral allene detected")

                # ======================
                # Stereoisomers
                # ======================
                mol_no = Chem.Mol(mol)

                for b in mol_no.GetBonds():
                    b.SetStereo(Chem.BondStereo.STEREONONE)

                for a in mol_no.GetAtoms():
                    a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

                isomers = list(EnumerateStereoisomers(mol_no))

                # 2D
                st.subheader("2D Isomers")

                labels=[]
                for i, iso in enumerate(isomers):
                    centers = Chem.FindMolChiralCenters(iso)
                    labels.append(f"Isomer {i+1} {centers}")

                if isomers:
                    img = Draw.MolsToGridImage(isomers, molsPerRow=3, legends=labels)
                    st.image(img, use_container_width=True)

                # 3D
                st.subheader("3D Isomers")

                if isomers:
                    cols = st.columns(3)
                    for i, iso in enumerate(isomers):
                        with cols[i%3]:
                            render_3d(iso, labels[i])

        except Exception as e:
            st.error(str(e))
