import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol
import numpy as np

# ==============================
# 1. دالة الـ Ra/Sa (بسيطة ومستقرة)
# ==============================
def get_allene_config(mol):
    # بنعمل نسخة عشان منبوظش الجزيء الأصلي
    m = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(m, AllChem.ETKDG()) == -1: return ""
    conf = m.GetConformer()
    
    for bond in m.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            for nb in a2.GetBonds():
                if nb.GetIdx() == bond.GetIdx(): continue
                if nb.GetBondType() == Chem.BondType.DOUBLE:
                    a3 = nb.GetOtherAtom(a2)
                    l_subs = sorted([n for n in a1.GetNeighbors() if n.GetIdx() != a2.GetIdx()], key=lambda x: x.GetAtomicNum(), reverse=True)
                    r_subs = sorted([n for n in a3.GetNeighbors() if n.GetIdx() != a2.GetIdx()], key=lambda x: x.GetAtomicNum(), reverse=True)
                    if l_subs and r_subs:
                        p1, p3 = np.array(conf.GetAtomPosition(a1.GetIdx())), np.array(conf.GetAtomPosition(a3.GetIdx()))
                        pl, pr = np.array(conf.GetAtomPosition(l_subs[0].GetIdx())), np.array(conf.GetAtomPosition(r_subs[0].GetIdx()))
                        dot = np.dot(np.cross(pl-p1, p3-p1), pr-p3)
                        return "Ra" if dot > 0 else "Sa"
    return ""

# ==============================
# 2. دالة الـ 3D (اللي كانت شغالة معاكي)
# ==============================
def render_3d(mol):
    mol_3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
    mblock = Chem.MolToMolBlock(mol_3d)
    view = py3Dmol.view(width=400, height=300)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    view.zoomTo()
    showmol(view, height=300, width=400)

# ==============================
# 3. واجهة البرنامج
# ==============================
st.set_page_config(layout="wide")
st.title("Isomer Analyzer 2.0")

name = st.text_input("Enter Compound Name:", "2,3-pentadiene")

if st.button("Analyze"):
    results = pcp.get_compounds(name, 'name')
    if results:
        # تحويل الـ SMILES لجزيء وتوليد الأيزومرات
        base_mol = Chem.MolFromSmiles(results[0].smiles)
        isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(base_mol))
        
        st.write(f"Found {len(isomers)} Isomers")
        
        cols = st.columns(2)
        for i, iso in enumerate(isomers):
            with cols[i % 2]:
                # إجبار حساب الـ R/S
                Chem.AssignStereochemistry(iso, force=True)
                centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                axial = get_allene_config(iso)
                
                # العنوان والـ Labels
                label = f"Isomer {i+1} | R/S: {centers}"
                if axial: label += f" | Axial: {axial}"
                
                st.subheader(label)
                
                # الـ 2D
                st.image(Draw.MolToImage(iso, size=(300,300)))
                
                # الـ 3D
                render_3d(iso)
                st.divider()
