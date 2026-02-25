import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol
import numpy as np

# دالة لحساب الـ Ra/Sa للألين
def get_allene_label(mol):
    m = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(m, AllChem.ETKDG()) == -1: return ""
    conf = m.GetConformer()
    results = []
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE:
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            for nb in a2.GetBonds():
                if nb.GetIdx() == b.GetIdx(): continue
                if nb.GetBondType() == Chem.BondType.DOUBLE:
                    a3 = nb.GetOtherAtom(a2)
                    l_subs = sorted([n for n in a1.GetNeighbors() if n.GetIdx()!=a2.GetIdx()], key=lambda x: x.GetAtomicNum(), reverse=True)
                    r_subs = sorted([n for n in a3.GetNeighbors() if n.GetIdx()!=a2.GetIdx()], key=lambda x: x.GetAtomicNum(), reverse=True)
                    if len(l_subs)>0 and len(r_subs)>0:
                        v_ax = np.array(conf.GetAtomPosition(a3.GetIdx())) - np.array(conf.GetAtomPosition(a1.GetIdx()))
                        v_l = np.array(conf.GetAtomPosition(l_subs[0].GetIdx())) - np.array(conf.GetAtomPosition(a1.GetIdx()))
                        v_r = np.array(conf.GetAtomPosition(r_subs[0].GetIdx())) - np.array(conf.GetAtomPosition(a3.GetIdx()))
                        dot = np.dot(np.cross(v_l, v_ax), v_r)
                        results.append("Ra" if dot > 0 else "Sa")
    return " | ".join(results)

# دالة عرض الـ 3D
def render_3d(mol):
    m = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m, AllChem.ETKDG())
    m_block = Chem.MolToMolBlock(m)
    view = py3Dmol.view(width=400, height=300)
    view.addModel(m_block, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    view.zoomTo()
    showmol(view, height=300, width=400)

st.title("Advanced Chemical Isomer Analyzer")

name = st.text_input("Enter Molecule Name:", "2,3-pentadiene")

if st.button("Analyze"):
    results = pcp.get_compounds(name, 'name')
    if results:
        smiles = results[0].smiles
        mol = Chem.MolFromSmiles(smiles)
        
        # أهم خطوة: السماح بفك التشفير لكل الاحتمالات
        opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True)
        isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
        
        cols = st.columns(2)
        for i, iso in enumerate(isomers):
            with cols[i % 2]:
                # إجبار RDKit على حساب الـ R/S
                Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                # نستخدم includeUnassigned=True للتأكد من فحص كل الذرات
                centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                axial = get_allene_label(iso)
                
                st.subheader(f"Isomer {i+1}")
                st.write(f"**Centers:** {centers}")
                if axial: st.success(f"**Axial:** {axial}")
                
                # عرض 2D
                st.image(Draw.MolToImage(iso, size=(300, 300)))
                
                # عرض 3D
                try:
                    render_3d(iso)
                except:
                    st.error("3D Render Failed")
