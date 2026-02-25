import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol
import numpy as np

# ==============================
# 1. Allene Detection Functions
# ==============================
def detect_allene_axes(mol):
    """كشف محاور الألين في الجزيء"""
    axes = []
    for bond1 in mol.GetBonds():
        if bond1.GetBondType() != Chem.BondType.DOUBLE: continue
        a1 = bond1.GetBeginAtom()
        a2 = bond1.GetEndAtom()
        for bond2 in a2.GetBonds():
            if bond2.GetIdx() == bond1.GetIdx(): continue
            if bond2.GetBondType() != Chem.BondType.DOUBLE: continue
            a3 = bond2.GetOtherAtom(a2)
            if a1.GetSymbol()=="C" and a2.GetSymbol()=="C" and a3.GetSymbol()=="C":
                l_subs = [n.GetIdx() for n in a1.GetNeighbors() if n.GetIdx()!=a2.GetIdx()]
                r_subs = [n.GetIdx() for n in a3.GetNeighbors() if n.GetIdx()!=a2.GetIdx()]
                if len(l_subs)==2 and len(r_subs)==2:
                    axes.append((a1.GetIdx(), a2.GetIdx(), a3.GetIdx()))
    return axes

def assign_allene_ra_sa(mol, axis):
    """حساب التماثل المحوري Ra/Sa"""
    try:
        left_idx, center_idx, right_idx = axis
        # إضافة هيدروجين وتوليد شكل 3D ضروري للحساب
        mol_3d = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG()) == -1: return "Error"
        
        conf = mol_3d.GetConformer()
        left_atom = mol_3d.GetAtomWithIdx(left_idx)
        right_atom = mol_3d.GetAtomWithIdx(right_idx)
        
        l_subs = sorted([n for n in left_atom.GetNeighbors() if n.GetIdx()!=center_idx], 
                        key=lambda x: x.GetAtomicNum(), reverse=True)
        r_subs = sorted([n for n in right_atom.GetNeighbors() if n.GetIdx()!=center_idx], 
                        key=lambda x: x.GetAtomicNum(), reverse=True)
        
        def vec(a_idx, b_idx):
            return np.array(conf.GetAtomPosition(b_idx)) - np.array(conf.GetAtomPosition(a_idx))
        
        axis_v = vec(left_idx, right_idx)
        l_v = vec(left_idx, l_subs[0].GetIdx())
        r_v = vec(right_idx, r_subs[0].GetIdx())
        
        dot = np.dot(np.cross(l_v, axis_v), r_v)
        return "Ra" if dot > 0 else "Sa"
    except:
        return "N/A"

# ==============================
# 2. Visualization Logic
# ==============================
def render_3d_safe(mol, title):
    """دالة عرض الـ 3D مع التأكد من توليد الـ Conformers"""
    try:
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        mblock = Chem.MolToMolBlock(mol_3d)
        view = py3Dmol.view(width=400, height=300)
        view.addModel(mblock, 'mol')
        view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
        view.zoomTo()
        st.write(f"**{title}**")
        showmol(view, height=300, width=400)
    except:
        st.error(f"Could not render 3D for {title}")

# ==============================
# 3. Streamlit UI
# ==============================
st.set_page_config(page_title="StereoMaster Pro", layout="wide")
st.title("🧪 Advanced Stereo Analysis (R/S, E/Z, Ra/Sa)")

compound_name = st.text_input("Enter Compound Name (e.g. 2,3-pentadiene):")

if st.button("Run Analysis"):
    if compound_name:
        with st.spinner("Processing..."):
            comps = pcp.get_compounds(compound_name, 'name')
            if not comps:
                st.error("Not found in PubChem.")
            else:
                base_smiles = comps[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # توليد الـ Isomers
                mol_unspec = Chem.Mol(mol)
                for b in mol_unspec.GetBonds(): b.SetStereo(Chem.BondStereo.STEREONONE)
                for a in mol_unspec.GetAtoms(): a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
                
                isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol_unspec))
                
                st.subheader(f"Found {len(isomers)} Isomers")
                
                cols_2d = st.columns(len(isomers) if len(isomers) < 4 else 4)
                
                for i, iso in enumerate(isomers):
                    # حساب الـ R/S العادي
                    Chem.AssignStereochemistry(iso, force=True)
                    rs_centers = Chem.FindMolChiralCenters(iso)
                    
                    # حساب الـ Ra/Sa للألينات
                    allene_axes = detect_allene_axes(iso)
                    allene_labels = []
                    for ax in allene_axes:
                        allene_labels.append(assign_allene_ra_sa(iso, ax))
                    
                    # تجميع الـ Label النهائي
                    label = f"Isomer {i+1}\n"
                    if rs_centers: label += f"Centers: {rs_centers}\n"
                    if allene_labels: label += f"Allene: {allene_labels}"
                    
                    # عرض 2D
                    with st.container():
                        st.image(Draw.MolToImage(iso, legend=label, size=(300, 300)))
                        # عرض 3D تحت كل صورة
                        render_3d_safe(iso, f"3D Model {i+1}")
    else:
        st.warning("Please enter a name.")           
