           import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol
import numpy as np

# ==============================
# 1. Allene Detection Logic
# ==============================
def get_allene_config(mol):
    """يكشف محاور الألين ويحدد Ra/Sa لكل محور"""
    axes_configs = []
    # لازم يكون عندنا 3D عشان نحسب الـ Ra/Sa
    mol_3d = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG()) == -1:
        return []

    for bond in mol_3d.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # نبحث عن ذرة كربون في النص بين رابطتين مزدوجتين
            for next_bond in a2.GetBonds():
                if next_bond.GetIdx() == bond.GetIdx(): continue
                if next_bond.GetBondType() == Chem.BondType.DOUBLE:
                    a3 = next_bond.GetOtherAtom(a2)
                    
                    # كشف الأطراف (Substituents)
                    left_subs = [n for n in a1.GetNeighbors() if n.GetIdx() != a2.GetIdx()]
                    right_subs = [n for n in a3.GetNeighbors() if n.GetIdx() != a2.GetIdx()]
                    
                    if len(left_subs) >= 2 and len(right_subs) >= 2:
                        # ترتيب حسب الأولوية (Atomic Number)
                        l_high = sorted(left_subs, key=lambda x: x.GetAtomicNum(), reverse=True)[0]
                        r_high = sorted(right_subs, key=lambda x: x.GetAtomicNum(), reverse=True)[0]
                        
                        conf = mol_3d.GetConformer()
                        def get_p(atom): return np.array(conf.GetAtomPosition(atom.GetIdx()))
                        
                        # المتجهات
                        v_axis = get_p(a3) - get_p(a1)
                        v_l = get_p(l_high) - get_p(a1)
                        v_r = get_p(r_high) - get_p(a3)
                        
                        dot = np.dot(np.cross(v_l, v_axis), v_r)
                        config = "Ra" if dot > 0 else "Sa"
                        axes_configs.append(f"Axis({a1.GetIdx()}-{a3.GetIdx()}): {config}")
    return axes_configs

# ==============================
# 2. 3D Renderer Fix
# ==============================
def make_3d_viewer(mol, width=400, height=300):
    """يضمن ظهور الـ 3D عبر تحويله لـ Block سليم"""
    m = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(m)
    m_block = Chem.MolToMolBlock(m)
    
    view = py3Dmol.view(width=width, height=height)
    view.addModel(m_block, 'mol')
    view.setStyle({'stick': {'colorscheme': 'Jmol', 'radius': 0.15}, 'sphere': {'scale': 0.25}})
    view.zoomTo()
    return view

# ==============================
# 3. Streamlit Interface
# ==============================
st.set_page_config(page_title="Stereo-Explorer 2026", layout="wide")
st.title("🔬 Stereo-Isomer Professional Analyzer")

name = st.text_input("Structure Name:", value="2,3-pentadiene")

if st.button("Analyze Now"):
    try:
        results = pcp.get_compounds(name, 'name')
        if not results:
            st.error("Compound not found.")
        else:
            smiles = results[0].smiles
            main_mol = Chem.MolFromSmiles(smiles)
            
            # توليد جميع الأيزومرات
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(main_mol))
            
            st.success(f"Successfully generated {len(isomers)} isomers.")
            
            # العرض في أعمدة
            cols = st.columns(2) # خليناها 2 عشان الـ 3D ياخد مساحته
            
            for i, iso in enumerate(isomers):
                with cols[i % 2]:
                    st.markdown(f"### Isomer {i+1}")
                    
                    # حساب الخصائص
                    Chem.AssignStereochemistry(iso, force=True)
                    rs_centers = Chem.FindMolChiralCenters(iso)
                    allene_data = get_allene_config(iso)
                    
                    # كتابة البيانات
                    st.write(f"**R/S Centers:** {rs_centers if rs_centers else 'None'}")
                    if allene_data:
                        st.info(f"**Axial Chirality:** {allene_data}")
                    
                    # عرض 2D
                    st.image(Draw.MolToImage(iso, size=(400, 400)))
                    
                    # عرض 3D
                    st.write("**Interactive 3D View:**")
                    view = make_3d_viewer(iso)
                    showmol(view, height=300, width=400)
                    st.divider()
                    
    except Exception as e:
        st.error(f"Error: {e}")
