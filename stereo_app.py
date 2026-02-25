import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, EnumerateStereoisomers
import py3Dmol
import numpy as np

# ==============================
# 1. الدالة السحرية لحساب Ra/Sa للألين
# ==============================
def get_allene_label(iso_mol):
    """بتحسب الـ Axial Chirality لكل أيزومر على حدة"""
    # لازم نعمل نسخة ونضيف هيدروجين ونثبت الشكل 3D
    m = Chem.AddHs(iso_mol)
    if AllChem.EmbedMolecule(m, AllChem.ETKDG()) == -1:
        return ""
    
    AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer()
    
    labels = []
    for bond in m.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # البحث عن المركز (C=C=C)
            for n_bond in a2.GetBonds():
                if n_bond.GetIdx() == bond.GetIdx(): continue
                if n_bond.GetBondType() == Chem.BondType.DOUBLE:
                    a3 = n_bond.GetOtherAtom(a2)
                    
                    # تحديد المجموعات ذات الأولوية
                    def get_high_priority(atom, exclude_idx):
                        subs = [n for n in atom.GetNeighbors() if n.GetIdx() != exclude_idx]
                        if not subs: return None
                        return sorted(subs, key=lambda x: x.GetAtomicNum(), reverse=True)[0]
                    
                    l_sub = get_high_priority(a1, a2.GetIdx())
                    r_sub = get_high_priority(a3, a2.GetIdx())
                    
                    if l_sub and r_sub:
                        # حساب المتجهات (Vector Math)
                        p1 = np.array(conf.GetAtomPosition(a1.GetIdx()))
                        p3 = np.array(conf.GetAtomPosition(a3.GetIdx()))
                        pl = np.array(conf.GetAtomPosition(l_sub.GetIdx()))
                        pr = np.array(conf.GetAtomPosition(r_sub.GetIdx()))
                        
                        v_axis = p3 - p1
                        v_l = pl - p1
                        v_r = pr - p3
                        
                        # ضرب اتجاهي لتحديد الالتفاف
                        dot = np.dot(np.cross(v_l, v_axis), v_r)
                        labels.append("Ra" if dot > 0 else "Sa")
    
    return " | ".join(labels) if labels else ""

# ==============================
# 2. عرض الـ 3D بطريقة مضمونة
# ==============================
def show_3d_render(mol):
    """بتعرض الـ 3D باستخدام py3Dmol داخل Streamlit"""
    m = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m, AllChem.ETKDG())
    m_block = Chem.MolToMolBlock(m)
    
    # بناء الـ HTML الخاص بـ py3Dmol
    xyzview = py3Dmol.view(width=400, height=300)
    xyzview.addModel(m_block, 'mol')
    xyzview.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
    xyzview.zoomTo()
    
    # تحويل العرض لـ HTML عشان يظهر في Streamlit
    obj = xyzview._make_html()
    st.components.v1.html(obj, height=350)

# ==============================
# 3. واجهة البرنامج
# ==============================
st.set_page_config(page_title="Chemical Isomer Pro 2026", layout="wide")
st.markdown("<h2 style='text-align: center; color: #4A90E2;'>Stereo Analysis System</h2>", unsafe_allow_html=True)

compound_name = st.text_input("Structure Name:", "2,3-pentadiene")

if st.button("Start Analysis"):
    try:
        results = pcp.get_compounds(compound_name, 'name')
        if not results:
            st.error("Compound not found.")
        else:
            base_mol = Chem.MolFromSmiles(results[0].smiles)
            
            # 1. تصفير أي استيريو موجود
            m_unspec = Chem.Mol(base_mol)
            for b in m_unspec.GetBonds(): b.SetStereo(Chem.BondStereo.STEREONONE)
            for a in m_unspec.GetAtoms(): a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            
            # 2. توليد الأيزومرات
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(m_unspec))
            
            st.info(f"Analyzed {compound_name}: Found {len(isomers)} potential isomers.")
            
            # 3. العرض
            cols = st.columns(2)
            for i, iso in enumerate(isomers):
                with cols[i % 2]:
                    st.subheader(f"Isomer #{i+1}")
                    
                    # حساب الـ R/S العادي
                    Chem.AssignStereochemistry(iso, force=True)
                    rs = Chem.FindMolChiralCenters(iso)
                    
                    # حساب الـ Ra/Sa للألين (هنا التعديل المهم)
                    axial = get_allene_label(iso)
                    
                    # عرض البيانات (إيه اللي موجود في الجزيء ده؟)
                    if rs: st.write(f"**Chiral Centers (R/S):** `{rs}`")
                    if axial: st.success(f"**Axial Chirality (Ra/Sa):** `{axial}`")
                    if not rs and not axial: st.write("No specific stereocenters detected.")
                    
                    # الـ 2D
                    st.image(Draw.MolToImage(iso, size=(400, 400)))
                    
                    # الـ 3D (الآن شغال بإذن الله)
                    show_3d_render(iso)
                    st.divider()

    except Exception as e:
        st.error(f"Something went wrong: {e}")
