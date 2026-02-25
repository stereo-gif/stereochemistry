import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, EnumerateStereoisomers
from stmol import showmol
import py3Dmol

# ==============================
# دالة عرض الـ 3D المستقرة جداً
# ==============================
def render_3d(mol):
    # إضافة الهيدروجين ضروري جداً لشكل الـ 3D
    mol_3d = Chem.AddHs(mol)
    # توليد الإحداثيات (الخطوة اللي بتخلي الـ 3D يظهر)
    AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
    # تحويل الجزيء لـ Block بيفهمه py3Dmol
    mblock = Chem.MolToMolBlock(mol_3d)
    
    view = py3Dmol.view(width=400, height=300)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    view.zoomTo()
    showmol(view, height=300, width=400)

# ==============================
# واجهة البرنامج
# ==============================
st.set_page_config(layout="wide")
st.title("Chemical Isomer Analyzer (Focus: 3D & R/S)")

compound_name = st.text_input("Enter Compound Name (e.g., Ibuprofen, Thalidomide):", "Ibuprofen")

if st.button("Generate Isomers"):
    try:
        # البحث في PubChem
        results = pcp.get_compounds(compound_name, 'name')
        
        if not results:
            st.error("Compound not found!")
        else:
            # تحويل الـ SMILES لجزيء
            base_smiles = results[0].smiles
            mol = Chem.MolFromSmiles(base_smiles)
            
            # توليد جميع الأيزومرات الفراغية الممكنة (R/S)
            opts = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=True)
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
            
            st.success(f"Found {len(isomers)} potential isomers for {compound_name}")
            
            # عرض النتائج في أعمدة
            cols = st.columns(2)
            for i, iso in enumerate(isomers):
                with cols[i % 2]:
                    # إجبار الـ RDKit على كشف الـ R/S
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                    
                    st.subheader(f"Isomer {i+1}")
                    st.write(f"**Chiral Centers (R/S):** `{centers}`")
                    
                    # عرض الـ 2D أولاً كمرجع
                    st.image(Draw.MolToImage(iso, size=(400, 400)))
                    
                    # عرض الـ 3D (الأساس)
                    st.write("**Interactive 3D Structure:**")
                    render_3d(iso)
                    st.divider()
                    
    except Exception as e:
        st.error(f"Error: {e}")
