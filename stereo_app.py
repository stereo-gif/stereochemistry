import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, EnumerateStereoisomers
import py3Dmol

# 1. إعدادات الصفحة
st.set_page_config(page_title="TriStereo - Chemical Isomer Analysis", layout="wide")

# 2. تصميم الواجهة
st.markdown("""
<style>
    .main-title { color: #800000; font-family: 'serif'; border-bottom: 2px solid #dcdde1; text-align: center; }
    .stButton>button { background-color: #800000; color: white; font-weight: bold; }
</style>
""", unsafe_allow_html=True)

st.markdown("<h1 class='main-title'>TriStereo: Chemical Isomer Analysis</h1>", unsafe_allow_html=True)

# 3. البحث
compound_name = st.text_input("🔍 Enter Compound Name:")

if compound_name:
    try:
        results = pcp.get_compounds(compound_name, 'name')
        if results:
            c = results[0]
            smiles = c.isomeric_smiles
            mol = Chem.MolFromSmiles(smiles)
            
            # --- تحليل الأيزومرات (Cis/Trans & R/S) ---
            Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            
            # كشف الروابط الثنائية (E/Z - Cis/Trans)
            has_double_bond_stereo = False
            for bond in mol.getBonds():
                if bond.getStereo() != Chem.BondStereo.STEREONONE:
                    has_double_bond_stereo = True
                    break

            # عرض النتائج
            m1, m2, m3 = st.columns(3)
            m1.metric("Formula", c.molecular_formula)
            
            # تحديد الحالة بناءً على التحليل
            if chiral_centers and has_double_bond_stereo:
                status = "R/S & E/Z Active"
            elif chiral_centers:
                status = "Chiral (R/S)"
            elif has_double_bond_stereo:
                status = "Geometric (E/Z)"
            else:
                status = "Achiral"
            
            m2.metric("Stereo Type", status)
            m3.metric("Chiral Centers", len(chiral_centers))

            # تنبيه المستخدم لنوع الأيزومر
            if has_double_bond_stereo:
                st.success("✅ This molecule has Geometric Isomers (Cis/Trans or E/Z)!")
            if chiral_centers:
                st.info(f"🧬 This molecule has Optical Isomers (R/S) at centers: {chiral_centers}")

            st.divider()

            tab1, tab2 = st.tabs(["🌐 Visuals", "📏 Projections Info"])

            with tab1:
                col1, col2 = st.columns(2)
                with col1:
                    st.image(Draw.MolToImage(mol, size=(400, 400)), caption="2D Structure")
                with col2:
                    # عرض الـ 3D
                    mol_3d = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol_3d)
                    mblock = Chem.MolToMolBlock(mol_3d)
                    view = py3Dmol.view(width=400, height=400)
                    view.addModel(mblock, 'mol')
                    view.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
                    view.zoomTo()
                    st.components.v1.html(view._make_html(), height=400)

            with tab2:
                st.subheader("Projection Analysis")
                if st.button("Analyze for Fischer/Newman"):
                    st.write("Checking if the molecule can be represented...")
                    if len(chiral_centers) > 0:
                        st.write("- **Fischer:** Possible (Chiral centers detected).")
                    if mol.getBondBetweenAtoms(0, 1): # مثال بسيط
                        st.write("- **Newman:** Look down C1-C2 bond.")
                    st.warning("Note: Dynamic drawing for Fischer/Newman is under development. Use 3D view to visualize angles.")

        else:
            st.error("Compound not found.")
    except Exception as e:
        st.error(f"Error: {e}")
