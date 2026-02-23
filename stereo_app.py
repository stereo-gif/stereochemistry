import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol

# 1. إعدادات الصفحة
st.set_page_config(page_title="TriStereo - Chemical Isomer Analysis", layout="wide")

# 2. تصميم الواجهة (CSS)
st.markdown("""
<style>
    .main-title { color: #800000; font-family: 'serif'; border-bottom: 2px solid #dcdde1; text-align: center; padding-bottom: 10px; }
    .guide-box { background-color: #f9f9f9; padding: 15px; border-left: 5px solid #800000; border-radius: 5px; margin-bottom: 20px; }
    .stButton>button { width: 100%; border-radius: 10px; background-color: #800000; color: white; height: 3em; font-weight: bold; }
    .stTabs [aria-selected="true"] { background-color: #800000 !important; color: white !important; }
</style>
""", unsafe_allow_html=True)

st.markdown("<h1 class='main-title'>TriStereo: Chemical Isomer Analysis</h1>", unsafe_allow_html=True)

# 3. الدليل المرجعي
with st.expander("📝 Stereoisomerism Quick Reference Guide"):
    st.markdown("""
    <div class="guide-box">
    1. <b>Cis / Trans:</b> Identical groups on same/opposite sides.<br>
    2. <b>E / Z:</b> Priority-based (CIP System). Together (Z) or opposite (E).<br>
    3. <b>R / S:</b> Absolute configuration of chiral centers.
    </div>
    """, unsafe_allow_html=True)

# 4. البحث
compound_name = st.text_input("🔍 Enter Compound Name:")

def render_3d(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mblock = Chem.MolToMolBlock(mol)
        view = py3Dmol.view(width=800, height=400)
        view.addModel(mblock, 'mol')
        view.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
        view.zoomTo()
        return view._make_html()
    except:
        return "<b>3D model generation failed.</b>"

if compound_name:
    try:
        results = pcp.get_compounds(compound_name, 'name')
        if results:
            c = results[0]
            smiles = c.isomeric_smiles
            mol = Chem.MolFromSmiles(smiles)
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            
            m1, m2, m3 = st.columns(3)
            m1.metric("Formula", c.molecular_formula)
            m2.metric("Status", "Chiral 🧬" if chiral_centers else "Achiral ✨")
            m3.metric("Chiral Centers", len(chiral_centers))

            st.divider()

            tab1, tab2, tab3 = st.tabs(["🖼️ Main Visuals", "📏 Projections", "📋 Molecular Data"])

            with tab1:
                col_left, col_right = st.columns(2)
                with col_left:
                    st.subheader("2D Structure")
                    st.image(Draw.MolToImage(mol, size=(400, 400)))
                with col_right:
                    st.subheader("3D Interactive Model")
                    st.components.v1.html(render_3d(smiles), height=420)

            with tab2:
                st.subheader("Structural Projections")
                p1, p2, p3 = st.columns(3)
                
                with p1:
                    if st.button("Generate Fischer"):
                        st.info("Fischer Projection View")
                        st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/6/6b/Fischer_projection_example.svg/250px-Fischer_projection_example.svg.png")

                with p2:
                    if st.button("Generate Newman"):
                        st.info("Newman Projection View")
                        st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/2/23/Staggered_conformation.svg/250px-Staggered_conformation.svg.png")

                with p3:
                    if st.button("Generate Sawhorse"):
                        st.info("Sawhorse Projection View")
                        st.write("Visualizing spatial arrangement from side-angle...")

            with tab3:
                st.json(c.to_dict())
        else:
            st.error("Compound not found.")
    except Exception as e:
        st.error(f"Error: {e}")
else:
    st.info("Enter a chemical name to begin.")
