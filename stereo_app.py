import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol

# 1. إعدادات الصفحة
st.set_page_config(page_title="TriStereo - Chemical Isomer Analysis", layout="wide")

# 2. تصميم الواجهة (CSS) لتنسيق الألوان والتخطيط
st.markdown("""
<style>
    .main-title { color: #800000; font-family: 'serif'; border-bottom: 2px solid #dcdde1; text-align: center; padding-bottom: 10px; }
    .guide-box { background-color: #f9f9f9; padding: 15px; border-left: 5px solid #800000; border-radius: 5px; margin-bottom: 20px; }
    .stButton>button { width: 100%; border-radius: 10px; background-color: #800000; color: white; height: 3em; font-weight: bold; }
    .stTabs [data-baseweb="tab-list"] { gap: 24px; }
    .stTabs [data-baseweb="tab"] { height: 50px; background-color: #f0f2f6; border-radius: 5px 5px 0px 0px; padding: 10px 20px; }
    .stTabs [aria-selected="true"] { background-color: #800000 !important; color: white !important; }
</style>
""", unsafe_allow_html=True)

st.markdown("<h1 class='main-title'>TriStereo: Chemical Isomer Analysis</h1>", unsafe_allow_html=True)

# 3. الدليل المرجعي (ملاحظاتك الشخصية)
with st.expander("📝 Stereoisomerism Quick Reference Guide"):
    st.markdown("""
    <div class="guide-box">
    1. <b>Cis / Trans (Relative):</b> Identical groups on same/opposite sides.<br>
    2. <b>E / Z (Absolute - CIP System):</b> High-priority groups (atomic number) are together (Z) or opposite (E).<br>
    3. <b>R / S (Optical):</b> Absolute configuration of chiral centers.<br>
    <i>*Note: E/Z is required when all 4 groups on the double bond are different.*</i>
    </div>
    """, unsafe_allow_html=True)

# 4. البحث عن المركب
compound_name = st.text_input("🔍 Enter Compound Name (e.g., 2-butanol, Lactic acid, Maleic acid):")

def render_3d(smiles):
    """دالة لتوليد عرض 3D تفاعلي باستخدام py3Dmol"""
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
        return "<b>3D model generation failed for this molecule.</b>"

if compound_name:
    try:
        # جلب البيانات من PubChem
        results = pcp.get_compounds(compound_name, 'name')
        if results:
            c = results[0]
            smiles = c.isomeric_smiles
            mol = Chem.MolFromSmiles(smiles)
            
            # فحص الكيرالية وتحديد المراكز
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            status = "Chiral 🧬" if chiral_centers else "Achiral ✨"
            
            # عرض المقاييس الأساسية
            m1, m2, m3 = st.columns(3)
            m1.metric("Formula", c.molecular_formula)
            m2.metric("Stereo Status", status)
            m3.metric("Chiral Centers", len(chiral_centers))

            st.divider()

            # --- التبويبات لمنع التشتت ---
            tab1, tab2, tab3 = st.tabs(["🖼️ Main Visuals (2D/3D)", "📏 Advanced Projections", "📋 Molecular Data"])

            with tab1:
                col_left, col_right = st.columns(2)
                with col_left:
                    st.subheader("2D Chemical Structure")
                    img = Draw.MolToImage(mol, size=(400, 400))
                    st.image(img, use_container_width=True)
                with col_right:
                    st.subheader("3D Interactive Model")
                    html_3d = render_3d(smiles)
                    st.components.v1.html(html_3d, height=420)

            with tab2:
                st.subheader("Structural Projections")
                st.write("Click a button to generate the specific projection for exam-style study:")
                p1, p2, p3 = st.columns(3)
                
                with p1:
                    if st.button("Generate Fischer"):
                        st.info("Displaying Fischer Projection")
                        # ملاحظة: استبدلي الرابط بالكود الفعلي للرسم لاحقاً
                        st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/6/6b/Fischer_projection_example.svg/250px-Fischer_projection_example.svg.png", caption="Fischer Representation")
                        

                with p2:
                    if st.button("Generate Newman"):
                        st.info("Displaying Newman Projection")
                        st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/2/23/Staggered_conformation.svg/250px-Staggered_conformation.svg.png", caption="Newman Representation")
                        

[Image of Newman projection of ethane]


                with p3:
                    if st.button("Generate Sawhorse"):
                        st.info("Displaying Sawhorse Projection")
                        st.write("*(Perspective side-view of the C-C bond)*")
                        

            with tab3:
                st.subheader("Full Metadata")
                st.json(c.to_dict())

        else:
            st.error("Compound not found in PubChem database. Please try IUPAC name.")
    except Exception as e:
        st.error(f"Something went wrong: {e}")
else:
    st.info("💡 Tip: Enter 'Glucose' or 'Alanine' to see complex stereocenters.")
