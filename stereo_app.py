import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

# 1. إعدادات الصفحة
st.set_page_config(
    page_title="Chemical Isomer Analysis", 
    layout="wide",
    initial_sidebar_state="collapsed"
)

# 2. تصميم الواجهة (Custom CSS لضمان ثبات الألوان في الـ Light Mode)
st.markdown("""
<style>
    .main {
        background-color: #ffffff;
    }
    .stTextInput > div > div > input {
        background-color: #f8f9fa;
    }
    .reportview-container .main .block-container {
        padding-top: 2rem;
    }
</style>

<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1; padding-bottom: 10px;'>
    Chemical Isomer Analysis System
</h2>

<div style="background-color: #fdf2f2; padding: 20px; border: 1px solid #ffcfcf; border-left: 5px solid #800000; border-radius: 5px; margin-bottom: 25px; font-family: sans-serif;">
    <strong style="color: #800000; font-size: 18px;">Stereoisomerism Reference Guide:</strong><br><br>
    <ul style="list-style-type: none; padding-left: 0;">
        <li>1. <b style="color: #b22222;">Cis / Trans (Relative):</b> Identical groups on the same or opposite sides of a double bond.</li>
        <li>2. <b style="color: #b22222;">E / Z (Absolute - CIP System):</b> <b>Z (Zusammen)</b> high-priority together, <b>E (Entgegen)</b> opposite.</li>
        <li>3. <b style="color: #b22222;">R / S (Optical):</b> Absolute configuration of chiral centers based on priority.</li>
    </ul>
    <small style="color: #666;">*Note: E/Z is required when all 4 groups on the double bond are different.</small>
</div>
""", unsafe_allow_html=True)

# 3. مدخلات المستخدم
compound_name = st.text_input("Enter Structure Name (e.g., 1,2-dichloroethene or Thalidomide):", "")

if st.button("Analyze Isomers"):
    if not compound_name:
        st.warning("Please enter a compound name first.")
    else:
        try:
            with st.spinner('Fetching data from PubChem and analyzing...'):
                results = pcp.get_compounds(compound_name, 'name')
            
            if not results:
                st.error(f"❌ No compound found for: {compound_name}")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # إزالة أي معلومات فراغية موجودة للبدء من الصفر في التوليد
                mol_no_stereo = Chem.Mol(mol)
                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
                
                # توليد الأيزومرات
                isomers = list(EnumerateStereoisomers(mol_no_stereo))
                
                # التحقق من الـ Achiral
                if len(isomers) <= 1:
                    st.info(f"✨ The compound **{compound_name}** is **Achiral**. It does not have geometric (E/Z) or optical (R/S) isomers.")

                labels = []
                for i, iso in enumerate(isomers):
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    stereo_info = []
                    
                    # فحص الروابط الثنائية (E/Z)
                    for bond in iso.GetBonds():
                        stereo = bond.GetStereo()
                        if stereo == Chem.BondStereo.STEREOE: stereo_info.append("E")
                        elif stereo == Chem.BondStereo.STEREOZ: stereo_info.append("Z")
                    
                    # فحص المراكز الكيرالية (R/S)
                    chiral_centers = Chem.FindMolChiralCenters(iso)
                    for center in chiral_centers:
                        stereo_info.append(f"({center[1]})")
                    
                    label = f"Isomer {i+1}: " + (", ".join(stereo_info) if stereo_info else "Achiral Structure")
                    labels.append(label)

                st.success(f"Analyzed Structure: **{compound_name}** | Total Forms: **{len(isomers)}**")
                
                # رسم الأيزومرات في شبكة
                img = Draw.MolsToGridImage(
                    isomers, 
                    molsPerRow=3, 
                    subImgSize=(350, 350), 
                    legends=labels
                )
                
                st.image(img, use_container_width=True)

        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")

st.markdown("---")
st.caption("Powered by RDKit, PubChemPy, and Streamlit.")
