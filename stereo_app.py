import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

# 1. إعدادات الصفحة
st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

# 2. واجهة المستخدم والتنسيق
st.markdown("""
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System</h2>
<div style="background-color: #ffffff; padding: 15px; border: 1px solid #e1e1e1; border-left: 5px solid #800000; margin-bottom: 20px; font-family: sans-serif;">
    <strong style="color: #800000;">Stereoisomerism Reference Guide:</strong><br>
    1. <b style="color: #b22222;">Cis / Trans:</b> Relative position of identical groups.<br>
    2. <b style="color: #b22222;">E / Z:</b> Based on CIP priority (Z = together, E = opposite).<br>
    3. <b style="color: #b22222;">R / S:</b> Absolute configuration of chiral centers.
</div>
""", unsafe_allow_html=True)

# 3. مدخلات المستخدم
compound_name = st.text_input("Enter Structure Name (e.g., But-2-ene):", "")

if st.button("Analyze Isomers"):
    if not compound_name:
        st.warning("Please enter a compound name first.")
    else:
        try:
            # البحث في قاعدة بيانات PubChem
            results = pcp.get_compounds(compound_name, 'name')
            
            if not results:
                st.error(f"❌ No compound found for: {compound_name}")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # إزالة الكيمياء الفراغية لتوليد كل الاحتمالات
                mol_no_stereo = Chem.Mol(mol)
                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
                
                # توليد الأيزومرات الفراغية
                isomers = list(EnumerateStereoisomers(mol_no_stereo))
                
                # التحقق إذا كان المركب Achiral
                # إذا وجدنا أيزومر واحد فقط، فهذا يعني أن المركب لا يملك تماثل فراغي متغير
                if len(isomers) <= 1:
                    st.info(f"✨ The compound **{compound_name}** is **Achiral** (It has no geometric or optical isomers).")
                
                labels = []
                for i, iso in enumerate(isomers):
                    # تعيين الكيمياء الفراغية (E/Z, R/S)
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    
                    stereo_info = []
                    # استخراج معلومات الروابط الثنائية E/Z
                    for bond in iso.GetBonds():
                        stereo = bond.GetStereo()
                        if stereo == Chem.BondStereo.STEREOE:
                            stereo_info.append("E")
                        elif stereo == Chem.BondStereo.STEREOZ:
                            stereo_info.append("Z")
                    
                    # استخراج معلومات المراكز الكيرالية R/S
                    chiral_centers = Chem.FindMolChiralCenters(iso)
                    for center in chiral_centers:
                        stereo_info.append(f"({center[1]})")
                    
                    # تسمية الصورة
                    label = f"Isomer {i+1}: " + (", ".join(stereo_info) if stereo_info else "Achiral Structure")
                    labels.append(label)

                if len(isomers) > 0:
                    st.success(f"Analyzed Structure: **{compound_name}** | Found **{len(isomers)}** forms.")
                    
                    # رسم الشبكة الصور
                    img = Draw.MolsToGridImage(
                        isomers, 
                        molsPerRow=3, 
                        subImgSize=(400, 400), 
                        legends=labels
                    )
                    st.image(img, use_container_width=True)

        except Exception as e:
            st.error(f"An error occurred during analysis: {e}")

# تذييل الصفحة
st.markdown("---")
st.caption("Chemistry Tool | Developed with RDKit & Streamlit")
