import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

# 1. إعدادات الصفحة
st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

# 2. تصميم الواجهة (مع دليل المراجعة الخاص بكِ)
st.markdown("""
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System</h2>
<div style="background-color: #ffffff; padding: 15px; border: 1px solid #e1e1e1; border-left: 4px solid #800000; margin-bottom: 20px; font-family: sans-serif;">
	<strong style="color: #800000;">Stereoisomerism Reference Guide:</strong><br>
	1. <b style="color: #b22222;">Cis / Trans (Relative):</b> Identical groups on same/opposite sides.<br>
	2. <b style="color: #b22222;">E / Z (Absolute - CIP System):</b> <b>Z (Zusammen)</b> together, <b>E (Entgegen)</b> opposite.<br>
	3. <b style="color: #b22222;">R / S (Optical):</b> Absolute configuration of chiral centers (and axial chirality like Allenes).
</div>
""", unsafe_allow_html=True)

# 3. مدخلات المستخدم
compound_name = st.text_input("Enter Structure Name (e.g., 1,3-dichloropropadiene or But-2-ene):", "")

if st.button("Analyze Isomers"):
    if not compound_name:
        st.warning("Please enter a compound name first.")
    else:
        try:
            # جلب البيانات من PubChem
            results = pcp.get_compounds(compound_name, 'name')
            
            if not results:
                st.error(f"❌ No compound found for: {compound_name}")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # إزالة أي معلومات فراغية مسبقة لإعادة حساب كل الاحتمالات (بما فيها الـ Allenes)
                mol_no_stereo = Chem.Mol(mol)
                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
                for atom in mol_no_stereo.GetAtoms():
                    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

                # تفعيل البحث عن مراكز الـ Stereo المحتملة (الروابط والمحاور)
                Chem.AssignStereochemistry(mol_no_stereo, force=True, cleanIt=True)
                
                # خيارات الـ Enumeration: 
                # tryEmbedding=True ضرورية جداً للـ Allenes لأنها تعتمد على الشكل ثلاثي الأبعاد
                options = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol_no_stereo, options=options))

                if len(isomers) <= 1:
                    st.info(f"✨ The compound {compound_name} appears to have no additional geometric or optical isomers in this configuration.")

                labels = []
                for i, iso in enumerate(isomers):
                    # حساب الـ Stereo لكل أيزومر ناتج
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    stereo_info = []
                    
                    # 1. كشف الـ E/Z للروابط الثنائية
                    for bond in iso.GetBonds():
                        st_bond = bond.GetStereo()
                        if st_bond == Chem.BondStereo.STEREOE: stereo_info.append("E")
                        elif st_bond == Chem.BondStereo.STEREOZ: stereo_info.append("Z")
                    
                    # 2. كشف الـ R/S للمراكز الكايرالية (والـ Axial Chirality في بعض الحالات)
                    chiral_centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                    for center in chiral_centers:
                        stereo_info.append(f"({center[1]})")
                    
                    label = f"Isomer {i+1}: " + (", ".join(set(stereo_info)) if stereo_info else "Achiral Structure")
                    labels.append(label)

                st.success(f"Analyzed Structure: {compound_name} | Found Forms: **{len(isomers)}**")
                
                # عرض النتائج في شبكة صور
                img = Draw.MolsToGridImage(
                    isomers,
                    molsPerRow=2 if len(isomers) <= 2 else 3,
                    subImgSize=(400, 400),
                    legends=labels
                )
                
                st.image(img, use_container_width=True)

        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")

st.markdown("---")
st.caption("Developed for Stereochemistry Analysis | Powered by RDKit & PubChemPy")
