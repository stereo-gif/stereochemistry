import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

# 1. إعدادات الصفحة
st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

st.markdown("""
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System</h2>
<div style="background-color: #ffffff; padding: 15px; border: 1px solid #e1e1e1; border-left: 4px solid #800000; margin-bottom: 20px;">
	<strong style="color: #800000;">Stereoisomerism Reference Guide:</strong><br>
	1. <b>Cis / Trans:</b> Relative position.<br>
	2. <b>E / Z:</b> CIP System for double bonds.<br>
	3. <b>R / S:</b> Absolute configuration (including Axial Chirality for Allenes).
</div>
""", unsafe_allow_html=True)

compound_name = st.text_input("Enter Structure Name (e.g., 1,3-dichloropropadiene):", "")

if st.button("Analyze Isomers"):
    if not compound_name:
        st.warning("Please enter a name.")
    else:
        try:
            # جلب المركب من PubChem
            results = pcp.get_compounds(compound_name, 'name')
            if not results:
                st.error("Compound not found.")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # --- التعديل الجذري لحل مشكلة الألّين ---
                # 1. إضافة هيدروجينات (مهم جداً للـ Allenes عشان الزوايا تبان)
                mol = Chem.AddHs(mol)
                
                # 2. مسح أي معلومات فراغية قديمة
                Chem.RemoveStereochemistry(mol)
                
                # 3. تفعيل البحث عن مراكز الكايرالية المحورية والروابط
                # نستخدم الـ Flag الخاص بالـ Allenes
                Chem.FindPotentialStereoBonds(mol)
                
                # 4. إعداد خيارات التوليد مع إجبار البرنامج على فحص الـ 3D
                options = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol, options=options))
                
                # لو لسه مش شايف أيزومرز، بنجرب طريقة "التلاعب بالروابط"
                if len(isomers) == 1:
                    # محاولة يدوية لتحديد روابط الألّين كـ Bonds قابلة للتغير
                    for bond in mol.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            bond.SetStereo(Chem.BondStereo.STEREONONE)
                    isomers = list(EnumerateStereoisomers(mol, options=options))

                labels = []
                final_mols = []
                for i, iso in enumerate(isomers):
                    # إزالة الهيدروجينات للرسم النظيف بعد الحساب
                    clean_iso = Chem.RemoveHs(iso)
                    Chem.AssignStereochemistry(clean_iso, force=True, cleanIt=True)
                    
                    info = []
                    # البحث عن R/S أو Axial Chirality
                    centers = Chem.FindMolChiralCenters(clean_iso, includeUnassigned=True)
                    for c in centers:
                        info.append(f"({c[1]})")
                    
                    # البحث عن E/Z
                    for b in clean_iso.GetBonds():
                        if b.GetStereo() == Chem.BondStereo.STEREOE: info.append("E")
                        elif b.GetStereo() == Chem.BondStereo.STEREOZ: info.append("Z")
                    
                    label = f"Isomer {i+1}: " + (", ".join(set(info)) if info else "Achiral")
                    labels.append(label)
                    final_mols.append(clean_iso)

                st.success(f"Found **{len(isomers)}** forms for {compound_name}")
                
                img = Draw.MolsToGridImage(
                    final_mols, 
                    molsPerRow=2, 
                    subImgSize=(400, 400), 
                    legends=labels,
                    useSVG=True # SVG بيخلي الخطوط أوضح في الـ Allenes
                )
                st.write(img, unsafe_allow_html=True)

        except Exception as e:
            st.error(f"Error: {e}")
