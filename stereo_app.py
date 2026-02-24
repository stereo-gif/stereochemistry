import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from stmol import showmol
import py3Dmol

# 1. إعدادات الصفحة
st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

# 2. تصميم الواجهة (مع دليل المراجعة الخاص بكِ)
st.markdown("""
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System</h2>
<div style="background-color: #ffffff; padding: 15px; border: 1px solid #e1e1e1; border-left: 4px solid #800000; margin-bottom: 20px; font-family: sans-serif;">
	<strong style="color: #800000;">Stereoisomerism Reference Guide:</strong><br>
	1. <b style="color: #b22222;">Cis / Trans:</b> Relative position.<br>
	2. <b style="color: #b22222;">E / Z:</b> CIP System.<br>
	3. <b style="color: #b22222;">R / S:</b> Absolute configuration (including Axial Chirality for Allenes).
</div>
""", unsafe_allow_html=True)

# 3. مدخلات المستخدم
compound_name = st.text_input("Enter Structure Name (e.g., 1,3-dichloropropadiene):", "")

if st.button("Analyze Isomers"):
    if not compound_name:
        st.warning("Please enter a name.")
    else:
        try:
            results = pcp.get_compounds(compound_name, 'name')
            if not results:
                st.error("❌ Compound not found.")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # --- [تجهيز الألّين] ---
                # البحث عن الكربون المركزي في الألّين وتنشيطه
                patt = Chem.MolFromSmarts("C=C=C")
                if mol.HasSubstructMatch(patt):
                    matches = mol.GetSubstructMatches(patt)
                    for match in matches:
                        mol.GetAtomWithIdx(match[1]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

                mol = Chem.AddHs(mol)
                opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol, options=opts))

                st.success(f"Found **{len(isomers)}** Isomers for {compound_name}")

                # عرض الأيزومرز
                cols = st.columns(len(isomers))
                for idx, iso in enumerate(isomers):
                    with cols[idx]:
                        # حساب الكيمياء الفراغية
                        Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                        centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                        label = f"Isomer {idx+1}: " + (f"({centers[0][1]})" if centers else "Axial Chirality")
                        st.subheader(label)

                        # --- رسم الـ 3D ---
                        # توليد إحداثيات 3D للألّين
                        AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                        mblock = Chem.MolToMolBlock(iso)
                        
                        view = py3Dmol.view(width=400, height=400)
                        view.addModel(mblock, 'mol')
                        view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
                        view.zoomTo()
                        showmol(view, height=400)

                        # رسم الـ 2D تحتها للتأكيد
                        img = Draw.MolToImage(Chem.RemoveHs(iso), size=(300, 300))
                        st.image(img, caption=f"2D Projection of {label}")

        except Exception as e:
            st.error(f"Error: {e}")

st.markdown("---")
st.caption("Powered by RDKit & py3Dmol")
