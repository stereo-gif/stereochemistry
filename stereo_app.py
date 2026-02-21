import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

# 1. إعدادات الصفحة
st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

st.markdown("""
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System</h2>
<div style="background-color: #ffffff; padding: 15px; border: 1px solid #e1e1e1; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000;">Stereoisomerism Reference Guide:</strong><br>
    1. <b>E / Z System:</b> High priority (1) vs Low priority (2) on each carbon.<br>
    2. <b>CIP Priority:</b> Higher atomic number gets priority (1).
</div>
""", unsafe_allow_html=True)

compound_name = st.text_input("Enter Structure Name:", "")

if st.button("Analyze Isomers"):
    if not compound_name:
        st.warning("Please enter a name.")
    else:
        try:
            results = pcp.get_compounds(compound_name, 'name')
            if not results:
                st.error(f"❌ No compound found.")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                mol_no_stereo = Chem.Mol(mol)
                
                # إزالة أي كيمياء فراغية لتوليد كل الاحتمالات
                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
                
                isomers = list(EnumerateStereoisomers(mol_no_stereo))
                
                if len(isomers) <= 1:
                    st.info(f"✨ The compound **{compound_name}** is **Achiral / No Geometric Isomers**.")

                labels = []
                for i, iso in enumerate(isomers):
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    
                    # حساب رتبة الأولوية (CIP Rank) لكل ذرة
                    # الرتبة الأعلى (أصغر رقم في الترتيب) تعني أولوية أعلى
                    ranks = list(Chem.CanvasV8.GetAtomPriorities(iso)) if hasattr(Chem, 'CanvasV8') else []
                    
                    # لو الرتب مش متاحة بالطريقة دي، بنستخدم الـ CIPRank المخفي
                    for atom in iso.GetAtoms():
                        if atom.HasProp('_CIPRank'):
                            # هنعرض الأولوية كـ (1) أو (2) لتسهيل الفهم
                            # لاحظي: RDKit بتدي رتبة صفر لأقل أولوية، فإحنا بنعكسها للتبسيط
                            rank_val = atom.GetProp('_CIPRank')
                            atom.SetProp('atomNote', f"p:{rank_val}")

                    stereo_info = []
                    for bond in iso.GetBonds():
                        stereo = bond.GetStereo()
                        if stereo == Chem.BondStereo.STEREOE: stereo_info.append("E")
                        elif stereo == Chem.BondStereo.STEREOZ: stereo_info.append("Z")
                    
                    centers = Chem.FindMolChiralCenters(iso)
                    for c in centers: stereo_info.append(f"({c[1]})")
                    
                    labels.append(f"Isomer {i+1}: " + (", ".join(stereo_info) if stereo_info else "Achiral"))

                # إعدادات الرسم: مش هنفعل الـ Indices عشان ميعملش زحمة
                options = Draw.MolDrawOptions()
                options.prepareMolsBeforeDrawing = True
                
                img = Draw.MolsToGridImage(
                    isomers, 
                    molsPerRow=3, 
                    subImgSize=(400, 400), 
                    legends=labels,
                    drawOptions=options
                )
                st.image(img)

        except Exception as e:
            st.error(f"Error: {e}")
