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
    1. <b>Cis / Trans:</b> Relative position.<br>
    2. <b>E / Z (CIP System):</b> Based on atomic number priority (1 > 2).<br>
    3. <b>R / S:</b> Absolute configuration of chiral centers.
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
                st.error(f"❌ No compound found for: {compound_name}")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # توليد الأيزومرات
                mol_no_stereo = Chem.Mol(mol)
                for bond in mol_no_stereo.GetBonds():
                    bond.SetStereo(Chem.BondStereo.STEREONONE)
                
                isomers = list(EnumerateStereoisomers(mol_no_stereo))
                
                # التحقق إذا كان المركب Achiral
                if len(isomers) <= 1:
                    st.info(f"✨ The compound **{compound_name}** is **Achiral** (No stereoisomers found).")
                
                labels = []
                for i, iso in enumerate(isomers):
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    
                    # تفعيل إظهار أرقام الذرات (Indices) للمساعدة في فهم الأولوية
                    for atom in iso.GetAtoms():
                        atom.SetProp('atomNote', str(atom.GetIdx()))
                    
                    stereo_info = []
                    # استخراج E/Z
                    for bond in iso.GetBonds():
                        stereo = bond.GetStereo()
                        if stereo == Chem.BondStereo.STEREOE: stereo_info.append("E")
                        elif stereo == Chem.BondStereo.STEREOZ: stereo_info.append("Z")
                    
                    # استخراج R/S
                    centers = Chem.FindMolChiralCenters(iso)
                    for c in centers: stereo_info.append(f"({c[1]})")
                    
                    label = f"Isomer {i+1}: " + (", ".join(stereo_info) if stereo_info else "Achiral")
                    labels.append(label)

                st.success(f"Analyzed: **{compound_name}**")
                
                # إعدادات الرسم لإظهار الأرقام
                draw_config = Draw.MolDrawOptions()
                draw_config.addAtomIndices = True # ده هيحط رقم كل ذرة جنبها
                
                img = Draw.MolsToGridImage(
                    isomers, 
                    molsPerRow=3, 
                    subImgSize=(400, 400), 
                    legends=labels,
                    useSVG=False
                )
                st.image(img)

        except Exception as e:
            st.error(f"Error: {e}")
