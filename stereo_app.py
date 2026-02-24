import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

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
            results = pcp.get_compounds(compound_name, 'name')
            if not results:
                st.error("Compound not found.")
            else:
                base_smiles = results[0].smiles
                mol = Chem.MolFromSmiles(base_smiles)
                
                # --- الخطوة السحرية للألّينات (Allenes) ---
                # بنورّر على ذرة الكربون اللي في النص (مرتبطة بـ 2 روابط ثنائية)
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == 'C':
                        double_bonds = [b for b in atom.GetBonds() if b.GetBondType() == Chem.rdchem.BondType.DOUBLE]
                        if len(double_bonds) == 2:
                            # بنعلم الذرة دي إنها محتمل تكون "Centric" للألّين
                            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW) 
                
                # تنظيف الـ Stereo لإعادة الحساب
                mol_no_stereo = Chem.Mol(mol)
                Chem.AssignStereochemistry(mol_no_stereo, force=True, cleanIt=True)
                
                # استخدام خيار الـ Embedding لمحاكاة الشكل الفراغي المتعامد للألّين
                options = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol_no_stereo, options=options))
                # ---------------------------------------

                labels = []
                for i, iso in enumerate(isomers):
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    info = []
                    
                    # كشف الـ R/S (بما في ذلك الـ Axial chirality للألّين)
                    centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                    for c in centers:
                        info.append(f"({c[1]})")
                        
                    for b in iso.GetBonds():
                        if b.GetStereo() == Chem.BondStereo.STEREOE: info.append("E")
                        elif b.GetStereo() == Chem.BondStereo.STEREOZ: info.append("Z")
                    
                    label = f"Isomer {i+1}: " + (", ".join(set(info)) if info else "Achiral")
                    labels.append(label)

                st.success(f"Found **{len(isomers)}** forms for {compound_name}")
                
                img = Draw.MolsToGridImage(isomers, molsPerRow=2, subImgSize=(400, 400), legends=labels)
                st.image(img, use_container_width=True)

        except Exception as e:
            st.error(f"Error: {e}")
