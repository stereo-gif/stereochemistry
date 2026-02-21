import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

# إعداد واجهة الموقع
st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

st.markdown("""
<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System</h2>
<div style="background-color: #f9f9f9; padding: 15px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000;">Stereoisomerism Reference Guide:</strong><br>
    1. <b>Cis / Trans:</b> Relative position.<br>
    2. <b>E / Z (CIP System):</b> <b>Z (Zusammen)</b> together, <b>E (Entgegen)</b> opposite.<br>
    3. <b>R / S (Optical):</b> Absolute configuration.
</div>
""", unsafe_allow_html=True)

compound_name = st.text_input("Enter Structure Name:", placeholder="e.g., 1,2-dichloroethene")

if st.button("Analyze Isomers"):
    try:
        results = pcp.get_compounds(compound_name, 'name')
        if not results:
            st.error(f"❌ No compound found: {compound_name}")
        else:
            base_smiles = results[0].smiles
            mol = Chem.MolFromSmiles(base_smiles)
            mol_no_stereo = Chem.Mol(mol)
            for bond in mol_no_stereo.GetBonds():
                bond.SetStereo(Chem.BondStereo.STEREONONE)
            
            isomers = list(EnumerateStereoisomers(mol_no_stereo))
            labels = []
            for i, iso in enumerate(isomers):
                Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                stereo_info = [bond.getStereo().name[-1] for bond in iso.GetBonds() if bond.getStereo() in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]]
                chiral_centers = Chem.FindMolChiralCenters(iso)
                for center in chiral_centers:
                    stereo_info.append(f"({center[1]})")
                labels.append(f"Isomer {i+1}: " + (", ".join(stereo_info) if stereo_info else "Achiral"))

            st.success(f"Analyzed Structure: {compound_name} | Total Isomers: {len(isomers)}")
            img = Draw.MolsToGridImage(isomers, molsPerRow=3, subImgSize=(300, 300), legends=labels)
            st.image(img)
    except Exception as e:
        st.error(f"Error: {e}")
