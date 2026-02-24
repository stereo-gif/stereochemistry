import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from stmol import showmol
import py3Dmol

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
                
                # --- التعديل الجذري للألّين ---
                # 1. بنحدد روابط الألّين ونخليها "قابلة للتغيير" (Unspecified Stereo)
                patt = Chem.MolFromSmarts("C=C=C")
                matches = mol.GetSubstructMatches(patt)
                for match in matches:
                    # بنمسك الروابط الثنائية في الألّين ونخليها نشطة فراغياً
                    bond1 = mol.GetBondBetweenAtoms(match[0], match[1])
                    bond2 = mol.GetBondBetweenAtoms(match[1], match[2])
                    bond1.SetStereo(Chem.BondStereo.STEREONONE)
                    bond2.SetStereo(Chem.BondStereo.STEREONONE)
                
                # 2. بنضيف هيدروجينات لأنها أساس الـ Axial Chirality
                mol = Chem.AddHs(mol)
                
                # 3. التوليد مع إجبار فحص الـ 3D (Embedding)
                opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol, options=opts))

                # لو لسه مطلع واحد (وده نادر مع الكود ده)، بنجرب نعكس الـ Stereo يدوياً
                if len(isomers) == 1:
                    st.info("Generating enantiomers via reflection...")
                    iso2 = Chem.Mol(isomers[0])
                    for atom in iso2.GetAtoms():
                        if atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                        elif atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                    isomers.append(iso2)

                st.success(f"Found **{len(isomers)}** Isomers")

                cols = st.columns(len(isomers))
                for idx, iso in enumerate(isomers):
                    with cols[idx]:
                        # بناء الـ 3D
                        AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                        mblock = Chem.MolToMolBlock(iso)
                        
                        # حساب النوع (R/S)
                        Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                        label = f"Isomer {idx+1}"
                        st.subheader(label)

                        # عرض الـ 3D
                        view = py3Dmol.view(width=400, height=300)
                        view.addModel(mblock, 'mol')
                        view.setStyle({'stick': {'radius':0.2}, 'sphere': {'scale': 0.3}})
                        view.zoomTo()
                        showmol(view, height=300)
                        
                        # رسم الـ 2D مع الـ Wedges
                        st.image(Draw.MolToImage(Chem.RemoveHs(iso), size=(300, 300), wedgeBonds=True))

        except Exception as e:
            st.error(f"Error: {e}")
