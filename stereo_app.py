import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="Chemical Isomer Analysis", layout="wide")

# (نفس واجهة الـ HTML بتاعتك)
st.markdown("<h2 style='color: #800000;'>Chemical Isomer Analysis System</h2>", unsafe_allow_html=True)

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
                
                # كشف الألّين (Allene) وتنشيط المحور الكايرالي
                patt = Chem.MolFromSmarts("C=C=C")
                if mol.HasSubstructMatch(patt):
                    for match in mol.GetSubstructMatches(patt):
                        mol.GetAtomWithIdx(match[1]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

                mol = Chem.AddHs(mol)
                opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
                isomers = list(EnumerateStereoisomers(mol, options=opts))

                st.success(f"Found **{len(isomers)}** Isomers for {compound_name}")

                cols = st.columns(len(isomers))
                for idx, iso in enumerate(isomers):
                    with cols[idx]:
                        # محاولة بناء الـ 3D
                        AllChem.EmbedMolecule(iso, AllChem.ETKDG())
                        mblock = Chem.MolToMolBlock(iso)
                        
                        # حساب الـ Stereo Label (R/S)
                        Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                        centers = Chem.FindMolChiralCenters(iso, includeUnassigned=True)
                        label = f"Isomer {idx+1} " + (f"({centers[0][1]})" if centers else "")
                        st.subheader(label)

                        # عرض الـ 3D (هنا السحر بيبان)
                        view = py3Dmol.view(width=400, height=300)
                        view.addModel(mblock, 'mol')
                        view.setStyle({'stick': {'colorscheme': 'Jmol', 'radius': 0.15}, 'sphere': {'scale': 0.25}})
                        view.zoomTo()
                        showmol(view, height=300)
                        
                        # رسم 2D للتأكيد تحتها
                        st.image(Draw.MolToImage(Chem.RemoveHs(iso), size=(300, 300)))

        except Exception as e:
            st.error(f"Error details: {e}")
