import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

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
                
                # --- الخطوة الحاسمة: تعريف الألّين يدوياً كمركز كايرالي ---
                for atom in mol.GetAtoms():
                    # بندور على الكربون اللي في نص الألّين
                    if atom.GetSymbol() == 'C':
                        # الكربون ده لازم يكون عامل رابطتين ثنائيتين
                        double_bonds = [b for b in atom.GetBonds() if b.GetBondType() == Chem.rdchem.BondType.DOUBLE]
                        if len(double_bonds) == 2:
                            # بنجبر RDKit يعتبره مركز "كايرالي" عشان الـ Enumerate يشوفه
                            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                            # بنعلم عليه عشان يبان في الحسابات
                            atom.SetProp("_ChiralityPossible", "1")

                # تنظيف وتجهيز
                mol = Chem.AddHs(mol)
                Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
                
                # التوليد (بإجبار البرنامج يفك كل المراكز اللي حددناها)
                options = StereoEnumerationOptions(onlyUnassigned=False, tryEmbedding=True)
                isomers = list(EnumerateStereoisomers(mol, options=options))

                labels = []
                final_mols = []
                for i, iso in enumerate(isomers):
                    # تنظيف الرسم
                    res_mol = Chem.RemoveHs(iso)
                    Chem.AssignStereochemistry(res_mol, force=True, cleanIt=True)
                    
                    # قراءة النوع (R/S) أو الـ Parity
                    centers = Chem.FindMolChiralCenters(res_mol, includeUnassigned=True)
                    info = []
                    for c in centers:
                        # في الألّين الـ R/S بيتحسب كـ Axial Chirality
                        info.append(f"({c[1]})")
                    
                    label = f"Isomer {i+1}: " + (", ".join(info) if info else "Achiral Structure")
                    labels.append(label)
                    final_mols.append(res_mol)

                st.success(f"Found **{len(isomers)}** forms for {compound_name}")
                
                # الرسم (مهم نستخدم SVG عشان الـ Wedges والـ Dashes تبان في الألّين)
                img = Draw.MolsToGridImage(final_mols, molsPerRow=2, subImgSize=(400, 400), legends=labels, useSVG=True)
                st.write(img, unsafe_allow_html=True)

        except Exception as e:
            st.error(f"Error: {e}")
