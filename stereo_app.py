import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import py3Dmol
from stmol import showmol

# إعدادات الصفحة
st.set_page_config(page_title="Allene Stereoisomer Visualizer", layout="wide")

st.title("🧪 Allene Stereoisomer Visualizer")
st.write("أدخلي اسم مركب من عائلة الألين (SMILES) لعرض أشكاله الفراغية.")

# مدخلات المستخدم
# مثال: 2,3-pentadiene -> CC=C=CC
smiles_input = st.text_input("أدخلي رمز SMILES للمركب:", value="CC=C=CC")

if smiles_input:
    try:
        # 1. إنشاء المركب الأساسي
        mol = Chem.MolFromSmiles(smiles_input)
        mol = Chem.AddHs(mol)
        
        # توليد الـ Stereoisomers الممكنة
        # ملحوظة: RDKit يتعامل مع الألينات كمحاور كيرالية
        isomers = list(Chem.EnumerateStereoisomers.EnumerateStereoisomers(mol))
        
        st.subheader(f"تم العثور على {len(isomers)} Stereoisomers")
        
        cols = st.columns(len(isomers))
        
        for idx, iso in enumerate(isomers):
            with cols[idx]:
                st.markdown(f"### Isomer {idx + 1}")
                
                # --- العرض 2D (مع Wedge & Dashed) ---
                # نستخدم التنسيق اللي بيوضح الأبعاد
                Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                img = Draw.MolToImage(iso, size=(300, 300), wedgeBonds=True)
                st.image(img, caption=f"2D Projection (Wedge/Hatched)")
                
                # --- العرض 3D ---
                # توليد إحداثيات 3D
                iso_3d = Chem.Mol(iso)
                AllChem.EmbedMolecule(iso_3d, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(iso_3d)
                
                # تحويل المركب لـ Block نصي لـ py3Dmol
                mblock = Chem.MolToMolBlock(iso_3d)
                
                viewer = py3Dmol.view(width=300, height=300)
                viewer.addModel(mblock, 'mol')
                viewer.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
                viewer.zoomTo()
                
                showmol(viewer, height=300, width=300)
                
    except Exception as e:
        st.error(f"حدث خطأ: تأكدي من كتابة SMILES بشكل صحيح. الخطأ: {e}")

# تذكير سريع من مرجعك الخاص
st.sidebar.info("""
**تذكير من دليل المراجع الخاص بكِ:**
* في الألينات، بنستخدم $R/S$ لوصف الترتيب الفراغي للمجموعات حول المحور.
* الـ **Wedge** (المثلث المظلل) يعني المجموعة خارجة ناحيتك.
* الـ **Hatched** (الشرط) يعني المجموعة داخلة لجوه الورقة.
""")
