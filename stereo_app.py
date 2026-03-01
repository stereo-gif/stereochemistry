import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, EnumerateStereoisomers
import pubchempy as pcp
import py3Dmol
from stmol import showmol

# إعدادات واجهة Streamlit
st.set_caption("Allene Chemistry Tool")
st.title("🧪 مستكشف الألينات الفراغي")
st.markdown("---")

# مدخل البحث بالاسم
compound_name = st.text_input("أدخلي اسم مركب الألين (بالإنجليزي):", placeholder="مثال: 2,3-pentadiene")

if compound_name:
    try:
        # 1. تحويل الاسم إلى SMILES باستخدام PubChem
        with st.spinner('جاري البحث عن بنية المركب...'):
            results = pcp.get_compounds(compound_name, 'name')
            
        if results:
            smiles = results[0].isomeric_smiles
            mol_base = Chem.MolFromSmiles(smiles)
            mol_base = Chem.AddHs(mol_base)
            
            # 2. توليد الـ Stereoisomers (الـ R والـ S للألين)
            opts = EnumerateStereoisomers.StereoisomerEnumerationOptions(tryEmbedding=True)
            isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol_base, options=opts))
            
            st.success(f"تم العثور على المركب! يوجد {len(isomers)} أشكال فراغية محتملة:")
            
            # عرض الأشكال في أعمدة
            cols = st.columns(len(isomers))
            
            for idx, iso in enumerate(isomers):
                with cols[idx]:
                    # تحديد الـ Stereochemistry (R/S) بناءً على مرجعك
                    Chem.AssignStereochemistry(iso, force=True, cleanIt=True)
                    
                    st.subheader(f"الشكل رقم {idx + 1}")
                    
                    # --- الرسم 2D (Hatched & Dark Wedge) ---
                    # الـ RDKit هيرسم الألين بجهة مسطحة وجهة فيها روابط Wedge/Dashed
                    img = Draw.MolToImage(iso, size=(400, 400), wedgeBonds=True, kekulize=True)
                    st.image(img, caption="رسم 2D يوضح الأبعاد (Wedge/Hatched)")
                    
                    # --- الرسم 3D التفاعلي ---
                    st.write("النموذج المجسم (3D):")
                    iso_3d = Chem.Mol(iso)
                    AllChem.EmbedMolecule(iso_3d, AllChem.ETKDG())
                    AllChem.MMFFOptimizeMolecule(iso_3d)
                    
                    mblock = Chem.MolToMolBlock(iso_3d)
                    viewer = py3Dmol.view(width=350, height=350)
                    viewer.addModel(mblock, 'mol')
                    viewer.setStyle({'stick': {'colorscheme': 'cyanCarbon'}, 'sphere': {'radius': 0.3}})
                    viewer.zoomTo()
                    showmol(viewer, height=350, width=350)
                    
        else:
            st.error("لم أتمكن من العثور على مركب بهذا الاسم. تأكدي من السبلنج.")
            
    except Exception as e:
        st.error(f"حدث خطأ أثناء المعالجة: {e}")

# تذكير من المرجع الخاص بك
with st.sidebar:
    st.header("📌 تذكير القواعد")
    st.info("""
    **حسب الـ Reference Guide الخاص بكِ:**
    * **Cis/Trans:** تستخدم للوصف النسبي.
    * **E/Z:** تستخدم لوصف الأولوية حسب (CIP System).
    * **R/S:** في الألينات، يتم ترتيب المجموعات الـ 4 حسب الأولوية وتحديد الدوران حول المحور.
    """)
