from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def render_smart_2d(mol):
    if mol is None:
        return None

    # 1. نسخة من الموليكيول
    m = Chem.Mol(mol)
    is_allene = m.HasSubstructMatch(Chem.MolFromSmarts("C=C=C"))
    
    # 2. إضافة الهيدروجين (مهم للألينات والـ Wedges)
    m = Chem.AddHs(m)
    
    # 3. محاولة الـ 3D Embedding بطريقة متوافقة مع كل النسخ
    # جربنا نستخدم EmbedMolecule مباشرة مع تمرير المحاولات كـ Parameter
    if AllChem.EmbedMolecule(m, maxAttempts=5000, randomSeed=42) != -1:
        AllChem.Compute2DCoords(m)
        Chem.WedgeMolBonds(m, m.GetConformer())
    else:
        # لو فشل الـ 3D، ارسم 2D عادي
        AllChem.Compute2DCoords(m)

    # 4. إعدادات الرسم لزيادة الوضوح ومنع الـ Ra/Sa
    d_opts = Draw.MolDrawOptions()
    d_opts.addStereoAnnotation = False  # حذف Ra, Sa, R, S والـ Isomer numbers
    
    # تحسينات شكلية
    if is_allene:
        d_opts.bondLineWidth = 3.0
        d_opts.minFontSize = 20
    else:
        d_opts.bondLineWidth = 1.8
        d_opts.minFontSize = 14

    # 5. توليد الصورة بدون Legend (لإخفاء Search Bitter)
    img = Draw.MolToImage(m, size=(500, 500), options=d_opts, legend="")
    
    return img
