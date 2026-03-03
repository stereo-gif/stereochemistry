from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from PIL import Image # عشان نضمن إن الصورة تفتح في أي نظام

def render_smart_2d(mol):
    # التأكد إن الموليكيول موجود أصلاً عشان ميديناش شاشة بيضا أو Error
    if mol is None:
        print("Error: Molecule object is None!")
        return None

    # 1. نسخة من الموليكيول
    m = Chem.Mol(mol)
    is_allene = m.HasSubstructMatch(Chem.MolFromSmarts("C=C=C"))
    
    # 2. إضافة الهيدروجين (مهم جداً للـ Wedges)
    m = Chem.AddHs(m)
    
    # 3. محاولة الـ 3D Embedding
    # استخدمنا 5000 محاولة (maxAttempts) عشان نضمن إنه ميفشلش ويطلع شاشة بيضا
    params = AllChem.ETKDG()
    params.maxAttempts = 5000 
    
    if AllChem.EmbedMolecule(m, params) != -1:
        AllChem.Compute2DCoords(m)
        Chem.WedgeMolBonds(m, m.GetConformer())
    else:
        # لو فشل الـ 3D، بنرسم 2D عادي كخطة بديلة
        AllChem.Compute2DCoords(m)

    # 4. إعدادات الرسم لزيادة الوضوح ومنع الـ Ra/Sa
    d_opts = Draw.MolDrawOptions()
    d_opts.addStereoAnnotation = False  # حذف Ra, Sa, R, S
    d_opts.prepareMolsBeforeDrawing = True # بيصلح الأخطاء الشائعة قبل الرسم
    
    if is_allene:
        d_opts.bondLineWidth = 3.0
        d_opts.minFontSize = 18
    else:
        d_opts.bondLineWidth = 1.6
        d_opts.minFontSize = 14

    # 5. توليد الصورة بدون Legend (عشان نشيل كلمة Search Bitter)
    img = Draw.MolToImage(m, size=(500, 500), options=d_opts, legend="")
    
    return img

# --- تجربة الكود ---
# جرب تحط الـ SMILES بتاعك هنا مكان ده:
smiles_input = "C/C=C\C" # مثال لمركب Cis-2-Butene (حسب دليلك المرجعي)
my_mol = Chem.MolFromSmiles(smiles_input)

image_result = render_smart_2d(my_mol)

if image_result:
    image_result.show() # هيفتح الصورة في عارض الصور الافتراضي بتاع الويندوز/الماك
