def detect_allene_axes(mol):
    """
    Robust C=C=C detection
    """
    axes = []
    patt = Chem.MolFromSmarts("C=C=C")
    matches = mol.GetSubstructMatches(patt)

    for match in matches:
        left_idx, center_idx, right_idx = match

        left = mol.GetAtomWithIdx(left_idx)
        right = mol.GetAtomWithIdx(right_idx)

        left_subs = [n for n in left.GetNeighbors() if n.GetIdx()!=center_idx]
        right_subs = [n for n in right.GetNeighbors() if n.GetIdx()!=center_idx]

        if len(left_subs)==2 and len(right_subs)==2:
            # تأكد إن substituents مختلفة فعلاً
            if left_subs[0].GetAtomicNum() != left_subs[1].GetAtomicNum() and \
               right_subs[0].GetAtomicNum() != right_subs[1].GetAtomicNum():
                axes.append((left_idx, center_idx, right_idx))

    return axes
