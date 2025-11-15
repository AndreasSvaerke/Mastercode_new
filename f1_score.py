def f1_score(tp, fp, fn):
    """
    Calculate the F1 score given true positives, false positives, and false negatives.

    Parameters:
    tp (int): True Positives
    fp (int): False Positives
    fn (int): False Negatives

    Returns:
    float: F1 score
    """
    if tp == 0:
        return 0.0
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1 = 2 * (precision * recall) / (precision + recall)
    return f1

print(f1_score(169,0,14))  # example usage
