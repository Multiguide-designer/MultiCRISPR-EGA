import warnings


def calcCcTopScore(guideSeq, otSeq):
    """
    calculate the CC top score
    see http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0124633#sec002
    # no mismatch -> most likely off-target
    >>> int(calcCcTopScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    224

    # mismatch in 5' part
    >>> int(calcCcTopScore("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    222

    # mismatch in 3' part
    >>> int(calcCcTopScore("GGGGGGGGGGGGGGGGGGGA","GGGGGGGGGGGGGGGGGGGG"))
    185

    # only mismatches -> least likely offtarget
    >>> int(calcCcTopScore("AAAAAAAAAAAAAAAAAAAA","GGGGGGGGGGGGGGGGGGGG"))
    0
    """
    if len(guideSeq)==23:
        guideSeq = guideSeq[:20]
        otSeq = otSeq[:20]
    elif len(guideSeq) < 20:
        warnings.warn('Running CCtop. The length of the target sequence is less than 20. This sequence will be skipped because it may not be valid for the scoring algorithm.')
        return None
    else:
        guideSeq = guideSeq[:20]
        otSeq = otSeq[:20]
        warnings.warn('Running CCtop. The length of the target sequence is not 23. Maybe your target sequence is not in the format of 20nt guide+3nt PAM, which may cause accuracy problems.')

    if not (len(guideSeq)==len(otSeq)==20):
        raise Exception("Not 20bp long: %s %dbp<-> %s %dbp" % (guideSeq, len(guideSeq), otSeq, len(otSeq)))
    score = 0.0
    for i in range(0, 20):
        if guideSeq[i]!=otSeq[i]:
            score += 1.2**(i+1)
    return 224.0-score

