import re
import pandas
from language import *

def fbf_from_peaks(input_filename='lib/combined_fbf.txt'):
    """Returns targetSet object."""
    #fbf1_targets = pandas.read_csv('lib/Table S2 FBF1 to FBF2 N2.txt', sep='\t')
    #fbf1_targets = set(fbf1_targets.Transcript)
    #fbf2_targets = pandas.read_csv('lib/Table S2 FBF2 to FBF2 N2.txt', sep='\t')
    #fbf2_targets = set(fbf2_targets.Transcript)
    #_fbf_targs = fbf1_targets | fbf2_targets
    joint_targets = pandas.read_csv(input_filename, sep='\t')
    _fbf_targs = set(joint_targets.transcript_id)
    fbf_targs = set()
    for t in _fbf_targs:
        s = t.split('.')
        if len(s) > 1:
            s[1] = re.sub('[a-z]$', '', s[1])
            fbf_targs.add('.'.join(s[0:2]))
        else:
            fbf_targs.add(t)
    print("""
FBF targets x.x.x loci: {n}
              x.x loci: {nn}""".format(
        n=len(_fbf_targs), nn=len(fbf_targs)))
    fbf_ts = targetSet(fbf_targs, 'Loci')
    return fbf_ts

