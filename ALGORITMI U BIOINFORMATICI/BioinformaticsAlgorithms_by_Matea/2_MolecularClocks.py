import os
import shutil

# Motif Finding Problem ==> Finding a motif that is responsible for circadian gene expression

# IMPLANTED MOTIF PROBLEM : Find all (k, d) motifs in a collection of strings
# Input : A collection of strings DNA and integers : k and d
# Output : All (k, d) - motifs in DNA

# BRUTE FORCE : - generating all k-mers and check which of them are (k, d) - motifs

def MotivEnumeration(DNA , k , d):
    Patterns=[]
    #for kmerPattern in DNA:

