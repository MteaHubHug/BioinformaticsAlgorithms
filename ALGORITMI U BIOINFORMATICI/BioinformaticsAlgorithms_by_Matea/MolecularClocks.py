import os
import shutil
from DnaReplication import Neighbors
# Motif Finding Problem ==> Finding a motif that is responsible for circadian gene expression

# IMPLANTED MOTIF PROBLEM : Find all (k, d) motifs in a collection of strings
# Input : A collection of strings DNA and integers : k and d
# Output : All (k, d) - motifs in DNA


######################## 2 A ####################################################
# BRUTE FORCE : - generating all k-mers and check which of them are (k, d) - motifs

def find_k_mers(Text, k):
    kmers=[]
    for i in range(len(Text) - k + 1 ):
        komad=""
        z=i
        for j in range(k):
            komad += Text[z]
            z+=1
        kmers.append(komad)
    return kmers



def HammingDistance(p,q):
    if(len(p)!=len(q)):
        print("Greska")
        return -1
    HD= sum(c1 != c2 for c1, c2 in zip(p, q))
    return HD

def differingPatterns(pattern1,pattern2,d):
    HD=HammingDistance(pattern1,pattern2)
    if(HD==-1): return False
    elif(HD>d): return False
    else : return True

def MotifInString(motif,DNA,d):
    for i in range(len(DNA)):
        substr=DNA[i : (i+len(motif))]
        if (len(substr) == len(motif)):
            if(differingPatterns(substr,motif,d)):
                return True
    return False

def MotifInDNA_Mismatches(motif, DNA_collection, d):
    flag=1
    for DNA in DNA_collection:
        if(MotifInString(motif,DNA,d)==False):
           flag=0
    if(flag==1): return True
    else: return False


def MotifEnumeration(DNA_collection , k , d):
    Patterns=[]
    all_kmerPatterns=[]
    for DNA in DNA_collection:
        kmerPatterns=find_k_mers(DNA,k)
        kmerPatterns=set(kmerPatterns)
        kmerPatterns=list(kmerPatterns)
        for kmerPattern in kmerPatterns:
            Neighborhood=Neighbors(kmerPattern,d)
            for neighbor in Neighborhood:
                all_kmerPatterns.append(neighbor)
    all_kmerPatterns=list(set(all_kmerPatterns))
    for kmerPattern in all_kmerPatterns:
        if(MotifInDNA_Mismatches(kmerPattern,DNA_collection,d)):
            Patterns.append(kmerPattern)
    Patterns=set(Patterns)
    return Patterns

DNA_kolekcija= ["AAGGCC","AAGGCT","AAGGTC","AAGGTT"]
k=4
d=1
Paterni= MotifEnumeration(DNA_kolekcija,k,d)
print(" ====> Moguci motivi : " ,Paterni)

