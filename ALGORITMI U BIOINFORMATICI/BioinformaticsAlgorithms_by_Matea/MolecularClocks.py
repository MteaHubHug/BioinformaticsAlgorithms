from reading_tasks_text import *



# BRUTE FORCE - ideja :
# Za bilo koji (k,d)-motiv mora postojati max d mismatches u dna za taj specifican k-mer
# mozemo definirati sve k-mere i onda provjeriti je li neki od njih (k,d)-motiv

#################################### 2 A
## ZADATAK : ako nam je dano vise stringova, pronadjimo moguce (k,d)-motive :
# (k,d)-motiv ==> oznaka da imamo substring duzine k koji je slican drugome stringu do maksimalno d mismatches
import math
import os
from itertools import chain
import random
import numpy
import numpy as np

import DnaReplication

sample_dataset = ["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]
k = 3
d = 1

from DnaReplication import Neighbors


def findKmers(Text, k):
    kmers = []
    for i in range(len(Text) - k + 1):
        komad = ""
        z = i
        for j in range(k):
            komad += Text[z]
            z += 1
        kmers.append(komad)
    return kmers


def HammingDistance(p, q):
    if (len(p) != len(q)):
        print("Ne moze se mjeriti Hammingova distanca izmedju unijetih stringova", p, " i ", q,
              " jer nisu iste duljine.")
        return -1
    HD = sum(c1 != c2 for c1, c2 in zip(p, q))
    return HD


def ApproximatePatternCount(dna, pattern, d):
    # nadji broj patternovih susjeda (s najvise d razlika) u tekstu "dna"
    cnt = 0
    kmeri = findKmers(dna, len(pattern))
    susjedi = Neighbors(pattern, d)
    susjedi = list(set(susjedi))
    for kmer in kmeri:
        if kmer in susjedi:
            cnt += 1
    return cnt


def WordsWithMismatches(text, k, d):  ## ==> trazi sve k-mere s najvise d razlika u stringu text
    kmeri_texta = findKmers(text, k)  # text = " ATTTGGC"
    #print(kmeri_texta)
    svi_susjedi = []
    kmeri_frekvencije = {}
    for kmer in kmeri_texta:
        susjedi = Neighbors(kmer, d)
        svi_susjedi.append(susjedi)
    svi_susjedi = list(chain.from_iterable(svi_susjedi))
    svi_susjedi = list(set(svi_susjedi))
    for susjed in svi_susjedi:
        kmeri_frekvencije[susjed] = ApproximatePatternCount(text, susjed, d)
    return svi_susjedi


def MotifEnumeration(DNA_kolekcija, k, d):
    Patterns = []
    Patterns_rez = []
    for dna in DNA_kolekcija:
        Patterns.append(WordsWithMismatches(dna, k, d))
    Patterns = list(chain.from_iterable(Patterns))
    is_in_all_dnas = len(DNA_kolekcija)
    Patterns2 = list(set(Patterns))
    Patterns_dict = {}
    for pattern in Patterns2:
        Patterns_dict[pattern] = Patterns.count(pattern)
        if (Patterns_dict[pattern] >= is_in_all_dnas): Patterns_rez.append(pattern)
    Patterns_rez = list(set(Patterns_rez))
    return Patterns_rez


# motivEnur=MotifEnumeration(sample_dataset,3,1)
# print(motivEnur)

#infile= open_file("rosalind_ba2a.txt")
#sample_dataset=[]
#for stringy in infile[1:]:
#    sample_dataset.append(stringy)
#nums=infile[0].split(" ")
#k=int(nums[0])
#d=int(nums[1])
#motivEnur=MotifEnumeration(sample_dataset,k,d)
#print(motivEnur)
#for mot in motivEnur:
#    print(mot)

#########################
"""
Mi jos ne znamo odakle dolazi taj originalni motiv => tj. ne znamo jos koji je "konsenzus" ::: 73. str u knjizi 
Mozemo ga naci pomocu matrica Scorea, Counta...

Entropija mjeri koliko su stvari balansirane. Sto je vise konzerviran skup, entropija je manja ==> Daje nam nacin da scoreamo MotifMatrices
Entropija matrice == suma svih entropija po stupcima 
"""

"""
Motif Finding Problem : 
Ako nam je dana kolekcija (niz) stringova, pronadjimo sve k-mere - jedan po stringu, koji minimiziraju Score(Rezultantni Motiv)
Input : skup stringova dna, k
Output : kolekcija motiva k-mera, jedan iz svakog stringa dna koji minimiziraju Score (tako da to ide po svim mogucim izborima kmers) 

Brute force : 
Imali bi (n-k+1)^t razlicitih nacina za sloziti motive 
n- duljina stringa, k - duljina kmera, t - broj sekvenci  ==> slozenost : O(n^t * k * t)  ==> treba brzi algoritam 

Umjesto da generiramo konsenzusa, ajmo krenut od konsenzusa ... 
Treba bolja strategija za racunati Score. 
Neka imamo kolekciju k-mer motiva i k-mer Pattern 
--> mozemo definirati udaljenost izmedju Patterna i skupa motiva kao sumu svih HammingDistances izmedju Patterna i svakoga motiva
--> to moze zapravo biti Score

==> Umjesto da trazimo kolekciju k-mera koji minimiziraju Score(Motiv),
 ajmo traziti kolekciju k-mera za potencijalni konsenzus string koji minimizira d(Pattern,Motifs)


 Ekvivaletni Motif Finding Problem : 
 Ako nam je dan skup stringova, pronadji Pattern i kolekciju k-mera (jednog iz svakog stringa) 
 koji minimiziraju udaljenost izmjedji svih mogucih Patterna i svih mogucih skupova k-mera

 ==> Ako nam je dan Pattern, ne moramo istraziti sve motive da bismo minimizirali funkciju d(Pattern,Text)
 ° definirajmo Motifs(Patterns,DNA) kao skup svih k-mera koji minimiziraju udaljenost d za dani Pattern i za sve moguce Motive
"""


#################################### 2 B
##### ZADATAK : PRonadji k-mer Pattern koji minimizira d(Pattern,Dna) preko svih k-mer Patterna
######### Takav k-mer zovemo Median string za Dna
######### Input :Niz Dnaova i k (integer)
######## Output : k-mer Pattern koji minimizira d(Pattern,Dna) izmedju svih Patterna koji se mogu postaviti


def distance_Pattern_dna_substring(Pattern, dna_substring):
    kmers = findKmers(dna_substring, len(Pattern))
    kmer_HA_dict = {}
    for kmer in kmers:
        HA = HammingDistance(kmer, Pattern)
        kmer_HA_dict[kmer] = HA
    motif = min(kmer_HA_dict, key=kmer_HA_dict.get)
    distance = kmer_HA_dict[motif]
    return distance


def Motif(Pattern, dna_substring):
    kmers = findKmers(dna_substring, len(Pattern))
    kmer_HA_dict = {}
    for kmer in kmers:
        HA = HammingDistance(kmer, Pattern)
        kmer_HA_dict[kmer] = HA
    motif = min(kmer_HA_dict, key=kmer_HA_dict.get)
    return motif


'''Pattern="GATTCTCA"
dnastring="GCAAAGACGCTGACCAA"
motiv=Motif(Pattern,dnastring)
print(motiv)
dist=distance_Pattern_dna_substring(Pattern,dnastring)
print(dist)

output : GACGCTGA , 3'''


def d_distance_Pattern_DnaCollection(Pattern, Dna_kolekcija):
    cnt = 0
    for dna in Dna_kolekcija:
        dist = distance_Pattern_dna_substring(Pattern, dna)
        cnt += dist
    return cnt


'''dna_kolekcija= ["TTACCTTAAC","GATATCTGTC","ACGGCGTTCG","CCCTAAAGAG","CGTCAGAGGT"]
Pattern="AAA"
print(d_distance_Pattern_DnaCollection(Pattern,dna_kolekcija))
# OUTPUT : 5 '''

k = 3
dna_kolekcija = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCGT", "GCTGAGCACCGG", "AGTACGGGACAG"]
# output bi trebao biti : GAC

'''Nadji k-mer Pattern koji minimizira d(Pattern,dna_kolekcija)
 tako da prodjes kroz sve k-mer Patterne i dnaove u Dna_kolekciji za koje d=min Hamming Distance'''


def MedianString(Dna_kolekcija, k):
    distance = float("inf")
    motivi = []
    for dna in Dna_kolekcija:
        kmeri = findKmers(dna, k)
        kmeri = list(set(kmeri))
        for kmer in kmeri:
            motif = Motif(kmer, dna)
            motivi.append(motif)
    motivi_distance_dict = {}
    for motiv in motivi:
        dist = d_distance_Pattern_DnaCollection(motiv, dna_kolekcija)
        motivi_distance_dict[motiv] = dist
    median = min(motivi_distance_dict, key=motivi_distance_dict.get)
    return median


### median=MedianString(dna_kolekcija,k)
### print(median)   # REZ = "GAC"

#infile= open_file("rosalind_ba2b.txt")
#dna_kolekcija=[]
#for stringy in infile[1:]:
#    dna_kolekcija.append(stringy)

#k=int(infile[0])
#median=MedianString(dna_kolekcija,k)
#print(median)

#################################################################################################################
####################################### 2 C  #########################################################################
#################################################################################################################
"""
Greedy Motif Search problem - greedy je pohlepni algoritam 
greedy algoritam je pohlepan jer odabire rjesenje koje mu se najvise svidja po nekoj definiciji svidjanja
"""
dna_string = "ACTAACCAAGCA"
k = 4


def GetMatrix(dna_string, k):
    dna = []
    dna[:0] = dna_string
    l = int(len(dna_string) / k)
    dna = numpy.array(dna).reshape(l, k)
    return dna


def GetCountMatrix(mca, k):
    shape = (4, k)
    mca1=np.array(mca)
    mca_shape = mca1.shape
    m = mca_shape[0]
    # n=mca_shape[1]
    count_matrix = numpy.empty(shape, dtype=int, order='C')
    for j in range(k):  # nekad ovdje treba ici n , zavisi od zadatka
        cntA = 0
        cntC = 0
        cntG = 0
        cntT = 0
        for i in range(m):
            if (mca[i][j] == "A"): cntA += 1
            if (mca[i][j] == "C"): cntC += 1
            if (mca[i][j] == "G"): cntG += 1
            if (mca[i][j] == "T"): cntT += 1
        count_matrix[0][j] = cntA
        count_matrix[1][j] = cntC
        count_matrix[2][j] = cntG
        count_matrix[3][j] = cntT
    return count_matrix


def ProfilMatrix(matrica, k):
    # l = int(len(matrica) / k)
    # matrica= GetMatrix(dna_string,k)
    Count_matrix = GetCountMatrix(matrica, k)
    sume = []
    for j in range(k):
        suma = 0
        for i in range(4):
            suma += Count_matrix[i][j]
        sume.append(suma)
    shape = (4, k)
    profil = numpy.empty(shape, dtype=float)
    for i in range(4):
        for j in range(k):
            profil[i][j] = round((Count_matrix[i][j] / sume[j]), 2)
    return profil


def Probability(motiv, profil, k):
    p = 1
    for i in range(k):
        if (motiv[i] == "A"):
            p *= profil[0][i]
        elif (motiv[i] == "C"):
            p *= profil[1][i]
        elif (motiv[i] == "G"):
            p *= profil[2][i]
        elif (motiv[i] == "T"):
            p *= profil[3][i]
    return p


# matrica= GetMatrix(dna_string,k)
# print(matrica)
#l = int(len(dna_string) / k)
# count=GetCountMatrix(matrica,k, l)
# print(count)

# profil=ProfilMatrix(count,k,l)
# print(profil)
# prob=Probability("ACTA",profil,k)
# print(prob)
"""
[['A' 'C' 'T' 'A']
 ['A' 'C' 'C' 'A']
 ['A' 'G' 'C' 'A']]

[[3 0 0 3]
 [0 2 2 0]
 [0 1 0 0]
 [0 0 1 0]]

[[1.   0.   0.   1.  ]
 [0.   0.67 0.67 0.  ]
 [0.   0.33 0.   0.  ]
 [0.   0.   0.33 0.  ]]
"""


def ProfileMostProbableKmerProblem(text, k, profil):
    kmeri = findKmers(text, k)
    probabilities_kmers = {}
    for kmer in kmeri:
        p = Probability(kmer, profil, k)
        probabilities_kmers[kmer] = p
    mostProbableKmer = max(probabilities_kmers, key=probabilities_kmers.get)
    return mostProbableKmer

#infile= open_file("rosalind_ba2c.txt")
#text=infile[0]
#k=int(infile[1])
#prob_rows=[]
#for stringy in infile[2:]:
#    splitted=stringy.split(" ")
#    probs=[]
#    for spl in splitted:
#        probs.append(float(spl))
#    prob_rows.append(probs)
#profil=prob_rows
#rez=ProfileMostProbableKmerProblem(text, k, profil)
#print(rez)
#############################################################################################
###################################### 2 D ##################################################
#############################################################################################
# mostProbableKmer=ProfileMostProbableKmerProblem(dna_string,k,profil)
# print(mostProbableKmer)

dna_kolekcija = ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"]


def makeMatrix(motifs):
    new_motifs = []
    for insideList in motifs:
        for motif in insideList:
            list1 = []
            list1[:0] = motif
            new_motifs.append(list1)
    new_motifs = numpy.array(new_motifs)
    return new_motifs


def Score(Motifs):
    # Motifs=makeMatrix(Motifs)
    Motifs = np.array(Motifs)
    score = 0
    shape = Motifs.shape
    m = shape[0]
    n = shape[1]
    for j in range(n):
        cntA = 0
        cntC = 0
        cntG = 0
        cntT = 0
        count_columns = {}
        for i in range(m):
            if (Motifs[i][j] == "A"):
                cntA += 1
            elif (Motifs[i][j] == "C"):
                cntC += 1
            elif (Motifs[i][j] == "G"):
                cntG += 1
            elif (Motifs[i][j] == "T"):
                cntT += 1
        count_columns["A"] = cntA
        count_columns["C"] = cntC
        count_columns["G"] = cntG
        count_columns["T"] = cntT
        most_popular = max(count_columns, key=count_columns.get)
        count = count_columns[most_popular]
        diff = m - count
        score += diff
    return score


t = len(dna_kolekcija)


def razlozi(input):
    lista = []
    lista[:0] = input
    return lista


def InitializeBestMotifs(Dna_kolekcija, k):
    mca = []
    for dna in Dna_kolekcija:
        motif = dna[0: (0 + k)]
        motif = razlozi(motif)
        mca.append(motif)
    mca = numpy.array(mca)
    return mca


def GreedyMotifSearch(Dna_kolekcija, k, t):
    BestMotifs = InitializeBestMotifs(Dna_kolekcija, k)
    kmeri1 = findKmers(Dna_kolekcija[0], k)
    shape = (t, 1)
    for kmer in kmeri1:
        Motifs = []
        motiv = razlozi(kmer)
        Motifs.append(motiv)
        for i in range(1, t):
            Profile = ProfilMatrix(Motifs, k)
            Motif = razlozi(ProfileMostProbableKmerProblem(Dna_kolekcija[i], k, Profile))
            Motifs.append(Motif)
        if (Score(Motifs) < Score(BestMotifs)):
            BestMotifs = Motifs
    return BestMotifs


# BestMotifs=GreedyMotifSearch(dna_kolekcija,k,t)
# print(BestMotifs)
########## OUTPUT : [['C', 'C', 'T', 'T'], ['G', 'A', 'T', 'A'], ['A', 'C', 'G', 'G'], ['C', 'C', 'T', 'A'], ['C', 'A', 'G', 'A']]
############ za k=3 : OUTPUT :: [['A', 'A', 'C'], ['G', 'A', 'T'], ['A', 'C', 'G'], ['A', 'A', 'G'], ['G', 'A', 'G']]
'''
infile= open_file("rosalind_ba2d.txt")
nums=infile[0].split(" ")
k=int(nums[0])
t=int(nums[1])
dna_kolekcija=[]
for stringy in infile[1:]:
    dna_kolekcija.append(stringy)

BestMotifs=GreedyMotifSearch(dna_kolekcija,k,t)
BestMotifs2=[]
for mot in BestMotifs:
    motif=""
    for character in mot:
        motif+=character
    BestMotifs2.append(motif)
    print(motif)
    '''
##################################################################################################################################
##################################################    2 E   #######################################################################
##################################################################################################################################

dna_kolekcija = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT"]

# Moze se dogoditi da vjerojatnost da se neki string razlikuje za samo jedan znak, ali da je vjerojatnost =0 ,
# Da bismo izbjegli "unfair scoring" cesto zamijenjujemo nule u Profilnoj matrici s nekim vrlo malim brojevima (blizu nule) koje zovemo *pseudocounts*
# Ili u Count matrici svakom elementu dodajemo 1
# MotifsMatrix= [["TAAC"],["GTCT"],["ACTA"],["AGGT"]]

def odvoji(matrica):
    nova_mca = []
    for red in matrica:
        novi_red = razlozi(red)
        nova_mca.append(novi_red)
    nova_mca = np.array(nova_mca)
    return nova_mca


def CountPlus1(Dna_kolekcija, k):
    mca = odvoji(Dna_kolekcija)
    mca = GetCountMatrix(mca, k)
    mca = np.array(mca)
    nova_mca = mca
    shape = mca.shape
    m = shape[0]
    n = shape[1]
    for i in range(m):
        for j in range(n):
            nova_mca[i][j] = mca[i][j] + 1
    return nova_mca


def ProfilMatrixPlus(Dna_kolekcija, k):
    Count_matrix = CountPlus1(Dna_kolekcija, k)
    sume = []
    for j in range(k):
        suma = 0
        for i in range(4):
            suma += Count_matrix[i][j]
        sume.append(suma)
    shape = (4, k)
    profil = numpy.empty(shape, dtype=float)
    for i in range(4):
        for j in range(k):
            profil[i][j] = round((Count_matrix[i][j] / sume[j]), 2)
    return profil


'''count_plus=CountPlus1(dna_kolekcija,k)
print(count_plus)
profile=ProfilMatrixPlus(dna_kolekcija,k)
print(profile)'''


def GreedyMotifSearch1(Dna_kolekcija, k, t):
    BestMotifs = InitializeBestMotifs(Dna_kolekcija, k)
    kmeri1 = findKmers(Dna_kolekcija[0], k)
    shape = (t, 1)
    for kmer in kmeri1:
        Motifs = []
        motiv = razlozi(kmer)
        Motifs.append(motiv)
        for i in range(1, t):
            Profile = ProfilMatrixPlus(Motifs, k)
            Motif = razlozi(ProfileMostProbableKmerProblem(Dna_kolekcija[i], k, Profile))
            Motifs.append(Motif)
        if (Score(Motifs) < Score(BestMotifs)):
            BestMotifs = Motifs
    return BestMotifs


# bestMotifs=GreedyMotifSearch1(dna_kolekcija,k,len(dna_kolekcija))
# print(bestMotifs) # OUTPUT : [['A', 'C', 'C', 'T'], ['A', 'T', 'G', 'T'], ['A', 'C', 'G', 'G'], ['A', 'C', 'G', 'A'], ['A', 'G', 'G', 'T']]
## JEDNAKO KAO U KNJIZI - TOCNO RJESENJE :)

def Consensus(
        Motifs):  # input: [['A', 'C', 'C', 'T'], ['A', 'T', 'G', 'T'], ['A', 'C', 'G', 'G'], ['A', 'C', 'G', 'A'], ['A', 'G', 'G', 'T']]
    Motifs = np.array(Motifs)
    shape = Motifs.shape
    m = shape[0]
    n = shape[1]
    consensus = ""
    for j in range(n):
        cnts = {}
        cntA = 0
        cntC = 0
        cntG = 0
        cntT = 0
        for i in range(m):
            znak = Motifs[i][j]
            if (znak == "A"):
                cntA += 1
            elif (znak == "C"):
                cntC += 1
            elif (znak == "G"):
                cntG += 1
            elif (znak == "T"):
                cntT += 1
            cnts["A"] = cntA
            cnts["C"] = cntC
            cnts["G"] = cntG
            cnts["T"] = cntT
        izabrani_znak = max(cnts, key=cnts.get)
        consensus += izabrani_znak
    return consensus

# bestMotifs=GreedyMotifSearch1(dna_kolekcija,k,len(dna_kolekcija))
# cons=Consensus(bestMotifs)
# print(cons) ########### OUTPUT : ACGT ==> TOCNO KAO U KNJIZI =))


#################################################################################################################################################
#######################################################   2 F        ############################################################################
#################################################################################################################################################
# Randomizirani algoritmi - neintuitivni => Las Vegas algoritmi, Monte Carlo algoritmi
# Las Vegas algoritmi - spada u klasu algoritama u kojoj je ono do cega dodjemo garantirano da ce biti tocno iako se temelji na random idejama
# Primjer : trazimo najvecu vrijednost funkcije, a funkcije je specificna po tome da ima samo jedan maksimum (lokalni maksimum koji je ujedno i globalni)
######## ==> Tada je svejedno ocemo li poceti traziti maksimum lijevo od njega ili desno od njega - garantirano je da cemo ga pronaci  (skica : ._._/\_._ )
# Ti algoritmi zahtjevaju puno pretpostavki i ne daju nuzno tocno rjesenje, brzi su u pronalazenju aproksimativnog rjesenja
# ==> Spadaju u klasu Monte Carlo algoritama
# Takodjer mozemo imati vise lokalnih maksimuma i jedan globalni, pa ako postavimo vise "tragaca maksimuma",
# svaki ce naci neki maksimum, ne mora znaciti da ce taj koji su nasli biti globalni, ali ako postavimo puno tragaca,
# vjerojatno ce jedan zavrsiti u globalnom maksimumu ili nekom lokalnom koji je vrlo blizu po vrijednosti globalnom

dna_kolekcija = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTCG", "CCCTAACGAG", "CGTCAGAGGT"]
k = 4
mca_dna_kolekcija = odvoji(dna_kolekcija)
Profil = [[4 / 5, 0, 0, 1 / 5], [0, 3 / 5, 1 / 5, 0], [1 / 5, 1 / 5, 4 / 5, 0], [0, 1 / 5, 0, 4 / 5]]
Profil = np.array(Profil)


def Motifs(Profile, Dna_kolekcija):
    # mca=odvoji(Dna_kolekcija)
    best_kmers = []
    k = Profile.shape[1]
    for dna in Dna_kolekcija:
        kmeri = findKmers(dna, k)
        kmers_probabilities = {}
        for kmer in kmeri:
            kmer1 = razlozi(kmer)
            p = 1
            for i in range(k):  # <==> for nukleotid in kmer1:
                if (kmer1[i] == "A"):
                    p *= Profile[0][i]
                elif (kmer1[i] == "C"):
                    p *= Profile[1][i]
                elif (kmer1[i] == "G"):
                    p *= Profile[2][i]
                elif (kmer1[i] == "T"):
                    p *= Profile[3][i]
            kmers_probabilities[kmer] = p
        # print(kmers_probabilities)
        best_kmer = max(kmers_probabilities, key=kmers_probabilities.get)
        best_kmers.append(best_kmer)
    return best_kmers


# best_kmers=Motifs(Profil,dna_kolekcija)
# print(best_kmers)  # OUTPUT : ['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT'] ==> isto kao u knjizi ==> tocno :)

def RandomMotifs(Dna_kolekcija, k):
    t = len(dna_kolekcija)
    kmeri = []
    for dna in Dna_kolekcija:
        randy = random.randint(0, (len(dna) - k - 1))  # nekad ne treba dna[0] vec samo dna
        kmer = dna[randy: (randy + k)]
        kmeri.append(kmer)
    return kmeri


def RandomizedMotifSearch(Dna_kolekcija, k, t):
    motifs = RandomMotifs(Dna_kolekcija, k)
    motifs = odvoji(motifs)
    BestMotifs = motifs
    while (True):
        Profile = ProfilMatrixPlus(motifs, k)
        motifs = Motifs(Profile, Dna_kolekcija)
        motifs = odvoji(motifs)
        if (Score(motifs) < Score(BestMotifs)):
            BestMotifs = motifs
        else:
            return BestMotifs


# t=len(dna_kolekcija)
# bestMotifs=RandomizedMotifSearch(dna_kolekcija,k,t)
# print(bestMotifs)

'''
OUTPUT : 
[['A' 'C' 'C' 'T']    ==> ok 
 ['A' 'T' 'G' 'T']    ==> ok
 ['A' 'C' 'G' 'G']    ==> nok :/ ==> vidimo da nas je brzina kostala malo tocnosti, ali je rjesenje i dalje dosta dobro (samo jedan pogresan kmer) 
 ['A' 'C' 'G' 'A']    ==> ok
 ['A' 'G' 'G' 'T']]   ==> ok'''

"""
k=8
t=5
Rosalind_dataset=["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
rez=RandomizedMotifSearch(Rosalind_dataset,8,5)
OUTPUT : 
[['T' 'C' 'T' 'C' 'G' 'G' 'G' 'G'] ==> ok 
 ['T' 'G' 'T' 'A' 'A' 'G' 'T' 'G'] ==> nok
 ['T' 'A' 'C' 'A' 'G' 'G' 'C' 'G'] ==> ok
 ['T' 'T' 'C' 'A' 'G' 'G' 'T' 'G'] ==>ok
 ['T' 'G' 'T' 'T' 'G' 'G' 'C' 'C']] ==> nok     ==> 3/5 tocno ...==> kada bismo pokrenuli 100000 puta, sigurno bismo dobili najvise puta tocan rezultat
                                                                te bismo onda samo izabrali one koji se najcesce pojavljuju i stvorili konsenzus
                                                                 ==> ja sam pokrenula samo jednom i vec je bilo 3/5 tocno
 """


'''
infile= open_file("rosalind_ba2f.txt")
nums=infile[0].split(" ")
k=int(nums[0])
l=int(nums[1])
dna_kolekcija=[]
for stringy in infile[1:]:
    dna_kolekcija.append(stringy)
rez=RandomizedMotifSearch(dna_kolekcija,k,l)

rez2=[]
for l in rez:
    motif=""
    for character in l:
        motif+=character
    rez2.append(motif)
    print(motif)

'''
######################################################################################################################################################
################################################################ 2 G #################################################################################
######################################################################################################################################################
# GIBBS SAMPLER - GIBBSOVO UZORKOVANJE
# Ideja : Kada radimo Randomized Motif Search, tada mozemo sve motive promijeniti => To je mozda prevelika promijena
# U Randomized Motif Search, mi iz kolekcije Dna konstruiramo Profilnu matricu
# i s tom profilnom matricom trazimo nove motive koji najvise odgovaraju toj profilnoj matrici
#####
# Gibbs Sampler pristup : mijenja se samo jedan k-mer u svakom iteracijskom koraku
# Npr, biramo na slucajan nacin red (jedan dna iz kolekcije)
# ==> recimo da smo izabrali drugi red, konstruirali bismo profilnu matricu od ostalih dna iz kolekcije i s tom bismo matricom promijenili
######## izabrali k-mer u drugom redu matrice dna kolekcije (drugom dna stringu u kolekciji) ... itd
# Dakle, Gibbs Sampler krece s nasumicno odabranim k-merima u sekvencama. U svakoj iteraciji na slucajan nacin bira jednu DNA sekvencu koju ce promijeniti
# Koristi nasumicno odabrane k-mere "Motifs (to je niz motiva)" za doci do novoga, (nadamo se boljega) skupa k-mera

# Primjer distribucije simetricne kocke : (1/6 , 1/6 , 1/6 , 1/6 , 1/6 , 1/6 ) ... Primjer distribucije pristrane kocke (0.1, 0.2, 0.3, 0.5, 0.1, 0.25)
#################################################################################### Suma vjerojatnosti nije nuzno =1, mozemo normalizirati vjerojatnosti
#################################################################################### ==> sumiramo vjerojatnosti i podijelimo svaku vjv s dobivenim zbrojem
dna_kolekcija = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT"]
rosalind_dataset = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
                    "ACCAGCTCCACGTGCAATGTTGGCCTA"]
k = 8
t = 5
N = 100


def getMotifMatrixWithoutRandomIndex(motifs, index):
    t = len(motifs)
    novi_motifs = []
    for i in range(t):
        if (i != index):
            novi_motifs.append(motifs[i])
    return novi_motifs


''' nepotrebno - za izbrisati
def NormalizeProfileMatrix(Profile_matrix):
    shape = np.array(Profile_matrix).shape
    m = shape[0]
    n = shape[1]
    nova_profilna_matrica=Profile_matrix
    sume=[]
    for j in range(n):
        suma=0
        for i in range(4):
            suma+=Profile_matrix[i][j]
        sume.append(suma)

    for j in range(n):
        for i in range(4):
            nova_profilna_matrica[i][j]=nova_profilna_matrica[i][j]/sume[j]
    return nova_profilna_matrica



def Profil_generate_most_probable_motif(Profile_matrix):
    Profile_matrix=NormalizeProfileMatrix(Profile_matrix)
    shape=np.array(Profile_matrix).shape
    m=shape[0]
    n=shape[1]
    generirani_kmer = ""
    for j in range(n):
        max=Profile_matrix[0][j]
        maxi=0
        znak="x"
        for i in range(4):
            if(max<Profile_matrix[i][j]):
                max=Profile_matrix[i][j]
                maxi=i
            if(maxi==0): znak="A"
            elif(maxi==1): znak="C"
            elif(maxi==2): znak="G"
            elif(maxi==3): znak="T"
        generirani_kmer+=znak
    return generirani_kmer
'''


def Score_as_distance(Motifs):
    consensus = Consensus(Motifs)
    score = d_distance_Pattern_DnaCollection(consensus, Motifs)
    return score


def ProfileMostProbableKmerProblemNormalization(text, k, profil):
    kmeri = findKmers(text, k)
    probabilities_kmers = {}
    for kmer in kmeri:
        p = Probability(kmer, profil, k)
        probabilities_kmers[kmer] = p
    suma = sum(probabilities_kmers.values())
    for kljuc in probabilities_kmers:
        probabilities_kmers[kljuc] = probabilities_kmers[kljuc] / suma
    mostProbableKmer = max(probabilities_kmers, key=probabilities_kmers.get)
    # r = random.randint(0, len(kmeri) - 1)
    # mostProbableKmer = kmeri[r]
    return mostProbableKmer


def GibbsSampler(Dna_kolekcija, k, t, N):
    motifs = RandomMotifs(Dna_kolekcija, k)
    motifs = odvoji(motifs)
    BestMotifs = motifs
    for j in range(1, N):
        i = random.randint(0, t - 1)
        motifs_without_i = odvoji(getMotifMatrixWithoutRandomIndex(motifs, i))
        Profile = ProfilMatrixPlus(motifs_without_i, k)
        the_kmer = list(ProfileMostProbableKmerProblemNormalization(Dna_kolekcija[i], k, Profile))
        # the_kmer= list(  Profil_generate_most_probable_motif(Profile))
        motifs[i] = the_kmer
        if (Score(motifs) < 10): print("*******************************", Score(motifs))
        if (Score(motifs) < Score(BestMotifs)):
            BestMotifs = motifs
    print("BEST SCORE :", Score(BestMotifs))
    return BestMotifs



''' OUTPUT 
# ovaj algoritam nekad treba izvrititi vise puta i onda uzeti u obzir samo one rezultate kada je score barem 10... :)) 
BEST SCORE : 10
[['C' 'T' 'C' 'G' 'G' 'G' 'G' 'G']
 ['C' 'C' 'A' 'A' 'G' 'G' 'T' 'G']
 ['T' 'A' 'C' 'A' 'G' 'G' 'C' 'G']
 ['T' 'T' 'C' 'A' 'G' 'G' 'T' 'G']
 ['T' 'C' 'C' 'A' 'C' 'G' 'T' 'G']]


 BEST SCORE : 10
[['C' 'T' 'C' 'G' 'G' 'G' 'G' 'G']
 ['C' 'C' 'A' 'A' 'G' 'G' 'T' 'G']
 ['T' 'A' 'C' 'A' 'G' 'G' 'C' 'G']
 ['T' 'T' 'C' 'A' 'G' 'G' 'T' 'G']
 ['T' 'C' 'C' 'A' 'C' 'G' 'T' 'G']]
'''

'''
infile= open_file("rosalind_ba2g.txt")
nums=infile[0].split(" ")
k=int(nums[0])
t=int(nums[1])
times=int(nums[2])
dna_kolekcija=[]
for stringy in infile[1:]:
    dna_kolekcija.append(stringy)
best_motivs_Gibbs = GibbsSampler(dna_kolekcija, k, t,times)
rez2=[]
for l in best_motivs_Gibbs:
    motif=""
    for character in l:
        motif+=character
    rez2.append(motif)
    print(motif)

'''

#################################################################################################################################################
########################################### 2 H ################################################################################################
#################################################################################################################################################

def DistanceBetweenPatternAndString(Pattern, Dna_collection):
    k = len(Pattern)
    distance = 0
    for dna in Dna_collection:
        HammingovaDistanca = float("inf")
        kmeri = findKmers(dna, k)
        for kmer in kmeri:
            if (HammingovaDistanca > HammingDistance(Pattern, kmer)):
                HammingovaDistanca = HammingDistance(Pattern, kmer)
        distance = distance + HammingovaDistanca
    return distance


def MedianString2(Dna_collection, k):
    distance = float('inf')
    Median = "x"
    for i in range(int(math.pow(4, k)) - 1):  # ovdje mozda treba jos -1
        Pattern = DnaReplication.NumberToPattern(i, k)
        if distance > DistanceBetweenPatternAndString(Pattern, Dna_collection):
            distance = DistanceBetweenPatternAndString(Pattern, Dna_collection)
            Median = Pattern
    return Median

# medijan_string= MedianString2(dna_kolekcija,4)
# print(medijan_string)
