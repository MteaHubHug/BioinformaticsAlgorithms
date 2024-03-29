from reading_tasks_text import *

'''
Antibiotici su kemijski spojevi koji ubijaju bakterije.
Bakterije se bore jedne protiv drugih.
 Mi iskoristavamo bakterije koje napadaju one bakterije koje napadaju nas ljude da bismo iz njih dobili antibiotike
 i pomogli organizmu kojeg zelimo zastiti da se obrani.
Kako bakterije stvaraju antibiotike ?
Bacillus brevis stvara Tyrocidine B1 (10 aminokiselina dug niz)
Te aminokiseline imaju svoja imena i reprezentiramo ih jednim slovom

Kako je Bacillus brevis uspio napraviti taj antibiotik tirocinin
1. Transkripcija DNA niti u RNA nit ( A=>A , C=>C , G=>G , T=>U -uracil)
2. Translacija RNA u aminokiseline, mini-proteine
==> RNA se razbije (particionira se) u kodone
Npr, RNA 3-meri (kodoni) se sastoje od 3 nukleotida, a moguci nukleotidi su A,C,G,U
==> To je 4^3 kombinacija za sastaviti kodon ,
 ali zapravo imamo samo 20 aminokiselina
 jer neke kombinacije (kodoni) rade iste stvari

Mehanizam koji kopira kodone treba znati gdje poceti, a gdje zavrsiti
=> Stoga imamo start- kodone i stop - kodone
=> Postoje 3 stop kodona koja se nikad ne translatiraju u aminokiseline

'''
from collections import OrderedDict
#############################################################################
############################ 4 A ############################################
###############################################################################
# PROTEIN TRANSLATION PROBLEM
# Pronadji translaciju RNA stringa u string aminokiseline
# input : RNA string "Pattern" , niz "GeneticCode"
# output : translacija "Patterna" u string aminokiseline "Peptide"

GeneticCode = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAU": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACU": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGU": "S",
    "AUA": "I",
    "AUC": "I",
    "AUG": "M",
    "AUU": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAU": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCU": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGU": "R",
    "CUA": "L",
    "CUC": "L",
    "CUG": "L",
    "CUU": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAU": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCU": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGU": "G",
    "GUA": "V",
    "GUC": "V",
    "GUG": "V",
    "GUU": "V",
    "UAA": "*",
    "UAC": "Y",
    "UAG": "*",
    "UAU": "Y",
    "UCA": "S",
    "UCC": "S",
    "UCG": "S",
    "UCU": "S",
    "UGA": "*",
    "UGC": "C",
    "UGG": "W",
    "UGU": "C",
    "UUA": "L",
    "UUC": "F",
    "UUG": "L",
    "UUU": "F",
}
RNA_pattern="AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"

def get_peptides(RNA_Pattern):
    peptids=[]
    for i in range(0,len(RNA_Pattern),3):
        peptid= RNA_Pattern[i : i+3]
        if(len(peptid)==3):
           peptids.append(peptid)
    return peptids

def ProteinTranslation(RNA_Pattern, geneticCode):
    Peptide=""
    peptids=get_peptides(RNA_Pattern)
    for peptid in peptids:
        aminokiselina= geneticCode[peptid]
        Peptide+=aminokiselina
    Peptide=Peptide[:-1]
    return Peptide

#peptid= ProteinTranslation(RNA_pattern,GeneticCode)
#print(peptid)
#MAMAPRTEINSTRING ==> OK :)

'''
infile= open_file("rosalind_ba4a.txt")
RNA_pattern=infile[0]
peptid= ProteinTranslation(RNA_pattern,GeneticCode)
print(peptid)

'''
###############################################################################
################################ 4 B ##########################################
###############################################################################
# Peptide Encoding Problem

DNA="ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
amoinoacid="MA"


def ReverseComplement(Pattern):
    Reverse=""
    Pattern = Pattern.upper()
    for nukleotid in Pattern:
        if(nukleotid=="A"): Reverse+="T"
        elif(nukleotid=="C"): Reverse+="G"
        elif(nukleotid=="G"): Reverse+="C"
        else : Reverse+="A"
    Reverse = Reverse[::-1]
    return  Reverse

def DNA2RNA(DNA):
    rna=""
    for nukleotid in DNA:
        if(nukleotid=="A"):
          rna+="A"
        elif(nukleotid=="C"):
          rna+="C"
        elif(nukleotid=="G"):
          rna+="G"
        else:
          rna+="U"
    return rna

def get_start_kodon(aminostring,aminoacid):
    l=len(aminoacid)
    indexi=[]
    #print(aminostring,"************")
    for i in range(len(aminostring)-l+1):
        am=aminostring[i:i+l]
        if(am==aminoacid):
            indexi.append(i)
    return indexi

def PeptideEncoding(DNA,aminoacid):
    dna1=DNA
    dna2=DNA[1:]
    dna3=DNA[2:]

    reverse1=ReverseComplement(dna1)
    reverse2=reverse1[1:]
    reverse3=reverse1[2:]
    rna1=DNA2RNA(dna1)
    rna2 = DNA2RNA(dna2)
    rna3 = DNA2RNA(dna3)
    rna4 = DNA2RNA(reverse1)
    rna5 = DNA2RNA(reverse2)
    rna6 = DNA2RNA(reverse3)

    pepts=get_peptides(rna5)
    amino1=ProteinTranslation(rna1,GeneticCode)
    amino2=ProteinTranslation(rna2,GeneticCode)
    amino3=ProteinTranslation(rna3,GeneticCode)
    amino4 = ProteinTranslation(rna4, GeneticCode)
    amino5 = ProteinTranslation(rna5, GeneticCode)
    amino6 = ProteinTranslation(rna6, GeneticCode)
    aminoacides=[amino1,amino2,amino3]
    aminoacides_from_reversed=[amino4,amino5,amino6]
    res=[]
    for amino in aminoacides:
        indexi=get_start_kodon(amino,aminoacid)
        if(len(indexi)>0):
            for ind in indexi:
                real_ind=int(ind/3)
                pept=DNA[real_ind : real_ind + len(aminoacid)*3]
                res.append(pept)
    aminoacid_reversed=aminoacid[::-1]
    for amino in aminoacides_from_reversed:
        indexi = get_start_kodon(amino, aminoacid_reversed)
        if (len(indexi) > 0):
            for ind in indexi:
                real_ind=len(amino)-ind
                real_ind=int(real_ind/3)
                pept = DNA[real_ind: real_ind + len(aminoacid) * 3]
                pept= ReverseComplement(pept)
                res.append(pept)
    return res

#peptid_encoded=PeptideEncoding(DNA,amoinoacid)
#print(peptid_encoded)  # OUTPUT : ['ATGGCC', 'ATGGCC', 'GGCCAT']  ==> OK
'''

infile= open_file("rosalind_ba4b.txt")
DNA=infile[0]
aminoacid=infile[1]
peptid_encoded=PeptideEncoding(DNA,aminoacid)
for p in peptid_encoded:
    print(p)
    

'''
'''
DNA="TTCGAATCCGTATTGCTTTGGTTGGTTGAGAGAGCGACCTACATTGCTAGTCAGAATAGGTGATTCACACAAAGCTTAGCACCTGGGCAGCACCCCGTGATGTAAACCTATGGGAACTAAGGGAGTCCTGCGGTTTTAGCCAGCAAGCGAGCCGGCAGGAACACTCATACCATCGGACGCGTTTGACGCCTCCCCGGAAAGGAAGTATTTGAGCCTCATTATTACGTATTGCCCGTTAGTCGACAAATCAAGCCCTCGTACGCAGCTTATTCGTACGACGTGGAGGCGTTCCCACGGGCCTAACACGATTGGAACACCACCATAGTAGTGTGGTTCAAATACCTCCTTTGGAGATCTAGAGCTTCACTCTGATTCTAGAGGCAACTTTACAATCGCTCTACGAAATTGTATGGACATCATCAACCGGATATTCTGGGGCGGTAGAATTTCTTTTGTTCGAATCGCTCTAGGCCAGGATCAAATTAATTGAATTGCGGACTCAAGGATCGCGATAGCCGACACATCGGACGCTGTAGAAAGCCAGTCTCTGGATTTAATCCACCCTCTATGTTTGACAAAGCACTAAAACGGGATAGTTTCGGGTGGTATAAGTTTCCCAAGACGATTGCATCGCAATTCATCAACAACCATGAACTTACTGTTTTAGTACTTCCACACACCTTGTTAAATTACGCCTTTACTTCATGTTGCGGTGTGTGTTAGATAGTGTGCAGCTACAAGTCTACCGCCATCGCAGCTCGGGATACCGGCAGATGAGATGGTCCTGAGCTCGTACCGGACTCAAACTTTTTCCTTTACTACCTAGGAATCGCCCATGCGAATTTGTCGGACACACACCATTACATTAACGTCACAACAGCTACTGTTAGAATTTTGCTCTTGCAAATCCTGGAAAGAGTTAAAAAAACTCTTCCGCGCGCCAATAGGGTAAATAATAGATAGCCAGACGGCTGTAAGAGGTGATGACATTTGCAACAATCATGCTGTCGCATCTTCCGCAAGTTCATGTCGCGCCTAGGCAATGGATCTGCGAATGGGGGCCACGGGGTATGAACTACGGAATTCTAAGAAAGTTGCCATCCAGAGTTAAGGGTTTGAGGCTAGTTGCATCGCTGGTAACGAACTACCTCATTACTTGGACGCGCAGTGTGACTTCACTCCTGTATAGCGATGATGCCAAGCAGGAATTAGCAAATCTGAAGAGCGTTTCCAAACTGGCCACTTGGACTGACACCTATCGCGGGGGATTTCAGCGCGTGTCGCTCTCACATGAGAGCTGCCGTCAGGAGCGGTAGAGTTTAGAGAGGAATGCGACAAACTCCCTATTCACCTCTCTGGTGATGTAAGGATATTTACGCTTAGTTCTATGCCAGGCTTAGGGCCTCTCGGAACTTTGGTGAGTCCTTATTAATTGATGCTACCTCTCCCTTACCTTCGCCCCAAGTCACGTAGAAGTACTCAATCCTGCTACATGATAATCAAATATTTCCAACGTTGGGAAATCGGTGACATCACATACTAGTTAAGAAACCACTGTCAGTGAACTTATATCCGGGGGAGAAAATCTACTAACTTACATACGCTGTGCGAGCAGTTTTCATTATAAGAAAATATACTCCCGAGGTACCGCATCAAGCACGACATTCCCGGAGAGCATAACATTTCGGTGCACCTGCTTTTGTGCGCTTGCTTGCGGTTATTTATAAACTACGCACAAGGCGCAAACCGCAGTGCGCATGTTTTCTCCGCCTGGCTAGAACTCGACATTCTCGTCAACGCCAATCTATGTGAGAGGATTTAGACCTCTGTGAAAACGAGTCCCTCTATAGAATAAATACCCAGATGCCAATGGGGGTTCTATCCGATGGCAGTGCATGGAGTGGTGGCTCCAGATTAAGGATGAGGAGAGGTAAAGATAACAGTTCGGTCGCCACGACGCGTTGCCAATCGAAATATCAGTACTAAAAGGCCCACCGCTCCGCTTTAGTCCGACTTTACATCCTGTGGAAATTGTCGAACGGAGGCTACATCGGGCTATATGAGTGTGAAAACCTATACTTCTCGCGTCGTTACTCAGTGCCGGTCTCCTGTTTCCCCCAGTCTTACGTACCCTTATTGATATTTGCTTCACGTTGAAACGTCCTAACGCAGCGTAAAGAGGTGTTTGAACCTCATTACTATAAAATCGCGATCGAAGGTAGACTACATACGCAAACGCCGAAACCCTCAGTTGGCCTTGTTGCAAGTATGGAACGTTGTAAAATTTTTCCTAGACGTTGAGCTATCGGTACAAGGTCGTTAGCGTCCTTACCCTTCACTTATATGCCCGACAAAACGCGGGTCCTAGTGCAGTGGTGGGAGCTTGGAATCCCGCAATACAAGGACAACCTGTATCTCGTTCGGCGTTCCGCGATCACTCGATCCCGAACCACTCCAAGCCTGGTTGATCAGCAAAAGCGGAAGGATGGATAAAGGGCTACTGGTTAATGGATGTAAACTTCCAATGATGAAATCCTGGAAACGAGGGATCGGGTTACGGTGGCGAACGGGGTACGGCAACGTGGCTATCTAGAGCCCGACGTTACGACTCATGTACATGCTGCTACGTGGTTGAAGCTGACGTTCAGATGAAGCAGTACTGAGTCCTAGGGCTTACTACTACTCCAATAGGTCTGGCCGGCCAGATACAAAAGTTCGTGGCGGCTCACCCCCTTTCTGGCGGGTGTAGCTTGCTGACCGGTTTGCTCGATAACACAGGCTAGCGAATAGTAATGAGGTTCGAAAACCTCTTTCCAACGACTGAAAGGGTCTACACGAACTATCTACATTTCCCCGCCCATGTCCTTCCGTCTGGTTGCTTCTGGAGATCCTTTCGCATTATACCGCAGCGTAGTGGCTCTGGCATATATGAAAAAATCCTTCTGTGGGTATTTGTGCCATCACTTATTGTTCGTACCGATATGGGATTACAAGTGCGATGTGATAATAAGCGAAGAAGCCAACATGTTACACTGTTCATGCGCTCCGGGTAATGTGCGGGCACCATGCTCAGTTCCCGCTCGCAGTTGTCACTGTCCCTGTTTCGGCACCATAATCAACATTTCCACGGCCACGCTGGTGAATAACCGAGGATACCGAAGTACAGCAAGAATGAGAGCGGGACTCCTCCATCTGCTTGTAATACGCCTTCAAGATAGTCCATAAAACGGTCGGGGTCTTGTTTCGGACTAGCCGCTTTGAAACGGTGCATAGTTGTGTCAAGTGTGGACATTGGCTTTCTATCCTCGTCAGCGATCCTCGGAAAGACTCGGGCAGTCGCCCCGAATCGTAATTAGGTAGTAGTGCGGCTCAAAAACTTCCTTCGACCTAACCGCTATAATGTTCGTAGATATAAATTTCGTTTCAGTATTAACAGGCGCACCGTATATATACGGAATGGTGTCGCCCCATTAGCTGCTCGCCAATATTTATCTAAGACCGCGCGCGTCTAGCGCCTTTAGTAGTTGCACCCGAGTATAGTAATGGGGTTCGAAGACTTCCTTCGCAAGGCTGCCATACTGTATCACAAGTACTGACGGAGCCCCGGAGGAGTGCAGGATACGGCAAAGGAGACCATTACCGGGGCATGAGTCCAAGTTAGCCCGTTAGGTGAAGGACGCTGATACAATAGTGAATCCGTTACTGAAAGGTTTAGAAGACCGGGGGGCTCGCACTAGGTCCAAATATTATGAACCCTACTCCTGCAACTGAATTGGCCGTCCAGGCGATATTTAAAAGGGGTTACTAGCAGGTTCATCGGAGCCCGTACTCCTTCCGGGCATAGTCGTTCGACGGGTAGAAATTCATCCAGTCGTGCCGGATACCCCGAGAATACCCCTATTTTTTGATCCTTCACCATCATCGTCCGCGGACTCATCTAAGTACCTCAGACCGAAACTGTTATCGTAGCGAAGAGCGAACTCGAATGACATCGCTTGTCCAACAGGGAAAATATGTAAAGTATATGCAGATTATTATAGGAAGATCACAAACTCCATCGCGCCTAGGCCAAAGACTTGCCAGAACAACATCTCTTCCAGAGCAAGGAAGTGTTTGAACCTCACTATTATCGAGAGAAGTCCCATGAATTTATAATAGTGAGGCTCAAAAACTTCCTTCATCGTCGGGCGCTGGGGCGAGCTAGGCTTCCCTAGCCGTCATTACTGTACCCACGCCAAATATTATCGAGTATGACTACGAAGCCTTCACAAGGGGCTGCAGTAACAGACTAACTCACGCACACCGCAACTACTATCCTAAATTGAGTAAGAAAACCCTCGACTATAGCCCTTGAATTAAATTCGCTATTTAAGGAAGACCGCGCTTCCGCTCGCCCGCGTATAGTTTACCAGTCCCAAGAAACATGTGTGGCCAGCCTACCTGAAGAGCTGTGTAAAGATAGATTCTGACATCCTCAAAAAGAAGTTTTCGAGCCGCACTACTACGCACGGAGCTCCGTTATTCAAGGCATGTCGAAGTACAGCGTGGTGGCCTCTCCTGTTCTCCACCCCAGCTAAACCCACCTTGTTCGAATTCGCGCAACTGTATATGACATGAACACTTACAGGGGTTAGAAGTTACTGCAACCAAGATAGAGCTCGTCGAAGTAATAGTGCGGTTCAAAAACTTCCTTCAATTGGTCTCATCACTTAAATTTAAGAGCTATTGTGAGTACAGGTACGGATGCGGCTTCAGTGGATCTTCAGCATTCATTCCTTGTAGTAATGGGGTTCGAAGACTTCCTTGCCAGGGTACCAAACAAGTCTTGCGCATCCTCCTCCCTAAGGAGGTATTTGAACCCCACTATTACCCACGATAGAACATGCAGGGTTTGATAGTGGAACACCTTTTAGAATCTGGGGATAAATTCCCAGGACTAATGTATGGCTGTAGTAATGAGGTTCAAAAACCTCCTTTTCAGGTGGATCGCAGGCCGTGCTGCCTCACAAGCTGGGACGCCGTCCACGGTATAGCCGGCGTCGGCAGTTACTGTGAAATAGCGGAAACTCGATCCCAATATATCATCTTACGTTTGGCGCCCAATAGTCGCCCAGTACCCGTTGACAGTTCTTTAACTCGGCTTAGAACTACTAGACAGGTTCAACCGAACCTTGCCCTAGTTCCCACTCCCGTAATTCATTTGGGTTTGCATCTAGAATACTGGAGGGTGCATAGACGAAACGTGTACGTCGGAGAAAACGTAAGAAAATGACCTAGACTCATAGTAGTGAGGTTCGAAGACTTCCTTTCAGTGAAATCGATCCACCACTCGCCGCGAAGAGATAATAGCATAGAGCACAAGTGCGCGAGTAGAGAAAAAGGCTATCCCAACCGGGCACGTCCTTCGTGTTTGGCGTTTACATACGGCACCCCGTTTCTGCACGTTAACCGTCTAGTATCCAACGGTGGATGGCGGACGCTAGACTATAGATATGAGATATCGAGACCTGGAGCTGGGTGTGGCTGCAGCCCGGGTCATTGCGGGCTGTGAATTCAAGGGCATGTAAACAAGCGTATATCGAACAGTGGATGGGCACCTGCAATACTCACGGTAGAGTTAGCTCACAGGATTCACGTTGAGGACTATGAGTCCCTCTTCGCTAGCAGTCTGGGGGGATATGGAGTTTAATAGCTTGACGTAGTAATGGGGCTCAAACACCTCTTTGTGTGAGCACAGCTACTTGCATTAAGAGATTCTAAACAGCGATCATCTCGGCTATTTCGGGCCAGCCTTTTCGGCAGGATGTTATGTAGCATTTCTGGAAGCTTCCCCCTCGAATCTACTAGTGGTGAGAAGATGCCCCACCGATATTACTCTTTAATCTTGAGAAACCTAAAACCGATCTGACCTCAGACGGGCGGCTCCCACCCGAGGATAAACTCGTCAATAATAATGTGGCTCGAACACTTCTTTTCTCACTAGGCTTTTACGACACGCCACATGTATTTAAGCATCTACCTAACTTGTGTCTGCTGATATACAGCGCATTCTACGCCCAACCTACCAATTACTTCAACGTAGTGCGTGGCTAAAATTCAGGGGAGCTTCATCTCTGTCTTAATTTGAAGGTTCTTCCGGGGCGTTTGGGAATCTTCGTGCCTTTTGCGAGGTTAAGGTATCAAAGAAGTTTTCGAACCACATTATTACCGCCTTAAGCACCGGCGCATCCTGCTCGTGACAACTCTACCCTGCCCTGATAAAGGCACTGAACGTTCCAGAGAGTGCATCATTGACACGCGAGCAGGCCACAGTAGCCACAAGACGTATGGGTGATTATAGAATTGGTGGAGGTGTTGTTAACGATCAGGAGGACATTAGTGGGAGTTAGGAAAGACCCTATGTTCTCTCTATCGCGGACTTGTAACTTGACAAGCAAAAGGGTAAGAGAGCTGCACACCGAAGCAGGCCCTTCCTATACCTGTTTTTCCTACGCGTAGAGAGGAATCCAGAAAGGTGATAATTGGCATTCGATGAAAAAACAGTGTGCCACTGACTTAGTTCTATATGTGAAGAGCCTGTTAGCACGTGACGGCGGCCTTGGTATAGAGCCCTTAATGGTCTCCATCGCGTAGTAATGGGGCTCGAAAACCTCCTTACTTGGGATTGCGTGGCCTCCTTGTGAGTCATACACAAGGCTTAGGGCTATGGGGCGATACACTCCTTTTCGCGGCGCATGGGGCGGTGATGCCTACATAGTAGTAGTGACTGCCTTTCTGGGGGGCTATTTGTGGATGACCAACACCTGACCAGCGATGCAATCGCTAGGGGAGGTACACCTCTCATATGTTACAACAATCACCGAATTGTGTTTCGAATTCGAATCAAGTTTGCGGTGTCGACCAGATCTGGTCTTGCTGCCATACCGGGTTCGCCGCCTCCGGTGGATAGAACTGCATCTTAAGACATCTGGACCCAGCGGTAAGTAGCGGGAAGAGTTTAGAGTCATTCGTACAACTACAGGCTAAGGGCTTACTGGGGAGTTGTTGTAGGGCATAAAGATCGCCCCATGACTTTTCGTACTTTCCCCGATAGTTCACTCGCAGCGAGCTGCGGCTGGGCTTCGCCACACGAGTACGGGCAACATTTATCTCCTCTAATCACTGGGCACCGCGCGAGGAAATAGAAAACCCTAATCAGTGCTCATGGGCGCATCTATTGGTCTCCGCATGCACGATGCCGCGGAGTGCTTAGTTGTCCCTGCATAATCTTCGTAGATGTATAAGAGATTACCTATTTATTCGGTTTCGGTTCTAGACGTACCTTGCCGCATGAGTATAGGCTAATGAACTGAGTTGGCGCCAGAGGGAAAGGCATAATAATGCGGCTCGAATACTTCCTTAAGGAAGTATTCGAACCACATTACTAT"
amoinoacid="KEVFEPHYY"
peptid_encoded=PeptideEncoding(DNA,amoinoacid)
for p in peptid_encoded:
    print(p)
'''
###############################################################################
############################### 4 C############################################
###############################################################################
'''
Tirocinini i gramicidini su ciklicki peptidi
=> Jedna ciklicka reprezentacija ima puno linearnih reprezentacija
jer ne znamo odakle krenuti u krugu
Ali, kada rijesimo Peptid Encoding Problem ne pronadjemo niti jedan 
30-mer u Bacillus brevis koji je enkodiran - slika 4.4 
Centralna dogma molekularne biologije
 ==> svi peptidi moraju biti enkodirani s genomom
 Ribosom radi protein translation =>
 Ako umrtvimo ribosom, onda nema tko praviti proteine, 
 medjutim, iako smo ubili radnika (ribosoma), 2 proteina su se nastavila proizviditi
 : tirocidin i gramicidin 
 => Postoji mehanizam koji ne ovisi o ribosomu,  a sastavlja te peptide
 => Tirocidin i gramicidini su sintetizirani NRP sintetazom
 => Taj nezim  slaze antibioticke peptide bez obzira na RNA i geneticki kod
 => Sastavlja peptide tako sto sastavlja jednu po jednu aminokiseline
 => Postoje posebni radnici koji sastavljaju aminokiseline i onda 
 drugi radnik bira aminokiseline i sastavlja peptid
 
 Buduci da NRP-ovi ne slijede centralnu dogmu, ne mozemo pronaci stvari iz samoga genoma
 ==> Kako onda mozemo peptide koje smo dobili istrazivanjem mozemo 
 analizirati od kojih se aminokiselina sastoje? 
 
   Glavni stroj koji se koristi u sekvencioniranju peptida : maseni spektrometar 
   Mjeri masu molekula u daltonima (Da) ... prije toga razbije molekulu na atome
   Mi cemo aproksimirati masu molekule tako da zbrajamo protone i neutrone
   u molekuli (njezinim atomima)
   C_2H_3ON ~> 57 ; 2*12 + 3*1 + 1*16 + 1*14 = 57 
   
   Mi znamo aminokiseline i njihove mase
   Maseni spektrometar generira puno fragmenata 
   
   npr : imamo VEM i dobijemo :
   M , VE , VEM , EM , V ... 
   => kako imamo i mase, imamo eksperimentalni spektar
   
   *Cyclopeptide Sequencing Problem*
   pretpostavka :Maseni spektrometar razbija kopije ciklickih peptida na svakoj mogucoj vezi
   tako da rezultat sardzi mase svih mogucih linearnih fragmenata peptida
    (potpeptidi; subpeptides)
    Npr: 
    Ciklicki peptid NQEL ima 12 subpeptida :
    N,Q,E,L, NQ, QE, EL , LN , NQE,QEL,ELN, LNQ 
    i pretpostavljamo da se svaki moze pojaviti vise puta 
    => npr : za ELEL imamo isto 12 potpeptida : 
    E, L , E, L, EL, LE, EL, LE, ELE, LEL, ELE, LEL
    
    Theretical spectrum (teoretski spektar ) ciklickog peptida je 
    kolekcija masa svih njegovih subpeptida uz dodatak mase 0 i mase cijeloga peptida
    i zapisujemo mase od minimalne do maksimalne 
    
    Cyclospectrum(Peptide) - teoretski spektar peptida
    - taj teoretski spektar moze sadrzavati duplikate npr. 
    u NQEL ima potpeptide NQ i EL koji imaju jednake mase : 242  
   
'''
# Generating Theoretical Spectrum Problem :
# Generiranje teoretskog spektra nekog ciklickog peptida
# input : Peptid aminokiseline (string)
# output : Cyclospectrum(Peptide)

'''
za input NQEL , trebam dobiti output : 
    L   N   Q   E   LN  NQ  EL  QE  LNQ ELN  QEL NQE NQEL
|----------------------------------------------------------|
 0 113 114 128 129 227 242 242  257 355 256  370 371 484


-treba iskorisiti tablicu masa 
-razbiti na sve moguce kombinacije
-iscitamo njihove mase 
- ako imamo dva slova, onda zbrajamo njihove odgovarajuce mase
- i onda to sve slozimo u rastucem redoslijedu
'''

mass = {
    "G": 57,
    "A": 71,
    "S": 87,
    "P": 97,
    "V": 99,
    "T": 101,
    "C": 103,
    "I": 113,
    "L": 113,
    "N": 114,
    "D": 115,
    "K": 128,
    "Q": 128,
    "E": 129,
    "M": 131,
    "H": 137,
    "F": 147,
    "R": 156,
    "Y": 163,
    "W": 186,
}
Aminoacid_peptid_example="NQEL"
Aminoacid_peptid_example_rosalind="LEQN"
def get_all_linear_fragments_of_cyclicPeptide(Aminoacid_peptid):
    l=len(Aminoacid_peptid)
    linear_fragments=[]
    for i in range(l):
        a= Aminoacid_peptid[i:l] + Aminoacid_peptid[0 : i]
        linear_fragments.append(a)
    return linear_fragments

def get_every_subset_family(Aminoacid_peptid):
    linear_fragments=get_all_linear_fragments_of_cyclicPeptide(Aminoacid_peptid)
    n = len(Aminoacid_peptid)
    arr = []
    for a in linear_fragments:
        for i in range(0, n):
            for j in range(i, n):
                fragment=a[i:(j + 1)]
                if(len(fragment)<n):
                   arr.append(fragment)
    arr.append(Aminoacid_peptid)
    arr=list(set(arr))
    arr.sort()
    return arr


def GeneratingTheoreticalSpectrum(Aminoacid_peptid, mass_table):
    subpeptides=get_every_subset_family(Aminoacid_peptid)
    theoretical_spectrum={}
    for subpeptid in subpeptides:
        theoretical_spectrum[subpeptid]=0

    for key in theoretical_spectrum:
        for letter in key:
            mass=mass_table[letter]
            theoretical_spectrum[key]+=mass
    theoretical_spectrum["zero"]=0
    theoretical_spectrum={k: v for k, v in sorted(theoretical_spectrum.items(), key=lambda item: item[1])}
    mase=[]
    for key in theoretical_spectrum:
        mase.append(theoretical_spectrum[key])
    return mase


teoretski_spektar=GeneratingTheoreticalSpectrum(Aminoacid_peptid_example_rosalind,mass)
#print(teoretski_spektar)
# moj putput za NQEL :
# : {'zero': 0, 'L': 113, 'N': 114, 'Q': 128, 'E': 129, 'LN': 227, 'EL': 242, 'NQ': 242, 'QE': 257, 'LNQ': 355, 'ELN': 356, 'QEL': 370, 'NQE': 371, 'NQEL': 484}

# ROSALIND example :
# LEQN
# 0 113 114 128 129 227 242 242 257 355 356 370 371 484
# moj  output : [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484] ==> OK

'''
infile= open_file("rosalind_ba4c.txt")
teoretski_spektar=GeneratingTheoreticalSpectrum(infile[0],mass)

print( ' '.join(map(str, teoretski_spektar))  )

'''
'''

rosalind_primjer = "IAQMLFYCKVATN"
teoretski_spektar=GeneratingTheoreticalSpectrum(rosalind_primjer,mass)
print( ' '.join(map(str, teoretski_spektar))  )

'''
##################################################################################################################
###################################### 4 D ########################################################################
##################################################################################################################
'''# Cyclopeptide Sequencing Problem  - obrisi
# Za dani Ideal Sprectrum, nadji ciklicki peptide čiji teoretski spektar matcha eksperimentalni spektar
# Input : kolekcija (moguća ponavljanja) integera "Spectrum" koji odgovaraju Ideal Spectrumu
# Output : String aminokiseline "Peptide" t.d. CycloSpectrum(Peptide)==Spectrum  (ako takav postoji)
'''
'''def CycloSpectrum(Peptide):
    res=GeneratingTheoreticalSpectrum(Peptide,mass)
    return res






def BFCyclopeptideSequencing(Spectrum):
    for peptide in Spectrum:
        if(Mass(peptide,mass) == ParentMass(Spectrum)):
            if(Spectrum== CycloSpectrum(peptide)):
                return peptide'''


def get_all_masses(masses):
    res= masses.values()
    res=list(set(res))
    res.sort()
    return res

def get_all_combos(mass, all_masses, ctr):
    if mass in ctr:
        return ctr[mass]
    if mass < 0:   return 0
    if mass == 0:  return 1

    ctr[mass] = 0
    for aa_mass in all_masses:
        ctr[mass] += get_all_combos(mass - aa_mass, all_masses, ctr)

    return ctr[mass]

def CountingPeptidesWithGivenMass(m):
    all_masses=get_all_masses(mass)
    counter={}
    res=1
    res= get_all_combos(m,all_masses,counter)
    return res
m=1024
#combosNum=CountingPeptidesWithGivenMass(m)
#print(combosNum) # 14712706211 ==> OK

'''
infile= open_file("rosalind_ba4d.txt")
m=int(infile[0])
print(m)
combosNum=CountingPeptidesWithGivenMass(m)
print(combosNum)

'''
'''
m=1307
combosNum=CountingPeptidesWithGivenMass(m)
print(combosNum)

'''
###########################################################################################################################
################################################### 4 E ####################################################################
##########################################################################################################################
'''CYCLOPEPTIDESEQUENCING(Spectrum)
        Peptides ← a set containing only the empty peptide
        while Peptides is nonempty
            Peptides ← Expand(Peptides)
            for each peptide Peptide in Peptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Cyclospectrum(Peptide) = Spectrum
                        output Peptide
                    remove Peptide from Peptides
                else if Peptide is not consistent with Spectrum
                    remove Peptide from Peptides'''
def Mass(Peptide): # total mass
    total_mass=0
    for p in Peptide:
        total_mass += p
    return total_mass

def ParentMass(Spectrum):
    max_mass=max(Spectrum)
    return max_mass

def expand(peptides,masses):
    return [peptide+[mass] for peptide in peptides for mass in masses]

def linearSpectrum(peptide):
    def get_pairs():
        return [(i,j) for i in range(len(peptide)) for j in range(len(peptide)+1) if i<j]
    result=[sum(peptide[i:j]) for (i,j) in get_pairs()]
    result.append(0)
    result.sort()
    return result

def find_cyclopeptide_sequence(spectrum):
    def isConsistent(peptide):
        # Determine number of items in spect that match specified element
        def count(element, spect):
            return len([s for s in spect if s == element])

        peptide_spectrum = linearSpectrum(peptide)
        for element in peptide_spectrum:
            if count(element, peptide_spectrum) > count(element, spectrum):
                return False
        return True

    def cycloSpectrum(peptide):
        def get_pairs(index_range):
            return [(i, j) for i in index_range for j in range(i, i + len(index_range)) if j != i]
        augmented_peptide = peptide + peptide
        result = [sum(augmented_peptide[a:b]) for (a, b) in get_pairs(range(len(peptide)))]
        result.append(0)
        result.append(sum(peptide))
        result.sort()
        return result

    peptides = [[]]
    output = []
    masses = list(set(mass.values()))

    while len(peptides) > 0:
        next_peptides = []
        for peptide in expand(peptides, masses):
            if Mass(peptide) == ParentMass(spectrum):
                if cycloSpectrum(peptide) == spectrum:
                    output.append(peptide)
            else:
                if isConsistent(peptide):
                    next_peptides.append(peptide)
        peptides = next_peptides

    return output

spektar=[0, 113, 128, 186, 241, 299, 314, 427]
#rez=find_cyclopeptide_sequence(spektar)
# 186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186
#print(rez)


'''
infile= open_file("rosalind_ba4e.txt")
spektar=infile[0].split(" ")
spektar = [int(i) for i in spektar]
rez=find_cyclopeptide_sequence(spektar)
strings=[]
for l in rez:
    li=[str(i) for i in l]
    stringy="-".join(li)
    strings.append(stringy)

print(" ".join(strings))

'''