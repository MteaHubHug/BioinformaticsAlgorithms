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
    print(aminostring,"************")
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
#print(peptid_encoded  # OUTPUT : ['ATGGCC', 'ATGGCC', 'GGCCAT']  ==> OK
###############################################################################
###############################################################################
###############################################################################
