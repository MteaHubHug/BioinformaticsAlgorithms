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
###############################################################################
###############################################################################
