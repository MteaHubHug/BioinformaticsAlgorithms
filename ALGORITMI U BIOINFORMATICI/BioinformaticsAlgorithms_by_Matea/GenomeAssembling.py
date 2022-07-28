# Newspaper problem - mozemo li od puno fragmenata novina sastaviti novine ?
# Teze nego sto se cini - imamo dijelice puno kopija istog izdanja jednog primjerka
# Neke informacije su izgubljene, nemamo kompletne novine, samo dijelove
# Ne mozemo "zalijepiti" stranice, trebamo koristiti ideju Overlappinga ==> preklapanje fragmenata iz razlicitih kopija da bismo rekonstruirali primjerak
## Paralela s biologijom : neka dna izgleda ovako : ACGCTTTACGT... itd ==> trebamo naci nacin da prepoznamo da je to dna nekg covjeka
# ==> Kako doci do tih nukleotida koji predstavljaju taj Dna? ==> *sekvencioniranje* - pokusanje citanja sekvence
# ==> Problem je sto imamo mogucnost citati samo odredjene komadice (tehnoloski) -- recimo da je nase ogranicenje da mozemo procitati samo 5 nukleotida
# Recimo da imamo samo komadice i to razbacane po prostoru od zadanog dna gore : npr : TTTAC  ; ... ; GT ; ... ; ACGC ... (nisu ni iste duljine)
############################################################################## => svi su duljine <=5 tako da ih mozemo procitati s obzirom na ogranicenje
# Dakle, uzme se stanica ili vise stanica te osobe, te stanice imaju isti Dna, onda se biokemijski taj Dna odvoji te se "PCR" (Polymerase chain reaction)
#### metodom replicira ==> stvore se kopije toga ==> onda se biokemijskim reakcijama razbiju kopije u komadice
# ==> pa se ti komadici mogu ( s obzirom na ogranicenje tehnologije ) procitati
# Tada se dobije mix svih tih komadica, ali mi onda ne znamo jesu li neki od tih komadica kopije, idu li prije ili poslije odredjenog komada itd..
# Zadatak je : Pokusati sastaviti originalan Dna
# CILJ : sastaviti komadice i dobiti pocetni dna string
### Overleap ==> Ako imamo komadic jedne kopije : ACGC i komadic druge kopije : GCTT , ta dva komadica se preklapaju na "GC".


# String Composition Problem - Generiraj k-mer kompoziciju stringa
# INPUT : string Text i integer k
# OUTPUT : Composition_k(Text) gdje su k-meri poredani po abecedi
import numpy

from DnaReplication import *
from MolecularClocksMotifs import *


################################################################################################################################
######################################### 3 A ##################################################################################
################################################################################################################################
dna_text="TATGGGGTGC"
k=3
rosalind_sample="CAATCCAAC"
k=5
def Composition(Text,k):
    kmeri=findKmers(Text,k)
    kmeri.sort()
    print(kmeri)
    return kmeri

#kmeri=Composition(rosalind_sample,k) # ['ATG', 'GGG', 'GGG', 'GGT', 'GTG', 'TAT', 'TGC', 'TGG'] ==> ok :)
# rosalind - output : ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA'] ==> ok :)

# Ono sto nas zanima nije dobiti kmere iz stringa, vec iz kmera dobiti Dna string , zapravo...==> String Reconstruction Problem
# Input : integer k, kolekcija k-mera "Patterns"
# Output  : string Text (dna) s kompozicijom k-mera jednakom kao sto je kolekcija "Patterns" iz inputa

# Recimo da imamo ove k-mere : AAT, ATG, GTT, TAA, TGT , pretpostavimo da se overleapaju u svim simbolima osim u max. jednom
# Dakle, ocekujemo da se AAT i ATG preklapaju u AT
# Krenuli bismo od TAA, jer se TA ne moze preklapati ni sa cim od ponudjenih, a onda se AA iz TAA moze preklapati s AAT recimo...
# PRIMJER :
    # TAA
    #  AAT
    #   ATG
    #    TGT
    #     GTT
# ==> TAATGTT

# Malo komliciraniji PRIMJER kada se moze vise kmera pojaviti kao odgovarajuci sljedbenik promatranog kmera :
# kmeri : "AAT", "ATG", "ATG", "ATG", "CAT", "CCA", "GAT", "GCC", "GGA", "GGG", "GTT", "TAA", "TGC", "TGG", "TGT"

    # TAA
    #  AAT
    #   ATG
#==>  TAATG
# sada imamo vise izbora koji pocinju s TG : to su : TGC, TGG, TGT
# izabiremo onaj koji nas navodi na to da potrosimo sve elemente koje imamo u danoj kolekciji kmera
# a)  TAA               TAA
#      AAT     ==>       AAT
#       ATG               ATG
#        TGT               TGT ==> nema dalje ==> NE

# b)  TAA               TAA
#      AAT     ==>       AAT
#       ATG               ATG
#        TGC               TGC
#                           GCC
#                            CCA
#                             CAT
#                              ATG
#                               TGG
#                                GGA
#                                 GAT
#                                  ATG
#                                   TGT
#                                    GTT
# ==>                   TAATGCCATGGATGTT  ==> Iskoristeno je 14/15 kmera ==> ne valja
# Pokusajmo naci nacin da iskoristimo i petnaesti
# Ali, vidimo da se ATG pojavljuje 3 puta ==> 3 opcije za nastavke ==> sad vidimo koliko je komplicirano i koliko nacina ima
# ===> pokusajmo rijesiti grafom


############################################################################################################################################
############################################### 3 B ########################################################################################
############################################################################################################################################
# Kako od stringa napraviti graf ?  " STRING SPELLED BY A GENOME PATH PROBLEM "
# Ideja sufixa i prefixa


kolekcija_dijelica=["ACCGA","CCGAA","CGAAG","GAAGC","AAGCT"]


def StringSpelledByGenomePathProblem(kolekcija):
    res=kolekcija[0]
    for i in range(1,len(kolekcija)):
        res+=kolekcija[i][-1:]
    return res

#res=StringSpelledByGenomePathProblem(kolekcija_dijelica)
#print(res) # # output : ACCGAAGCT ==> OK
# kad imam veliki dataset - puno charactera, mogu provjeriti rjesenje preko ascii tj. numericke reprezentacije charactera

########################################################################################################################################
################################################### 3 C ##################################################################################
########################################################################################################################################

#kolekcija_dijelica=["AAT", "ATG", "ATG", "ATG", "CAT", "CCA", "GAT", "GCC", "GGA", "GGG", "GTT", "TAA", "TGC", "TGG", "TGT"]

def get_prefix(kmer):
    prefix=kmer[:-1]
    return prefix

def get_sufix(kmer):
    sufix=kmer[1:]
    return sufix

def AdjacencyList(Graph,kolekcija):
    sljedbenici_dict={}
    for i in range(len(kolekcija)):
        sljedbenici=[]
        for j in range(len(kolekcija)):
            if(i!=j):
                if(Graph[i][j]==1):
                    sljedbenici.append(kolekcija[j])
            sljedbenici_dict[kolekcija[i]]=sljedbenici
    #print(sljedbenici_dict)
    return sljedbenici_dict


def Overleap(kolekcija):
    shape=(len(kolekcija),len(kolekcija))
    Graph= numpy.zeros(shape, dtype=int, order='C')
    for i in range(len(kolekcija)):
        for j in range(len(kolekcija)):
            if(i!=j):
                kmer1=kolekcija[i]
                kmer2=kolekcija[j]
                sufix=get_sufix(kmer1)
                prefix = get_prefix(kmer2)
                if(prefix==sufix):
                    Graph[i][j]=1
    return Graph

#graf=Overleap(kolekcija_dijelica)
#print(graf)
#adj_list=AdjacencyList(graf,kolekcija_dijelica)
#print(adj_list)
'''OUTPUT ZA : ["ACCGA","CCGAA","CGAAG","GAAGC","AAGCT"] ==> ACCGAAGCT
[[0 1 0 0 0]
 [0 0 1 0 0]
 [0 0 0 1 0]
 [0 0 0 0 1]
 [0 0 0 0 0]]
{'ACCGA': ['CCGAA'], 'CCGAA': ['CGAAG'], 'CGAAG': ['GAAGC'], 'GAAGC': ['AAGCT'], 'AAGCT': []}''' # ==> OK

########################################################################################################################################
####################################################### 3 D ############################################################################
########################################################################################################################################

# Pokusajmo sloziti usmjerene grafove ==> trazi se najduzi graf, ali se kroz svaki svor smije proci samo jedanput
## odnosno, svaki kmer mogu iskoristiti samo jedan put


# Hammiltonov put - Hamiltonian path Problem :
# input : usmjereni graf
# output : put koji posjecuje svaki cvor u grafu i to samo jednom ( ako takav put postoji )

# Ideja ==> De Brujin
k=4
text="AAGATTCTCTAC"

def DeBrujinGraph(text,k):
    #print(text)
    edges=findKmers(text,k)
    #print(edges)
    nodes=[]
    for edge in edges:
        node=edge[:-1]
        nodes.append(node)
    nodes.append( edges[len(edges)-1][1:]  )
    #print(nodes)
    nodes_dict={}
    unique_nodes = list(set(nodes))
    for unique_node in unique_nodes:
        nodes_dict[unique_node]=[]
    for i in range(len(nodes)-1):
        nodes_dict[nodes[i]].append(  nodes[i+1] )
    #print("*******************")
    #print(nodes_dict)
    return nodes_dict

#deBrujin=DeBrujinGraph(text,k)
#print(deBrujin)
# OUTPUT : {'CTC': ['TCT'], 'CTA': ['TAC'], 'TAC': [], 'AGA': ['GAT'], 'TTC': ['TCT'], 'GAT': ['ATT'], 'ATT': ['TTC'], 'TCT': ['CTC', 'CTA'], 'AAG': ['AGA']} => OK

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



kolekcija_dijelica=["AAT", "ATG", "ATG", "ATG", "CAT", "CCA", "GAT", "GCC", "GGA", "GGG", "GTT", "TAA", "TGC", "TGG", "TGT"]

def get_sufixes(kmer):
    prefixi=[]
    for i in range(1,len(kmer)):
        prefix=kmer[i:]
        prefixi.append(prefix)


    #print(prefixi)
    return prefixi

def get_prefixes(kmer):
    sufixi=[]
    for i in range(1,len(kmer)):
        prefix=kmer[:i]
        sufixi.append(prefix)
    sufixi.reverse()
    #print(sufixi)
    return sufixi
#get_prefixes("ATGTCGT")
#get_sufixes("ATGTCGT")


