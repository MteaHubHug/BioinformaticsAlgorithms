# Finding origin of replication problem
# U originu pocinje replikacija DNA, zato je vazno pronaci gdje se nalazi.
# Ova funkcija broji koliko se puta dani Pattern (k-mer nukleotida) nalazi u danom Tekstu (DNA Genomu)
# 1A
# Input : A DNA string Genome
# Output : The location of oriC in Genome
import math


def PatternCount(Text, Pattern):
    count = 0
    for i in range(0, len(Text) - len(Pattern)):        # promatramo svaki znak od pocetka do kraja teksta (duljini teska oduzmemo duljinu patterna
        if (Text[i :  (i + len(Pattern))] == Pattern):  #### jer cemo je u if uvjetu dodati - ne zelimo "index out of range"
            count += 1
    return count

#Text = "ACAACTATGCATACTATCGGGAACTATCCT"
#Pattern= "ACTAT"
#res = PatternCount(Text,Pattern)
#print(res) # res = 3

# Frequent Words Problem - Find the most frequent k-mers in a string
# Ova funkcija nalazi najfrekventnije k-mere u DNA genomu -
# 1B
# Input : A string Text and an integer k
# Output : All most frequent k-mers in Text
def removeDuplicates(lista):
    result = []
    for i in lista:
        if i not in result:
            result.append(i)
    return result

def FrequentWords(Text, k):
    FrequentPatterns=[]
    COUNTS={}
    for i in range (0, len(Text) - k):
        Pattern = Text[i : i + k]
        count_i= PatternCount(Text, Pattern)
        COUNTS[i] = count_i
    maxCount =  max(COUNTS.values()) # build-in funkcija koja trazi maksimum u dictionary
    for i in range(0, len(Text) - k):
        if ( COUNTS[i] ==maxCount ):
            FrequentPatterns.append(  Text[i : i + k])
    FrequentPatterns= removeDuplicates(FrequentPatterns)  # Pretvaranje liste u Set kako bismo se rijesili duplikata
    return FrequentPatterns

#Text= "atcaatgatcaacgtaagcttctaagcATGATCAAGgtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagATGATCAAGagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaATGATCAAGctgctgctcttgatcatcgtttc"
#Text = Text.lower()
#k = 9
#res = FrequentWords(Text, k)
#print(res)    # ['atgatcaag', 'ctcttgatc', 'tcttgatca', 'cttgatcat']

###################################################################################################

# Reverse Complement
# Kada dobijemo moguće k-mere, bilo bi dobro provjeriti koji je njihov inverzni pattern.
# Naime, jedna DNA nit ima svoj smjer, a druga koja joj je komplementarna, kreće se u suprotnom smjeru.
# Algoritam ReverseComplement() nam pomaže naći inverzni pattern.

# 1C

def ReverseComplement(Pattern):
    Reverse=""
    Pattern = Pattern.upper()
    for nukleotid in Pattern:
        if(nukleotid=="A"): Reverse+="T"
        elif(nukleotid=="C"): Reverse+"G"
        elif(nukleotid=="G"): Reverse+="C"
        else : Reverse+="A"
    Reverse = Reverse[::-1]
    return  Reverse

#Pattern = "atgatcaag"
#Reverse = ReverseComplement(Pattern)
#print(Reverse)


#####################################################################
# Clump finding problem
# 1E

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

def ClumpFinding(Genome, k, L, t):
    clums=[]
    k_mers= find_k_mers(Genome,k)
    for i in range(len(Genome) - k +1 ):
        print("k-mer : ", k_mers[i])

    cnts=[]

    for i in range(len(Genome) - k +1):
        cnts.append(1)

    l=0
    for i in range(len(k_mers)-1):
        for j in range(i+1, len(k_mers)):
            if(k_mers[i]==k_mers[j]):
                cnts[l]+=1
        l+=1

    for i in range(len(cnts)):
        print(" cnt :: ", cnts[i])

    for i in range(len(cnts)):
        if(cnts[i]==t): clums.append(k_mers[i])

    return clums

#genom = "CTGCAATGCATGACAAGCCTGCAGTTGCAACCGTAGCATGACGCCCCCCTGAT"
#print("genom :  ", genom);
#k = 4;
#L = 25;
#t = 3;

#klamz = ClumpFinding(genom, k, L, t);

#for clum in klamz:
#    print(clum)    ### TGCA, CCCC

####################################################################

'''
1F - C# 
###Skew bi trebao doseći minimum na poziciji gdje Reverse polunit završava, a Forward polunit počinje. Točno to je lokacija Origina oriC.
using System;
using System.Collections.Generic;
/*
 SKEWi(Genome) as the difference between the total number
of occurrences of G and the total number of occurrences of C
 */

namespace zadF
{
    class Program
    {



        static int Gs(string genom)
        {
            int cnt = 0;
            for (int i = 0; i < genom.Length; i++)
            {
                if (genom[i] == 'G') cnt++;
            }
            return cnt;
        }

        static int Cs(string genom)
        {
            int cnt = 0;
            for (int i = 0; i < genom.Length; i++)
            {
                if (genom[i] == 'C') cnt++;
            }
            return cnt;
        }

        static int skew(string genom)
        {
            int cntG = Gs(genom);
            int cntC = Cs(genom);
            int diff = cntG - cntC;
            return diff;
        }

        static int[] tocke(string genom)
        {
            int[] razlike = new int[genom.Length];
            string komad = "";
            for (int i = 0; i < genom.Length; i++)
            {
                komad += genom[i];
                razlike[i] = skew(komad);
            }

            for(int i=0;i<razlike.Length;i++)
            {
                Console.WriteLine("pozicija : {0} ::: razlika G-C : {1}", i, razlike[i]);
            }


            return razlike;

        }

        static List<int> mins (string genom)
        {
            List<int> positions = new List<int>();
            int[] razlike = tocke(genom);
            int min = razlike[0];
            for(int i=0;i<razlike.Length;i++)
            {
                if (min > razlike[i]) min = razlike[i];
            }
            
            for (int i = 0; i < razlike.Length; i++)
            {
                if(min==razlike[i])
                {
                    positions.Add(i);
                }
            }
            return positions;
        }


        static void Main(string[] args)
        {

            string genom = "CATGGGCATCGGCCATACGCC";
            List<int> positions = mins(genom);
            foreach (int el in positions)
            {
                Console.WriteLine("minimum Skew dijagrama je u tocki : {0}",el);
            }

           
        }
    }
}

'''

####################################################################
# 1 G
# DnaA se zapravo može spojiti i na sitne varijacije DnaA boxova, a ne samo na "savršene" boxove.
# Kada imamo dva k-mera p i q, broj mismatchova između njih se naziva Hamming distance.

def HammingDistance(p,q):
    if(len(p)!=len(q)):
        print("Greska")
        return -1
    HD= sum(c1 != c2 for c1, c2 in zip(p, q))
    return HD

p='AGATGG'
q='AGATAA'
HammingDist=HammingDistance(p,q)
#print(HammingDist)


# 1 H Approximate Pattern Matching Problem
# Find all approximate occurrences of a pattern in a string
#### u C#-u ...

####################################################
# 1K - Frequency Array
# funkcije koje su potrebne u algoritmu:

def SymbolToNumber(symbol):
    num=-1
    symbol=symbol.upper()
    if(symbol=="A"): num=0
    elif(symbol=="C"): num=1
    elif(symbol=="G"): num=2
    elif(symbol=="T"): num=3
    return num

def LastSymbol(Pattern):
    patt= Pattern[-1:]
    return patt

def Prefix(Pattern):
    patt=Pattern[:-1]
    return  patt

def PatternToNumber(Pattern):
    if(len(Pattern)==0) : return 0
    symb=LastSymbol(Pattern)
    prefix= Prefix(Pattern)
    return 4 * PatternToNumber(prefix) + SymbolToNumber(symb)

### print(PatternToNumber("GT")) # =11

def NumberToSymbol(number):
    symb="X"
    if(number==0): symb="A"
    elif(number==1): symb="C"
    elif(number==2): symb="G"
    elif(number==3): symb="T"
    return symb

def NumberToPattern(index,k):
    if k==1 : return NumberToSymbol(index)

    prefixIndex= int(index//4) #kvocijent
    r = int(index%4)
    symbol = NumberToSymbol(r)
    PrefixPattern = NumberToPattern(prefixIndex, k-1)
    return PrefixPattern + symbol

###print(NumberToPattern(11,2)) # = GT

###### PRIMJER :
text= "AAGCAAAGGTGGG"
k=2
k_mers=find_k_mers(text,k)
print(k_mers)  # 2-meri : ['AA', 'AG', 'GC', 'CA', 'AA', 'AA', 'AG', 'GG', 'GT', 'TG', 'GG', 'GG']
indexi=[]
for kmer in k_mers:
    idx=PatternToNumber(kmer)
    kmer_idx_tuple= (kmer,idx)
    indexi.append(kmer_idx_tuple)

print(indexi)

def SortTuples(indexi):
    sorted_indices=sorted(indexi, key=lambda x: x[1])
    return sorted_indices

sorted_indices=SortTuples(indexi)
print(sorted_indices)

def FindingFrequentWordsBySorting(Text,k):
    FrequentPatterns=[]
    INDEX=[]
    COUNT=[]
    for i in range(len(Text)-k+1):
        Pattern=Text[i: (i+k)]
        INDEX.append(PatternToNumber(Pattern) )
        COUNT.append(1)
    SORTED_INDEX=INDEX
    SORTED_INDEX.sort()
    for i in range(1,len(Text)-k+1):
        if(SORTED_INDEX[i]==SORTED_INDEX[i-1]): # ovo je glavna fora di se iskoristava sortiranost
            COUNT[i]=COUNT[i-1] +1
    maxCount=max(COUNT)
    for i in range(0,len(Text)-k+1):
        if(COUNT[i]==maxCount):
            Pattern=NumberToPattern(SORTED_INDEX[i],k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns

frekventni_patterni=FindingFrequentWordsBySorting(text,k)
print(frekventni_patterni)  # ['AA', 'GG']

###############################
# CLUMP FINDING  poboljsani -  to ima u knjizi

################################
# NEIGHBOURS

def ImmediateNeighbours(Pattern):
    Pattern=Pattern.upper()
    nukleotids=["A","B","C","D"]
    Neighbourhood= [Pattern]
    for i in range(1,len(Pattern)):
        symbol=Pattern[i]
        for nukleotid in nukleotids:
            if(nukleotid!=symbol):
                Neighbour= Pattern[0:i] + nukleotid.lower() + Pattern[i+1: len(Pattern)]
                # mijenjamo simbol za nukleotid da bismo generirali susjedstvo
                Neighbourhood.append(Neighbour)
    return Neighbourhood

print("*********** Neposredno Susjedstvo :  ****************")
print(ImmediateNeighbours(text)) # ['AAGCAAAGGTGGG', 'AbGCAAAGGTGGG', 'AcGCAAAGGTGGG', 'AdGCAAAGGTGGG', 'AAaCAAAGGTGGG', 'AAbCAAAGGTGGG', 'AAcCAAAGGTGGG', 'AAdCAAAGGTGGG', 'AAGaAAAGGTGGG', 'AAGbAAAGGTGGG', 'AAGdAAAGGTGGG', 'AAGCbAAGGTGGG', 'AAGCcAAGGTGGG', 'AAGCdAAGGTGGG', 'AAGCAbAGGTGGG', 'AAGCAcAGGTGGG', 'AAGCAdAGGTGGG', 'AAGCAAbGGTGGG', 'AAGCAAcGGTGGG', 'AAGCAAdGGTGGG', 'AAGCAAAaGTGGG', 'AAGCAAAbGTGGG', 'AAGCAAAcGTGGG', 'AAGCAAAdGTGGG', 'AAGCAAAGaTGGG', 'AAGCAAAGbTGGG', 'AAGCAAAGcTGGG', 'AAGCAAAGdTGGG', 'AAGCAAAGGaGGG', 'AAGCAAAGGbGGG', 'AAGCAAAGGcGGG', 'AAGCAAAGGdGGG', 'AAGCAAAGGTaGG', 'AAGCAAAGGTbGG', 'AAGCAAAGGTcGG', 'AAGCAAAGGTdGG', 'AAGCAAAGGTGaG', 'AAGCAAAGGTGbG', 'AAGCAAAGGTGcG', 'AAGCAAAGGTGdG', 'AAGCAAAGGTGGa', 'AAGCAAAGGTGGb', 'AAGCAAAGGTGGc', 'AAGCAAAGGTGGd']

print("*********** Susjedstvo :  ****************")
def Suffix(Pattern):
    patt=Pattern[1:]
    return  patt

def FirstSymbol(Pattern):
    patt= Pattern[0]
    return patt
#### 1 N
def Neighbors(Pattern,d):
    if(d==0): return [Pattern]
    if(len(Pattern)==1): return["A","C","G","T"]
    nukleotids=["A","C","G","T"]
    Pattern=Pattern.upper()
    Neighborhood=[]
    SuffixNeighbors= Neighbors(Suffix(Pattern),d)
    for Text in SuffixNeighbors:
        if(HammingDistance(Suffix(Pattern),Text)<d):
            for nukleotid in nukleotids:
                Neighborhood.append( nukleotid+Text  )
        else:
            Neighborhood.append(FirstSymbol(Pattern)+ Text)
    return Neighborhood


#Susjedstvo=Neighbors(text,2)
#for susjed in Susjedstvo:
#    print(susjed)

# Na ovaj način možemo generirati Susjedstvo svih k-mera Hammingove distance najviše d od Patterna.
# Modificirajmo Funkciju Neighbours tako da generira Susjedstvo svih k-mera koji su udaljeni od Patterna točno za d ::

'''def IterativeNeighbors(Pattern,d):
    Neighbourhood=set()
    Neighbourhood.add(Pattern)
    for j in range(1,d):
        for Pattern_ in Neighbourhood:
            Neighbourhood.add(ImmediateNeighbours(Pattern_))
    return Neighbourhood'''


# Pomocu ovog algoritma mozemo naci Frekventne rijeci u stringu s mismatchovima pomocu sortiranja ==> poboljsani i brzi algoritam
def FindingFrequentWordsWithMismatchesBySorting(Text,k,d):
    FrequentPatterns=set()
    Neighborhoods=[]
    Index=[]
    Count=[]
    for i in range(len(Text)-k):  # ovdje mozda treba dodati -1 ==> len(Text)-k-1
        Neighborhoods.append(Neighbors(Text[i : (i+k)],d))
    NeighborhoodArray=[]
    for l in Neighborhoods:
        for string in l:
            NeighborhoodArray.append(string)

    for i in range(len(Neighborhoods)-1): # Ovdje mozda treba dodati -1
        Pattern= NeighborhoodArray[i]
        Index.append(  PatternToNumber(Pattern) )
        Count.append(1)
    SortedIndex=Index
    SortedIndex.sort()
    for i in range(len(Neighborhoods)-1-1):
        if(SortedIndex[i]==SortedIndex[i+1]):
            Count[i+1]=Count[i]+1
    maxCount=max(Count)
    for i in range(len(Neighborhoods)-1): # Ovdje mozda treba dodati -1
        if(Count[i]==maxCount):
            Pattern=NumberToPattern(SortedIndex[i],k)
            FrequentPatterns.add(Pattern)
    return FrequentPatterns

