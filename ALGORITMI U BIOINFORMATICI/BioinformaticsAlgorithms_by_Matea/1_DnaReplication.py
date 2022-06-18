# Finding origin of replication problem
# U originu pocinje replikacija DNA, zato je vazno pronaci gdje se nalazi.
# Ova funkcija broji koliko se puta dani Pattern (k-mer nukleotida) nalazi u danom Tekstu (DNA Genomu)
# 1A
# Input : A DNA string Genome
# Output : The location of oriC in Genome

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
#### u C#-u ...





