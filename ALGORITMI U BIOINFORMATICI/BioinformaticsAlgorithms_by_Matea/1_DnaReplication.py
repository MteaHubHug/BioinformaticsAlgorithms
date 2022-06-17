# Finding origin of replication problem
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