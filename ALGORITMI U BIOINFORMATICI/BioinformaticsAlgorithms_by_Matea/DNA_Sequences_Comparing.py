from reading_tasks_text import *


'''
Recimo da imamo 3 duga niza, segmenta za koja znamo da otprilike rade iste stvari
... npr kao automobili, svi imaju 4 kotaca, volan itd, i svi nam koriste za voznju,
ali nije svaki automobil jednak drugome
==> Dakle, dani segmenti generiraju neke aminokiseline, ali se razlikuju u tome koju aminokiselinu
generiraju
=> Ako oni rade tu slicnu stvar, a rade jer svi generiraju aminokiseline,
kako mi to mozemo prepoznati ako je zapisano samo tim slovima (oznake za aminkiseline)
=> Nema nista isto u njima, ali postoje A-domene (dijelovi koji dodaju aminkiselinu na odredjeni peptid)
=> Na tri mjesta se A-domene podudaraju

Pokusajmo narediti nase stringove da se malo vise podudaraju, a ne samo na tri mjesta
...ako imamo ovakvu situaciju :
YAFDLGYTCMFPV*L*LG*G*GELHIVQKETYTAPDEIAHYIKEHGITYIKLTPSLFHTIVNTASFAFDANFESLRLIVLGGEKIIPIDVIAFRKMYGHTEFINHYGPTEATIGA
AFDVSAGDFARAL*L*TG*G*QLIVCPNEVKMDPASLYAIIKKYDITIFEATPALVIPLMEYIYEQKLDISQLQILIVGSDSCSMEDFKTLVSRFGSTIRIVNSYGVTEACIDS
IAFDASSWEIYAP*L*LN*G*GTVVCIDYYTTIDIKALEAVFKQHHIRGAMLPPALLKQCLVSAPTMISSLEILFAAGDRLSSQDAILARRAVGSGVYNAYGPTENTVLS

mozemo drugi string pomaknut za 1 udesno da dobijemo vise podudaranja

YAFDLGYTCMFPV*LL*G*GG*ELHIVQKETYTAPDEIAHYI*K*EHG*I*TYIKLT*P*S*L*FHTIVNTASFAFDANFESLRLIVLGGEKIIPIDVIAFRKMYGHTEFINHYGPTEATIGA
-AFDVSAGDFARA*LL*T*GG*QLIVCPNEVKMDPASLYAII*K*KYD*I*TIFEAT*P*A*L*VIPLMEYIYEQKLDISQLQILIVGSDSCSMEDFKTLVSRFGSTIRIVNSYGVTEACIDS
IAFDASSWEIYAP*LL*N*GG*TVVCIDYYTTIDIKALEAVF*K*QHH*I*RGAMLP*P*A*L*LKQCLVSAPTMISSLEILFAAGDRLSSQDAILARRAVGSGVYNAYGPTENTVLS

==> naredimo to jos bolje, umetnimo neke praznine u treci string
( postupak se zove *alignanje(poravnanje A-domena* )

YAFDLGYTCMFPV*LL*G*GG*ELHIVQKETYTAPDEIAHYI*K*EHG*I*TYIKLT*P*S*L*FHTIVNTASFAFDANFES*L*RLIVLGGEKIIPI*D*VIAFRKMY*G*HTEFINHYGPTEATIGA
-AFDVSAGDFARA*LL*T*GG*QLIVCPNEVKMDPASLYAII*K*KYD*I*TIFEAT*P*A*L*VIPLMEYI-YEQKLDISQ*L*QILIVGSDSCSME*D*FKTLVSRF*G*STIRIVNSYGVTEACIDS
IAFDASSWEIYAP*LL*N*GG*TVVCIDYYTTIDIKALEAVF*K*QHH*I*RGAMLP*P*A*L*LKQCLVSA----PTMISS*L*EILFAAGDRLSSQ*D*AILARRAV*G*SGVYNAYGPTENTVLS

=> Sada nam se cini da se malo bolje podudaraju
=> Dovodi nas do pretpostavke da su svi ti stringovi bili otprilike jednaki,
ali se mutacijom/ prepisivanjem mozda izgubilo tih nekoliko slova koja su bila nebitna
=> kada bismo jos namjestali i umetali praznine, dobili bismo 19 stupaca koji se podudaraju
=> sada mozemo pretpostaviti da su tih 19 stupaca opisuje neku funkciju, a ostali stupci opisuju aminokiselinu

=> 19 stupaca koji se podudaraju : "ja sam A-domena"
=> ostali stupci : "sad kad znam da sam A-domena, ja specificno gradim odredjeni tip aminokiseline"

=> tih 19 stupaca koji se podudaraju se zovu *KONZERVIRANA BAZA* - svaka domena to ima
=> mozemo pretpostaviti da ostala slova, koja nisu konz. baza, kodiraju nekako ostale aminokiseline
=> kroz odredjene eksperimente se zakljucuje da  LTKVLGHIG, VGEIVGSID, and AWMFAAAVL kodiraju aminokiseline
=> oznaceni su ovako : #slovo# ::
Y*AFD*#L#GY#T#CMFPV*LL*G*GG*ELHIVQKETYTAPDEIAHYI*K*EHG*I*TYI#K#LT*P*S*L*FHTIVNTASFAFDANFES*L*RLIVLGGEKIIPI*D*VIAFRKMY*G*HTEFINHYGPTEATIGA
-*AFD*#V#SA#G#DFARA*LL*T*GG*QLIVCPNEVKMDPASLYAII*K*KYD*I*TIF#E#AT*P*A*L*VIPLMEYI-YEQKLDISQ*L*QILIVGSDSCSME*D*FKTLVSRF*G*STIRIVNSYGVTEACIDS
I*AFD*#A#SS#W#EIYAP*LL*N*GG*TVVCIDYYTTIDIKALEAVF*K*QHH*I*RGA#M#LP*P*A*L*LKQCLVSA----PTMISS*L*EILFAAGDRLSSQ*D*AILARRAV*G*SGVYNAYGPTENTVLS
...
ostala slova sluze samo formaciji, obliku i nema kodirajuce svojstvo (non-coding dio)

LTKVGHIG -> Asp
VGEIGSID -> Orn
AWMFAAVL -> Val

=> Vazno je uociti da bez konstruiranja konzervirane jezgre, istrazivaci ne bi mogli
zakljuciti sto je non-ribosomal code
=> jer se 24 aminokiseline ne slazu ako ne umetnemo praznine i ne naredimo to da vidimo podudaranje
=> kako mozemo alignati / poravnati stringove sto bolje?

~~~~~~~~~~~~~~~~~~~~~~INTRODUCTION TO SEQUENCE ALIGNMENT~~~~~~~~~~~~
Probajmo napraviti sto vise matchova poravnavanjem ova dva stringa :
ATGCATGC
TGCATGCA

==>
ATGCATGC-
-TGCATGCA

............
probaj rijesiti isti zadatak za ova dva stringa :  ATGCTTA, TGCATTAA :
ATGC-TTA-
-TGCATTAA

Def. ALIGNMENT od sequences v i w :
matrica od 2 reda t.d. prvi red sadrzi simbole iz v ,
a drugi red sadrzi simbole iz w
i dodani su "space" simboli
=> vazno! : ne dopustaju se dvije praznine jedna iznad(ispod) druge
 npr:
 A T - G T T A T A
 A T C G T - C - C

 slucajevi :
  -match - stupac u kojem se gornji i donji simbol podudaraju
  -mismatches - stupac u kojem se gornji i donji simbol ne podudaraju
  -indels - stupac koji sadrzi "space" symbol -
  -insertion - stupac koji ima "space" simbol u gornem retku
  -deletion - stupac koji ima "space" simbol u donjem retku

=> ono sto se podudara definira zajednicki podniz dva stringa

=> poravnanje dva stringa koje maksimizira broj matcheva odgovara
 najduljem zajednickom podnizu tih stringova
 (to je tezak problem)
 => jedan string moze imat vise od jednog odgovarajuceg poniza

 PROBLEM :
 LONGEST COMMON SUBSEQUENCE PROBLEM :
 pronadji najduzi zajednicki podniz dva stringa
 input : dva stringa
 output : najduzi zajednicki podniz ta dva stringa

 Ako se limitiramo na dvije A-domene (npr. za Asp i Orn),
 sukladno s onih 19 koje smo nasli dodavajuci "-",
 mozemo naci jos 10 matcheva ==> dobivamo najduzi zajednicki podniz duljine 29


 za
 Mateja i
 Matea =>

 Mateja
 Mate-a ==> 4 matcha ==> Mata je najduzi zajednicki niz ta dva stringa


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ THE MANHATTAN TOURIST PROBLEM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Koja je najbolja sightseeing strategija?

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DINAMICKO PROGRAMIRANJE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Cuvamo medjurezultate i koristimo rekurzije za racunanje kako treba vratiti novac
Umjesto da puno puta racunamo jednu te isto vrijednost, sacuvamo taj rezultat i koristimo ga u daljnjem racunanju
'''

#####################################################################################################################################################
############################################################## 5 A ##################################################################################
#####################################################################################################################################################
import numpy as np
# slozenost :: O( money * |Coins| )

def DPChange(money, Coins):
    MinNumCoins={}
    MinNumCoins[0]=0
    for m in range(1,money+1):
        MinNumCoins[m]=np.Infinity
        for i in range(len(Coins)):
            if(m>=Coins[i]):
                if(MinNumCoins[m - Coins[i]]+1 < MinNumCoins[m]):
                    MinNumCoins[m]=MinNumCoins[m - Coins[i]] +1
    return MinNumCoins[money]

money=40
coins=[1,5,10,20,25,50]
#minNumCoins=DPChange(money,coins)
#print(minNumCoins) #### => rez =2 ==> OK
'''
infile= open_file("rosalind_ba5a.txt")
money=int(infile[0])
coins=infile[1].split(",")
coins = [int(i) for i in coins]
minNumCoins=DPChange(money,coins)
print(minNumCoins)


'''
#####################################################################################################################################################
############################################################## 5 B ##################################################################################
#####################################################################################################################################################
'''
MANHATTAN TOURIST PROBLEM REVISITER
sljedeci pseudokod izracunava udaljenost najduljeg puta do nekog cvira (i,j)
u nasoj pravokutnoj mrezi ;
 zasniva se na ideji da je jedini nacin do nodea (i,j) je da se moramo spustiti dolje
 za (i-1) puta, i moramo se pomaknuti desno za (j-1) puta
'''
n=4
m=4
Down=[
[1, 0, 2, 4, 3],
[4, 6, 5, 2, 1],
[4, 4, 5, 2, 1],
[5, 6, 8, 5, 3]]
Right=[
[3, 2, 4, 0],
[3, 2, 4, 2],
[0, 7, 3, 3],
[3, 3, 0, 2],
[1, 3, 2, 2],
]


def ManhattanTourist(n,m,Down,Right):
    s=np.zeros((n+1, m+1))
    for i in range(1,n+1):
        s[i][0]= s[i-1][0] + Down[i-1][0]
    for j in range(1,m+1):
        s[0][j]= s[0][j-1] + Right[0][j-1]
    for i in range(1,n+1):
        for j in range(1,m+1):
            s[i][j]= max( s[i-1][j] + Down[i-1][j] , s[i][j-1] + Right[i][j-1])
    return s[n][m]

#put=ManhattanTourist(n,m,Down,Right)
#print(put)  # 34.0 ==> ok


'''

infile= open_file("rosalind_ba5b.txt")
nums=infile[0].split(" ")
n=int(nums[0])
m=int(nums[1])
flag=0
Down=[]
Right=[]
for l in infile[1:]:
    if(l!="-" and flag==0):
        row=l.split(" ")
        row = [int(i) for i in row]
        Down.append(row)
    elif(l=="-"):
        flag=1
    elif(l!="-" and flag==1):
        row = l.split(" ")
        row = [int(i) for i in row]
        Right.append(row)

put=ManhattanTourist(n,m,Down,Right)
print(put)

'''

#####################################################################################################################################################
############################################################## 5 E ##################################################################################
#####################################################################################################################################################
'''
Sto je krivo s LCS (Longest Commont Sequence) scoring modelom?
Moze se dogoditi da pretjeramo s alignmentiranjem / poravnavanjem tako sto dodajemo previse indela ( praznine) i tako se gubi originalna bioloska struktura 
Kod LCS scoring modela, nista nas ne sprjecava da dodamo ogroman broj indela kako bismo izgradili sto bolji LCS model
=> U originalnoj alignment igri smo nagradjivali bodovima sobole koji se matchaju, ali trebali bismo uvesti i neku vrstu penala za neke indele i mismatchove (nepodudarnosti)

~~~SCORING MATRICE~~~~
Da bismo generalizirali alignment scoring model, jos uvijek nagradjujemo match s +1,
ali takodjer penaliziramo mismatch s konstantom  µ (mismatch penalty) 
te penaliziramo indele s konstantom simga (indel penalty) 
==>
ALIGNMENT SCORE = #matches - µ* #mismatches - sigma * #indels
Primjer. Za µ=1, sigma=2 , alignment score= -4
A   T  -  G  T  T  A  T  A
A   T  C  G  T  -  C  -  C
---------------------------
+1 +1 -2 +1 +1 -2 -1 -2 -1 = -4

Biolozi su kasnije popravili cost funkciju tako da dopusta da su neke mutacije
vjerojatnije od drugih
=> Zbog toga se mismatch penali i indel penali razlikuju za razlicite simbole
SCORING MATRIX -  (k+1)×(k+1) matrica  (+1 zato sto ukljucujemo u indel)
koja sadrzi nagrade i penale 
Pr : za k=4 : sigma=s 
   A  C  G  T  - 
A +1 -µ -µ -µ -s
C -µ +1 -µ -µ -s
G -µ -µ +1 -µ -s
T -µ -µ -µ +1 -s
- -s -s -s -s
'''
# GLOBAL ALIGNMENT PROBLEM
# Nadji highest-scoring alignment dva stringa gdje su bodovi definirani Scoring matricom
# input : dva stringa i scoring matrica "Score"
# output : poravnanje dva stringa ciji je alignment score maksimiziran (prema matrici "Score")
       # uzevsi u obzir sva moguca poravnanja ta dva stringa

import pandas as pd
blossom_table = pd.read_csv('blossom_table.txt', sep="\t", header=None)
blossom_table=np.array(blossom_table)
#print(blossom_table)



#####################################################
def global_alignment(seq1, seq2, score_matrix, sig):
    len1, len2 = len(seq1), len(seq2)
    s = [[0] * (len2 + 1) for i in range(len1 + 1)]
    backtrack = [[0] * (len2 + 1) for i in range(len1 + 1)]
    for i in range(1, len1 + 1):
        s[i][0] = - i * sig
    for j in range(1, len2 + 1):
        s[0][j] = - j * sig
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            score_list = [s[i - 1][j] - sig, s[i][j - 1] - sig, s[i - 1][j - 1] + score_matrix[seq1[i - 1], seq2[j - 1]]]
            s[i][j] = max(score_list)
            backtrack[i][j] = score_list.index(s[i][j])
    indel_insert = lambda seq, i: seq[:i] + '-' + seq[i:]
    align1, align2 = seq1, seq2
    a, b = len1, len2
    max_score = str(s[a][b])
    while a * b != 0:
        if backtrack[a][b] == 0:
            a -= 1
            align2 = indel_insert(align2, b)
        elif backtrack[a][b] == 1:
            b -= 1
            align1 = indel_insert(align1, a)
        else:
            a -= 1
            b -= 1
    for i in range(a):
        align2 = indel_insert(align2, 0)
    for j in range(b):
        align1 = indel_insert(align1, 0)
    return max_score, align1, align2


if __name__ == '__main__':
    path5e= r"C:\Users\Matea\Desktop\ALGORITMI U BIOINFORMATICI\BioinformaticsAlgorithms_by_Matea\tekstovi_zadataka\rosalind_ba5e.txt"
    with open(path5e) as f:
        seq1 = f.readline().strip()
        seq2 = f.readline().strip()
    with open('BLOSUMM64.txt') as f1:
        lines = [line.strip().split() for line in f1.readlines()]
        score_matrix = {(i[0], i[1]): int(i[2]) for i in lines}
    penalty = 5
    #alignment = '\n'.join(global_alignment(seq1, seq2, score_matrix, penalty))
    #print(alignment)
# odkomentiraj print i alignment za output
'''output : 
8
PLEASANTLY
-MEA--N-LY
'''

