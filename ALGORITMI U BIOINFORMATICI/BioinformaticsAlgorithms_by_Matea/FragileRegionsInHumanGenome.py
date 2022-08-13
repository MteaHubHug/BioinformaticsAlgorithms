##############################################################################################################################################################
####################################################################  6 A #################################################################################
##############################################################################################################################################################

# Implementation
def Reversing(seq):
    seq.reverse()
    seq = list(map(lambda x: -x, seq))
    return seq


def GreedySorting(Permutation):
    change = []
    ARD = 0  # ApproxReversalDistance (ARD)

    for i in range(len(Permutation)):
        if Permutation[i] != i+1:
            if i+1 in Permutation:
                id_i = Permutation.index(abs(i+1))
            else:
                id_i = Permutation.index(-abs(i+1))
            Permutation[i: id_i + 1] = Reversing(Permutation[i: id_i + 1])
            ARD += 1
            change.append(Permutation[:])
        if Permutation[i] == -(i+1):
            Permutation[i] = i+1
            ARD += 1 # u knjizi se trazi return ARD, a u Rosalindu sortirana permutacija kao output
            change.append(Permutation[:])  # Permutation[:] mogu koristiti da bih napravila kopiju
    return change

P=[-3, +4, +1, +5, -2]
#rez=GreedySorting(P)
#print(rez)
# [[-1, -4, 3, 5, -2], [1, -4, 3, 5, -2], [1, 2, -5, -3, 4], [1, 2, 3, 5, 4], [1, 2, 3, -4, -5], [1, 2, 3, 4, -5], [1, 2, 3, 4, 5]]

##############################################################################################################################################################
####################################################################  6 B #################################################################################
##############################################################################################################################################################
# BREAKPOINTS
# Pronadji broj breakpoints od P
P =[+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14]

def breakpoint_count(permutation):
    a=permutation+[len(permutation)+1] # ovo pise u knjizi da se dodaje n+1 na kraj i 0 na pocetak permutacije da bi se moglo ovako checkirati breakpoints
    b= [0]+permutation
    c=list(map(lambda x,y: x - y != 1,a ,b))
    #print(c)  # [True, False, False, True, True, False, False, True, False, True, True, True, True, False, False]
    return sum(c)

#brojBreakpoints=breakpoint_count(P)
#print(brojBreakpoints) ### =8