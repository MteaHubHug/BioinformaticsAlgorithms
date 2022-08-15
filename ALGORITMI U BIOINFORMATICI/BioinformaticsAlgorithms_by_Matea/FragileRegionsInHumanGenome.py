##############################################################################################################################################################
####################################################################  6 A #################################################################################
##############################################################################################################################################################

# Implementation
import numpy


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

##############################################################################################################################################################
####################################################################  6 C #################################################################################
##############################################################################################################################################################

from collections import defaultdict


def Two_Break_Distance(P, Q):
    '''Returns the 2-Break Distance of Circular Chromosomes P and Q.'''

    # Construct the break point graph of P and Q.
    graph = defaultdict(list)
    for perm_cycle in P+Q:
        n = len(perm_cycle)
        for i in range(n):
            # Add the edge between consecutive items (both orders since the breakpoint graph is undirected).
            # Note: Modulo n in the higher index for the edge between the last and first elements.
            graph[perm_cycle[i]].append(-1*perm_cycle[(i+1) % n])
            graph[-1*perm_cycle[(i+1) % n]].append(perm_cycle[i])

    # Traverse the breakpoint graph to get the number of connected components.
    component_count = 0
    remaining = set(graph.keys())
    while remaining:
        component_count += 1
        queue = {remaining.pop()}  # Undirected graph, so we can choose a remaining node arbitrarily.
        while queue:
            # Select an element from the queue and get its remaining children.
            current = queue.pop()
            new_nodes = {node for node in graph[current] if node in remaining}
            # Add the new nodes to the queue, remove them from the remaining nodes.
            queue |= new_nodes
            remaining -= new_nodes

    # Theorem: d(P,Q) = blocks(P,Q) - cycles(P,Q)
    return sum(map(len,P)) - component_count


'''
(+1 +2 +3 +4 +5 +6)
(+1 -3 -6 -5)(+2 -4)
'''

P=[[+1, +2, +3, +4, +5, +6]]
Q=[[+1, -3, -6, -5],[+2, -4]]
#rez=Two_Break_Distance(P,Q)
#print(rez) # rez=3 ==> OK


##############################################################################################################################################################
####################################################################  6 F #################################################################################
##############################################################################################################################################################
'''
Chromosome To Cycle Problem (Od kromosoma do ciklusa tj. ciklicki graf)
input : kromosom "Chromosome" koji sadrzi synteny blokove
output : sekvenca "Nodes" integera izmedju 1 i 2n 
'''


def ChromosomeToCycle(Chromosome):
    Nodes = []
    for i in Chromosome:
        if i > 0:
            Nodes.append(2 * i - 1)
            Nodes.append(2 * i)
        else:
            Nodes.append(-2 * i)
            Nodes.append(-2 * i - 1)

    return Nodes


kromosom=[+1, -2, -3, +4]
#rez=ChromosomeToCycle(kromosom)
#print(rez) # [1, 2, 4, 3, 6, 5, 7, 8]  ok

##############################################################################################################################################################
####################################################################  6 G #################################################################################
##############################################################################################################################################################