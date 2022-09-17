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
'''
Cycle To Chromosome Problem (Ciklicki graf => Kromosom) ovo je inverzna funkcija prethodne
input : sekvenca "Nodes" koja se sastoji od integera izmedju 1 i 2n
output : kromosom "Chromosome" koji sadrzi n synteny blokova sto su rezultat funkcije CycleToChromosome(Nodes)
'''
# ovo samo prema pseudokodu iz knjige
def CycleToChromosome(Nodes):
    Chromosome = []
    iteracija = iter(Nodes)
    for i in iteracija:
        a, b = (i, next(iteracija))
        if a < b:
            Chromosome.append(b // 2)
        else:
            Chromosome.append(-a // 2)
    return Chromosome

sekvencaNodes=[1, 2, 4, 3, 6, 5, 7, 8]
#rez=CycleToChromosome(sekvencaNodes)
#print(rez)  # [1, -2, -3, 4] -- OK

##############################################################################################################################################################
####################################################################  6 H #################################################################################
##############################################################################################################################################################

'''
Colored Edges Problem (Problem obojanih bridova) - pronadji Obojene Bridove u genomu
input: genom P
output : kolekcija obojenih bridova u grafu genoma P u obliku (x,y)
'''
def ColoredEdges(P):
    Edges = []
    for Chromosome in P:
        Nodes = ChromosomeToCycle(Chromosome)
        for j in range(len(Chromosome)):
            edge=(Nodes[2*j-1],Nodes[2*j])
            Edges.append(edge)
    Edges.sort()
    return Edges

P=[[+1, -2, -3],[+4, +5, -6]]
#rez=ColoredEdges(P)
#print(rez)  # (2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)  =>OK

'''def ColoredEdges(P):
    Edges = []
    for Chromosome in P:
        Nodes = ChromosomeToCycle(Chromosome)
        it = iter(Nodes[1:]+[Nodes[0]])
        for i in it:
            Edges.append((i,next(it)))
    return Edges'''

##############################################################################################################################################################
####################################################################  6 I #################################################################################
##############################################################################################################################################################
'''
Graph To Genome Problem
input: Obojani bridovi grafa genoma
output : genom koji odgovara grafu genoma
'''

def split_edges_to_chromosomes(edges):
    final_pairs = []
    for i, pair in enumerate(edges):
        if pair[0] > pair[1]:
            final_pairs.append(i)

    chromosomes = []
    previous = 0
    for end in final_pairs:
        chromosomes.append(edges[previous : (end + 1)])
        previous = end + 1
    return chromosomes

def colored_edges_to_cycle(edges):
    cycle = []
    for i in range(len(edges) - 1):
        cycle.extend([edges[i][1], edges[i + 1][0]])
    final = [edges[-1][1], edges[0][0]]
    final.extend(cycle)
    return final

def GraphToGenome(GenomeGraph):
    P=[]
    chromosome_edges = split_edges_to_chromosomes(GenomeGraph)
    for edge in chromosome_edges:
        cycle=colored_edges_to_cycle(edge)
        Chromosome=CycleToChromosome(cycle)
        P.append(Chromosome)
    return P

grafGenoma=[(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]
#rez=GraphToGenome(grafGenoma)  # (+1 -2 -3)(-4 +5 -6)
#print(rez)


##############################################################################################################################################################
####################################################################  6 J #################################################################################
##############################################################################################################################################################

'''
2-Break On Genome Graph Problem
input: obojani bridovi grafa genoma "GenomeGraph" - s indeksima i  i'  j  j'
output: obojani bridovi grafa genoma koji su rezultat funkcije 2-break()

Sample Dataset
(2, 4), (3, 8), (7, 5), (6, 1)
1, 6, 3, 8
'''

def get2BreakOnGenomeGraph(GenomeGraph,i0,i1,j0,j1):
    def eq(x,y):
        u,v=x
        w,z=y
        return (u ==w and v==z) or (w==v and u==z)
    removed = [x for x in GenomeGraph if not eq(x,(i0,i1)) and not eq(x,(j0,j1))]
    return removed + [(i0,j0)] + [(i1,j1)]

grafGranoma=[(2, 4), (3, 8), (7, 5), (6, 1)]
i=1
i_=6
j=3
j_=8
rez=get2BreakOnGenomeGraph(grafGranoma,i,i_,j,j_)
print(rez)  #(2, 4), (3, 1), (7, 5), (6, 8)