import numpy as np
from numpy import matlib as mt
import json
import matplotlib.pyplot as plt
import networkx as nx

def proba():
  #isprobaj na matrici kad je broj vrhova= 3, a broj bridova= 2
  q=[[0,0,1],[0,0,1],[1,1,0]]
  w={}
  for i in range(3):
    w[i]=[]
  for i in range(3):
    for j in range(3):
      if(q[i][j]==1): w[i].append(j)



# create a row vector of given size

def create_matrix(size):
  M = mt.rand(1,size)

  # create a symmetric matrix size * size
  symmA = M.T * M

  A = symmA
  for i in range(size):
    for j in range(size):
      if (symmA[i, j] > 0.5):
        A[i, j] = int(1)
      else:
        A[i, j] = int(0)

  np.savetxt('matrica.txt', A, fmt='%.0f')
  return A

def open_matrix_fromTxt(filename):
  with open(filename, 'r') as f:
    l = [[int(num) for num in line.split(' ')] for line in f]
    return l

def create_adj_list(matrix):
  adjList = {}
  for i in range(size):
    adjList[i] = []
  for i in range(size):
    for j in range(size):
      if (matrix[i][j] == 1): adjList[i].append(j)
  return adjList

def save_dict_as_json(dict,filename):
  with open(filename, 'w') as fp:
    json.dump(dict, fp)

def open_dict_from_json(filename):
  f = open(filename)
  adj_list = json.load(f)
  f.close()
  return adj_list

def show_graph(adjacency_matrix):
    matrix=np.array(adjacency_matrix)
    G = nx.from_numpy_array(matrix)
    nx.draw(G, node_size=5)
    plt.savefig("plot.png")
    plt.show()

def edges_from_adjList(adj_list):
  edges=[]
  for key in adj_list:
    if(len(adj_list[key])>0):
      for el in adj_list[key]:
        k=int(key)
        tup=(k,el)
        atup=(el,k)
        if tup not in edges and atup not in edges:
          edges.append(tup)
  return edges

def degrees(matrix):
    degs=[]
    matrix=np.array(matrix)
    res = list(map(sum, matrix))
    res=list(map(int,res))
    return res

def hist_degrees(degrees):
  plt.hist(degrees)
  plt.savefig("degrees.png")
  plt.show()

size  = 12
matrix=create_matrix(size)
#filename='matrica.txt'
#matrix=open_matrix_fromTxt(filename)
#print(matrix)
#show_graph( matrix)
#adj_list=create_adj_list(matrix)
adj_list_name='adjList.json'
#save_dict_as_json(adj_list,adj_list_name)
adj_list=open_dict_from_json(adj_list_name)
print(adj_list)
edges=edges_from_adjList(adj_list)
print(edges)
m=len(edges) # broj bridova
print(m)
degs=degrees(matrix)
print(degs)
n=size # broj vrhova
k=(2*m)/n #prosjecni stupanj
print(k)
ro=k/(n-1) #gustoca mreze
print(ro)
hist_degrees(degs)



