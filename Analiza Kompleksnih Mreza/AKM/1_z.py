'''
Zadatak
1. Konstruirajte neki neusmjereni netežinski graf sa 10 do 20 vrhova i reprezentirajte
ga matricom susjedstva (ne put, ciklus, potpuni itd,...) u .txt file.
2. Prikažite ga grafički u nekom alatu.
3. Ispišite u novi .txt file listu susjeda za taj graf, listu bridova i niz stupnjeva te broj
bridova, prosječni stupanj i gustoću.
4. Nacrtajte histogram za distribuciju stupnjeva.
5. Uklonite neki od vrhova grafa.
6. Ispišite matricu susjedstva za novi graf i spremite u novi file.
7. Konstruirajte usmjereni graf sa oko 10 vrhova i riješite za njega zadatke 2 do 6.
'''


import numpy as np
from numpy import matlib as mt
import json
import matplotlib.pyplot as plt
import networkx as nx
import pickle
import pprint

def proba():
  #isprobaj na matrici kad je broj vrhova= 3, a broj bridova= 2
  q=[[0,0,1],[0,0,1],[1,1,0]]
  w={}
  for i in range(3):
    w[i]=[]
  for i in range(3):
    for j in range(3):
      if(q[i][j]==1): w[i].append(j)

class Matrix(object):

    # Initialize the matrix
    def __init__(self, size):
      M = mt.rand(1, size)

      # create a symmetric matrix size * size
      symmA = M.T * M

      self.A = symmA
      for i in range(size):
        for j in range(size):
          if (symmA[i, j] > 0.5):
            self.A[i, j] = int(1)
          else:
            self.A[i, j] = int(0)

    # Add edges
    def add_edge(self, v1, v2):
      if v1 == v2:
        print("Same vertex %d and %d" % (v1, v2))
      self.A[v1,v2] = 1
      self.A[v2,v1] = 1

    # Remove edges
    def remove_edge(self, v1, v2):
      if self.A[v1,v2] == 0:
        print("No edge between %d and %d" % (v1, v2))
        return
      self.A[v1,v2] = 0
      self.A[v2,v1] = 0


    def save_to_file(self):
      np.savetxt('matrica.txt', self.A, fmt='%.0f')
    def save_to_file2(self):
      np.savetxt('matrica2.txt', self.A, fmt='%.0f')


    # Print the matrix
    def print_matrix(self):
      for row in self.A:
        for val in row:
          print('{:4}'.format(val)),



# create a row vector of given size

def create_matrix1(size):
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

def save_info(matrix,infoname):
  adj_list = create_adj_list(matrix)
  #adj_list_name = 'adjList.json'
  # save_dict_as_json(adj_list,adj_list_name)
  # adj_list=open_dict_from_json(adj_list_name)
  edges = edges_from_adjList(adj_list) # lista bridova
  m = len(edges)  # broj bridova
  degs = degrees(matrix)
  n = len(matrix)  # broj vrhova
  k = (2 * m) / n  # prosjecni stupanj
  ro = k / (n - 1)  # gustoca mreze

  lista=[adj_list,edges,degs,k,ro]
  print(lista)

  with open(infoname, 'w') as f:
    for el in lista:
      f.write(str(el) + "\n")

size  = 12
mtx=Matrix(size)
filename='matrica.txt'
mtx.save_to_file()
matrix=open_matrix_fromTxt(filename)
print(matrix)
#show_graph( matrix)
infoname="info.txt"
save_info(matrix,infoname)
#hist_degrees(degs)


mtx.remove_edge(0,0)
mtx.save_to_file2()
matrix2=open_matrix_fromTxt("matrica2.txt")
save_info(matrix2,"info2.txt")