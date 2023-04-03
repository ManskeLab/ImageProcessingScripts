import numpy as np
from collections import defaultdict

"""
Algorithm by Keven Qiu
"""

class DisjointSet:
  def __init__(self) -> None:
    self.parent = {}
    self.rank = {}

  def inSet(self, x) -> bool:
    return x in self.parent
  
  def addSet(self, x) -> None:
    if (not self.inSet(x)):
      self.parent[x] = x
      self.rank[x] = 0
  
  def find(self, x):
    if (not self.inSet(x)):
      self.parent[x] = x
      self.rank[x] = 0

    if (x != self.parent[x]):
      self.parent[x] = self.find(self.parent[x])
    return self.parent[x]

  def union(self, x, y) -> None:
    parentX = self.find(x)
    parentY = self.find(y)

    if (x != y):
      if (self.rank[parentX] < self.rank[parentY]):
        self.parent[parentX] = parentY
      elif (self.rank[parentX] > self.rank[parentY]):
        self.parent[parentY] = parentX
      else:
        self.parent[parentY] = parentX
        self.rank[parentX] += 1

def findConnectedComponents(a: np.array, label: int) -> dict:
  size_x, size_y, size_z = a.shape
  
  def get6ConnectedNeighbours(x, y, z):
    neighbours = []
    if (x-1 >= 0):
      neighbours.append((x-1,y,z))
    if (x+1 <= size_x-1):
      neighbours.append((x+1,y,z))
    if (y-1 >= 0):
      neighbours.append((x,y-1,z))
    if (y+1 <= size_y-1):
      neighbours.append((x,y+1,z))
    if (z-1 >= 0):
      neighbours.append((x,y,z-1))
    if (z+1 <= size_z-1):
      neighbours.append((x,y,z+1))
    
    return neighbours
  
  labelToIndices = defaultdict(list)

  s = DisjointSet()

  for i in range(size_x):
    for j in range(size_y):
      for k in range(size_z):
        if (a[i,j,k] == label):
          s.addSet((i,j,k))
          neighbours = get6ConnectedNeighbours(i, j, k)
          for n in neighbours:
            if (a[n] == label):
              s.union((i,j,k), n)

  parentToKey = defaultdict(list)
  currLabel = 1
  for key in s.parent:
    if (s.find(key) not in parentToKey):
      parentToKey[s.find(key)] = currLabel
      currLabel += 1
    labelToIndices[parentToKey[s.find(key)]].append(key)

  return labelToIndices

def findMultipleLabelConnectedComponents(a: np.array, labels: list) -> dict:
  size_x, size_y, size_z = a.shape
  
  def get6ConnectedNeighbours(x, y, z):
    neighbours = []
    if (x-1 >= 0):
      neighbours.append((x-1,y,z))
    if (x+1 <= size_x-1):
      neighbours.append((x+1,y,z))
    if (y-1 >= 0):
      neighbours.append((x,y-1,z))
    if (y+1 <= size_y-1):
      neighbours.append((x,y+1,z))
    if (z-1 >= 0):
      neighbours.append((x,y,z-1))
    if (z+1 <= size_z-1):
      neighbours.append((x,y,z+1))
    
    return neighbours
  
  allLabelsToLabelToIndices = {}
  for label in labels:
    allLabelsToLabelToIndices[label] = defaultdict(list)

  s = {}
  for label in labels:
    s[label] = DisjointSet()

  for i in range(size_x):
    for j in range(size_y):
      for k in range(size_z):
        if (a[i,j,k] in labels):
          s[a[i,j,k]].addSet((i,j,k))
          neighbours = get6ConnectedNeighbours(i, j, k)
          for n in neighbours:
            if (a[n] == a[i,j,k]):
              s[a[i,j,k]].union((i,j,k), n)

  for label in labels:
    parentToKey = defaultdict(list)
    currLocalLabel = 1
    for key in s[label].parent:
      if (s[label].find(key) not in parentToKey):
        parentToKey[s[label].find(key)] = currLocalLabel
        currLocalLabel += 1
      allLabelsToLabelToIndices[label][parentToKey[s[label].find(key)]].append(key)

  return allLabelsToLabelToIndices
