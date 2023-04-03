import numpy as np
from collections import deque
from collections import defaultdict
"""
Algorithm by Keven Qiu
"""

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

  stack = deque()
  visited = np.zeros_like(a, dtype=bool)

  currLabel = 1
  labelToIndices = defaultdict(list)

  for i in range(size_x):
    for j in range(size_y):
      for k in range(size_z):
        if (not visited[i,j,k]):
          if (a[i,j,k] == label):
            labelToIndices[currLabel].append((i,j,k))
            visited[i,j,k] = True
            stack.append((i,j,k))

            while (len(stack)):
              s = stack[-1]
              stack.pop()

              if (not visited[s]):
                visited[s] = True 

              neighbours = get6ConnectedNeighbours(s[0], s[1], s[2])

              for n in neighbours:
                if (not visited[n] and a[n] == label):
                  labelToIndices[currLabel].append(n)
                  visited[n] = True
                  stack.append(n)

            currLabel += 1

  return labelToIndices
