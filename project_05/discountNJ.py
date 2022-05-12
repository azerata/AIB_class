from dataclasses import dataclass
import numpy as np
import sys


class Leaf(object):
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return f'{self.name}'
    __repr__=__str__


class Node(object):
    def __init__(self, children = []):
        self.children = children
    def __str__(self):
        children = [str(child) for child in self.children]
        return '({})'.format(','.join(children))
    __repr__=__str__
    def get_leaves(self, out = [])->list:
        for child in self.children:
            if type(child) == Leaf:
                out.append(child.name)
            elif type(child) == Node:
                child.get_leaves(out)
        return out

@dataclass
class pair:
    value: float
    i: int
    j: int

def read_dist(path):
    '''
   
    '''
    file = open(path)
    tmp = int(file.readline())

    t = {}
    score = np.zeros((tmp,tmp))
    scores = [file.readline() for i in range(tmp)]

    i = 0
    for val in scores:
        tmp = [_ for _ in val.split()]
        t[tmp[0]]=i
        for j in range(len(tmp[1:])):
            score[i,j] = float(tmp[j +1])
        i+=1
    file.close()

    return t, score


def nj(t, matrix):
    s = len(t)
    itl = [l for l in t]
    T = Node()
    for i in t:
        T.children.append(Leaf(i))
    while s > 3:
        best = None
        r = [(1.0/(s-2)) * np.sum(matrix[i]) for i in range(s)]
        N = np.zeros((s,s))
        for i in range(s):
            for j in range(s):
                if i != j:
                    N[i,j] = matrix[i,j] - (r[i] + r[j])
                    if best is None:
                        best = pair(N[i,j], i, j)
                    elif N[i,j] < best.value:
                        best = pair(N[i,j], i, j)
                    
                    #Add new note to tree
        if best.i > best.j:
            T.children.append( Node(children=[T.children.pop(best.i), T.children.pop(best.j)]))
        else:
            T.children.append( Node(children=[T.children.pop(best.j), T.children.pop(best.i)]))
        #Delete row/column in D
        ND = np.pad(matrix, [(0, 1), (0, 1)], mode='constant')

        for m in range(s):
            if m == best.j or m == best.i:
                pass
            else:
                ND[m, s] = 1/2 * (matrix[best.i, m] + matrix[best.j, m] - matrix[best.i, best.j] )
                ND[s, m] = 1/2 * (matrix[best.i, m] + matrix[best.j, m] - matrix[best.i, best.j] )

        ND = np.delete(ND, [best.i, best.j], 0)
        ND = np.delete(ND, [best.i, best.j], 1)

        matrix = ND
        s -= 1
    #Final distances.

    return T




if __name__ == "__main__":
    tmp = read_dist(sys.argv[1])
    t = tmp[0]
    q = tmp[1]

    print(nj ( t, q))




