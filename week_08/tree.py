from __future__ import annotations
from dataclasses import dataclass
import re
import sys

@dataclass
class Link:
    element: str | int
    next: None | Link

class Stack(object):
    def __init__(self, stack = None):
        self.stack = stack

    def __bool__(self):
        return self.stack is not None

    def __repr__(self) -> str:
        return f"Stack({repr(self.stack)}"

    def push(self, x):
        self.stack = Link(x, self.stack)

    def top(self):
        return self.stack.element

    def pop(self):
        out = self.stack.element
        self.stack = self.stack.next
        return out

    def is_empty(self):
        return self.stack is None

class Leaf(object):
    def __init__(self, name, distance = float):
        self.name = name
        self.dist = distance
    def __str__(self):
        return f'{self.name}:{self.dist}'
    __repr__=__str__
        

class Node(object):
    def __init__(self, children = (), distance = 0.0):
        self.children = children
        self.dist = distance
    def __str__(self):
        children = [str(child) for child in self.children]
        if self.dist != 0.0:
            return '({}):{}'.format(','.join(children),self.dist)
        else: return '({})'.format(','.join(children))
    __repr__=__str__
    def get_leaves(self, out = [])->list:
        for child in self.children:
            if type(child) == Leaf:
                out.append(child.name)
            elif type(child) == Node:
                child.get_leaves(out)
        return out


def Newick_parser(file):
    stack = Stack()
    lst = read_newick(file)
    for item in lst:
        if item == '(':
            stack.push(item)
        elif item[0] == ':':
            tmp = item.split(':')
            stack.top().dist = float(tmp[1].strip(', '))
            continue
        elif item == ')':
            children = []
            while (x:= stack.pop()) != '(':
                children.append(x)
            stack.push( Node(children ))
        else:
            tmp = item.split(':')
            stack.push(Leaf(tmp[0],float(tmp[1].strip(', '))))
    return stack.pop()

def read_newick(file):
    f = open(file)
    lst = [line.strip('\n') for line in f.readlines()]
    f.close()
    out = []
    for item in lst:
        if ')' in item:
            i = 0
            while i < len(item):
                if item[i] == ')':
                    out.append(item[0:i])
                    out.append(')')
                    break
                i +=1
        else:
            out.append(item)
    return out

@dataclass
class Split:
    node: Node
    leaves: node.get_leaves()

def rfdist(u:Node, v:Node):
    out = []

def get_splits(tree, out = []):
    for child in tree.children:
        if type(child) == Node:
            get_splits(child, out)
            out.append(child)
        else:
            out.append(child)
    return out


o = read_newick(sys.argv[1])
out = Newick_parser(sys.argv[1])

o = [child for child in out.children]
for child in o:
    if type(child) == Node:
        print(child.get_leaves())

out = get_splits(out)

print(out)