from __future__ import annotations
from dataclasses import dataclass
import re
import sys

@dataclass
class Link:
    element: str | int
    next: None | Link
    
def drop(lnk:Link, k:int)-> Link:
    if k == 0:
        return lnk
    if k > 0:
        return drop(lnk.next, k-1)

def take(lnk:Link, k:int)->Link|None:
    if k == 0:
        return None
    if k > 0:
        return Link(lnk.element, take(lnk.next, k-1))

def reverse(lnk:Link, acc = None)->Link:
    if lnk is None: 
        return acc
    acc = Link(lnk.element, acc)
    return reverse(lnk.next, acc )

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

def tokenize(tree):
    return re.findall(r'[()]|\w+', tree)
   
class Leaf(object):
    def __init__(self, name, distance = float):
        self.name = name
        self.dist = distance
    def __str__(self):
        return self.name
    __repr__=__str__
        

class Node(object):
    def __init__(self, children = (), distance = float):
        self.children = children
        self.dist = distance
    def __str__(self):
        children = [str(child) for child in self.children]
        return '({})'.format(','.join(children))
    __repr__=__str__

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

o = read_newick(sys.argv[1])
out = Newick_parser(sys.argv[1])
print(o)
print(out)
print(out.children[0].dist)