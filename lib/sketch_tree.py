# import sys
# import random
# import pandas as pd
import os
import subprocess
import math

DASHINGLOC="/home/jbonnie1/lib/dashing/dashing"


class SketchObject:
    def __init__(self, kval, sketch):
        self.k = kval
        self.sketch = sketch
        self.card = None


class SketchTreeNode:
    ''' A SketchTree node
    Public members:
        - parent: a SketchTreeNode object for a non-root node; None for the root node
        #- code: the path from root to this node
        - storage: directory(?)... file prefix for sketches
        - k_min: current minimum value of k of node sketches
        - k_max: current maximum value of k of node sketches
        - tree_row: level in binary tree
        - row_index: index in row
        - sketches: array of sketch objects ??? (would contain sketch name, sketch max d_k, sketch k)

        if this node is a leaf:
        - leaf_input: file path to input sequencing file

        if this node is a non-leaf:
        - left: a SketchTreeNode pointed leftward from this node (with edge label '0')
        - right: a SketchTreeNode pointed rightward from this node (with edge label '1')

    Private members: (x)
        - alphabets: the alphabet set of `sequence`. Each element should be unique
        - sequence: the corresponding sequence of this node
        - l_or_r: 0 or 1 depending on if left or right child
    '''

    def __init__(self, parent, row_index, tree_row, kmin=4, kmax=0, storage: str = '', leaf_input: str = '') -> None:
        self.parent = parent
        self.storage = storage
        self.leaf_input = leaf_input
        self.row_index = row_index
        self.tree_row = tree_row
        self.kmin = kmin
        self.kmax = kmax
        self.left = None
        self.right = None
        self.sketches = {}
        self.delta = 0
        self.deltak = None

        if leaf_input != '':
            if leaf_input != 'dummy':
                self._find_delta()
        else:
            pass

    #         if parent is None:
    #             # The only situation where parent doesn't exist is when this node is the first in the row
    #             self.row_index = 0
    #         else:
    #             self.row_index = (2 * parent.row_index) + l_or_r

    # if at bottom layer of nodes
    # if leaf_input != '':
    #     self.tree_row = 0
    # else:
    #     self.tree_row = self.left.tree_row + 1
    def is_leaf(self) -> bool:
        ''' Return True if this node is a leaf. '''
        return self.leaf_input != '' and self.leaf_input != 'dummy'

    def is_dummy(self) -> bool:
        ''' Return True if this node is a leaf. '''
        return self.leaf_input == 'dummy'

    def _get_leaf_sketch(self, kval):
        sketchdir = os.path.join(self.storage, "k" + str(kval))
        os.makedirs(sketchdir, exist_ok=True)
        sketchloc = os.path.join(sketchdir, str(abs(hash(self.leaf_input))))
        sketch = SketchObject(kval=kval, sketch=sketchloc)
        sketch_call = subprocess.run(
            ["/home/jbonnie1/lib/dashing/dashing", "sketch", "-k" + str(kval), "-p10", "-o", str(sketchloc),
             self.leaf_input])
        print(sketch_call)
        # cmd = ["/home/jbonnie1/lib/dashing/dashing", "card", "--presketched", sketchloc]
        # card_proc = subprocess.run(cmd, capture_output=True, text=True)
        # card_lines=card_proc.stdout.strip().split('\n')
        # ##TODO: this could be how we create a bunch of sketches from the same output table!!!
        # for line in card_lines[1:]:
        #     y = line.strip().split("\t")
        #     sketch.card=float(y[1])
        # print("sketch: ", sketch)
        return sketch

    def _get_union_sketch(self, kval):

        kstr = "k" + str(kval)
        sketchdir = os.path.join(self.storage, kstr)
        os.makedirs(sketchdir, exist_ok=True)
        left_sketch = self.left.sketches[kstr].sketch
        right_sketch = self.right.sketches[kstr].sketch
        ## RIGHT NOW ASSUMING LEFT AND RIGHT EXIST
        # sketchdir = os.path.join(self.storage,"k"+ str(kval))
        # os.makedirs(sketchdir, exist_ok=True)
        sketchloc = os.path.join(sketchdir, str(abs(hash(left_sketch + right_sketch))))
        sketch = SketchObject(kval=kval, sketch=sketchloc)
        sketch_call = subprocess.run(
            ["/home/jbonnie1/lib/dashing/dashing", "union", "-k" + str(kval), "-p10", "-o", str(sketchloc), left_sketch,
             right_sketch])
        print(sketch_call)
        return sketch

    def fill_sketches(self, kmin, kmax):
        '''
        when you already have the sketches you need for delta, but your tree needs more in case user asks
        '''
        while self.kmin > kmin:
            self.get_sketch(self.kmin - 1)
            self.kmin -= 1
        while self.kmax < kmax:
            self.get_sketch(self.kmin + 1)
            self.kmin += 1

    def _get_sketch(self, kval):
        sketchdir = os.path.join(self.storage, "k" + str(kval))
        os.makedirs(sketchdir, exist_ok=True)

        if self.is_leaf:
            new_sketch = self._get_leaf_sketch(kval)
        if self.is_dummy:
            pass
        else:
            new_sketch = self._get_union_sketch(kval)

        cmd = ["/home/jbonnie1/lib/dashing/dashing", "card", "--presketched", new_sketch.sketch]
        card_proc = subprocess.run(cmd, capture_output=True, text=True)
        card_lines = card_proc.stdout.strip().split('\n')
        ##TODO: this could be how we create a bunch of sketches and cards from the same output table!!!
        for line in card_lines[1:]:
            y = line.strip().split("\t")
            new_sketch.card = float(y[1])
        print("sketch: ", new_sketch)
        return new_sketch

    def _find_delta(self):
        last_delta = 0
        k = self.kmin - 1

        # go right from min k
        while self.delta <= last_delta:
            self.deltak = k
            self.delta = last_delta
            k += 1
            #           #TODO: can make this return a list of sketches my dudes!!!
            new_sketch = self._get_sketch(k)
            # self.sketches.append(new_sketch)
            self.sketches['k' + str(k)] = new_sketch
            last_delta = new_sketch.card / new_sketch.k
        while k < self.kmax:
            k += 1
            new_sketch = self._get_sketch(k)
            # self.sketches.append(new_sketch)
            self.sketches['k' + str(k)] = new_sketch
        self.kmax = k


class SketchTree:
    '''
    Build a sketch tree.
    input:
        - files: list of sequence file locations
        -
    properties:
        - real_leaves: list of real (non-dummy) terminal nodes
        - kmin
        - kmax
        - nodes(?)

    '''

    def __init__(self, sequence_files: list, kmin=7, kmax=18):
        self.kmin = kmin
        self.kmax = kmax
        self.leaves = []
        self.n = len(sequence_files)
        self.depth = math.ceil(math.log2(self.n))
        self.nleaves = math.pow(2, self.depth)
        self.st = [[None] * (2 ** n) for n in range(self.nleaves, -1, -1)]
        pass

        return

    def _build_tree():
        pass

    def _traverse_tree(self, node, depth, code):
        ''' Traverse a SketchTree recursively. '''
        if node.left:
            self._traverse_tree(node.left, depth + 1, code + '0')
        if node.right:
            self._traverse_tree(node.right, depth + 1, code + '1')
        if not (node.left or node.right):
            pass
            # self._code_lengths.append(depth)
            # self._code_words.append(code)
            # self._symbols.append(node.symbol)

