#============================================================================#
#                             Data Tree Module                               #
#============================================================================#

class tree_entry(object):
    def __init__(self, in_dict):
        self._dict = in_dict
        for key in self._dict:
            setattr(self, key, self._dict[key])

    def __str__(self):
        return self._dict.__str__()

    def __repr__(self):
        return self.__str__()

class data_tree(object):
     def __init__(self, name=None, branches=[]):
         '''
         Data Tree class for storing properties of a data event

         Arguments
         ---------

         name      (string)       Name of Data Tree

        
         Optional Arguments
         ------------------

         branches  (list of strings)   Names of Branches to create
         '''
         self._tree = {}
         self.name = name
         if branches:
             for b in branches:
                 self._tree[b] = []

     def add_branch(self, name, content=[]):
         '''
         Add a new branch to the data tree

         Arguments
         ---------

         name   (string)    Name of new branch to create

        
         Optional Arguments
         ------------------
 
         content (list of any type)      Data to store in branch
         '''
         self._tree[name] = [i for i in content] if not isinstance(content, list) else content

     def fill(self, values):
         '''
         Fill values into branches
 
         Arguments
         ---------

         values  (list of values)       List of values for each branch
         '''

         assert len(values) == len(list(self._tree.keys())), "Number of Values must match Number of Branches"

         for v, branch in zip(values, self._tree):
             self._tree[branch].append(v)

     def getN(self):
         '''Get Number of Entries in Data Tree'''
         try:
             return len(list(self._tree.values())[0])
         except IndexError:
             return None

     def getEntry(self, i):
         '''Return the values for the i-th entry in the Data Tree'''
         _output = {}
         for branch in self._tree:
             _output[branch] = self._tree[branch][i]
         return tree_entry(_output)

     def __str__(self):
         _str='''
====================================
  DataTree : {}
  Entries  : {}
====================================
         '''.format(self.name, self.getN())

         for branch in self._tree:
             _str += ' {}       {}\n'.format(branch, len(self._tree[branch]))

         return _str

     def draw(self, branch, bins=100, *args, **kwargs):
         import matplotlib.pyplot as plt
         plt.hist(self._tree[branch], bins=bins, *args, **kwargs)
         plt.title('A Histogram of {}'.format(branch))
         plt.xlabel(branch)
         plt.show()
