import pickle
from os import path
import re

class cluster():

    def __init__(self, clustname:str, rootdir:str, coldefsname=None):
        
        pathname = path.join(rootdir, clustname+'.pkl')
        with open(pathname,'rb') as f:
            self.objs = pickle.load(f)

        defname = path.join(rootdir,
            coldefsname+'.pkl' if coldefsname is not None else clustname+ '_definition.pkl')
        with open(defname,'rb') as f:
            self.coldefs = pickle.load(f)

    def findcolumns(self, rex):
        target = '.*' + rex + '.*'
        cols =   [(c, d['desc'])for c,d in self.coldefs.items() if re.match(target,d['desc'],flags=re.IGNORECASE)]
        return cols
         


if __name__ == '__main__':

    alpha_per = cluster('alpha_per','./data/J_A+A_628_A66')

    print(f'Alpha-Persei has {len(alpha_per.objs)} members')

    cols = alpha_per.findcolumns('Gaia')
    for c in cols:
        print(f'Column: {c[0]}, Description: {c[1]}')    
