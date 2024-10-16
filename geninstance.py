import numpy as np
import matplotlib as plt
import networkx as nx
import sys
import math
import random
from datetime import datetime

class genInput:
    def __init__(self):
        self.N = 0 # num nodes
        self.G = nx.empty_graph()
        #self.nodes = [] #node set
        #self.edges = [] #edge set
        self.bU = 0 #upper budget
        self.bL = 0 #lower budget
        self.cU = [] #upper variable costs
        self.cL = [] #lower variable costs
        self.ciL = [] #lower variable cost increase

    #def initialise_data(self, numnodes:"int", graph_type:"str"="er", p:"float"=0.3, edge_density:"float"=0.2, *args):
    def initialise_data(self, **kwargs): 
        
        inp = {'numnodes':100, 'gtype':"er", 'p':0.3, 'edge_density':0.2, 'ws_param':0.05, 'itemcosts':50, 'budgetprop':0.6, 'test':True}
        for key,value in kwargs.items():
            inp[key] = value

        if inp['test']: 
            inp['numnodes'] = 5
            inp['gtype'] = 'wheel'
        self.N = inp['numnodes']
        if inp['gtype'] == 'er': 
            self.G = nx.erdos_renyi_graph(self.N,inp['p'])
        elif inp['gtype'] == 'ws':
            self.G = nx.empty_graph(inp['numnodes'])
            while not nx.is_connected(self.G):
                print ("GRAPH IS CONNECTED*******:", nx.is_connected(self.G))
                self.G = nx.connected_watts_strogatz_graph(inp['numnodes'],math.ceil(self.N*inp['edge_density']),
                                                           random.normalvariate(inp['p'],inp['ws_param']))
        else:
            self.G = nx.wheel_graph(inp['numnodes'])

        print ("DRAWING GRAPH")
        nx.draw(self.G)

        self.cU = np.random.randint(inp['itemcosts'],inp['itemcosts']*2,self.N)
        self.cL = np.random.randint(inp['itemcosts'],inp['itemcosts']*2,self.N)
        self.cU = list(map(lambda x: int(x), self.cU))
        self.cL = list(map(lambda x: int(x), self.cL))
        self.ciL = list(map(lambda x: x*2, self.cL))
        self.bU = math.ceil(sum(self.cU)*inp['budgetprop'])
        self.bL = math.ceil(sum(self.cL)*inp['budgetprop'])
