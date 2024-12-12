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
        self.inp = {'numnodes':100, 'gtype':"er", 'p':0.3, 'edge_density':0.2, 'ws_param':0.05, 'itemcosts':50, 'budgetprop':0.6, 'test':True}

    #def initialise_data(self, numnodes:"int", graph_type:"str"="er", p:"float"=0.3, edge_density:"float"=0.2, *args):
    def initialise_data(self, **kwargs): 
        
        #inp = {'numnodes':100, 'gtype':"er", 'p':0.3, 'edge_density':0.2, 'ws_param':0.05, 'itemcosts':50, 'budgetprop':0.6, 'test':True}
        for key,value in kwargs.items():
            self.inp[key] = value

        inp = self.inp
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
        if inp['itemcosts'] == 1:
            self.ciL = [math.ceil(sum(self.cU)*inp['budgetprop'])]*self.N
        else:
            self.ciL = list(map(lambda x: x*2, self.cL))
        self.bU = math.ceil(sum(self.cU)*inp['budgetprop'])
        self.bL = math.ceil(sum(self.cL)*inp['budgetprop'])

        
    def write_file(self, os):
        inp = self.inp
        if os[0] == "m":
            filename = "data/{0}_{1}_{2}_{3}_{4}-{5}-{6}.txt".format(inp['gtype'], inp['numnodes'],int(inp['p']*100), int(inp['edge_density']*100), int(inp['budgetprop']*100),int(inp['itemcosts']), os[1:])
        else:
            filename = "data\{0}_{1}_{2}_{3}_{4}-{5}-{6}.txt".format(inp['gtype'], inp['numnodes'], int(inp['p']*100), int(inp['edge_density']*100), int(inp['budgetprop']*100), int(inp['itemcosts']), os[1:])

        file = open(filename, "w")
        b = " ".join(["{"+str(i)+"}" for i in range(7)])
        b = "p " + b+"\n"
        a = [i for i in list(self.inp.values())[:-1]]
        file.write(b.format(*a))
        for i,j in self.G.edges:
            file.write("e {0} {1}\n".format(i,j))
        file.write("cU "+" ".join(map(str,self.cU)))
        file.write("\ncL "+" ".join(map(str, self.cL)))
        file.write("\niL "+" ".join(map(str, self.ciL)))
        file.write("\nbU "+str(self.bU))
        file.write("\nbL "+str(self.bL))

        file.close()

    def read_file(self, **kwargs):
        
        for key,value in kwargs.items():
            self.inp[key] = value

        inp = self.inp
        os = inp['os']
        if os[0] == "m":
            filename = "data/{0}_{1}_{2}_{3}_{4}-{5}-{6}.txt".format(inp['gtype'], inp['numnodes'], int(inp['p']*100), int(inp['edge_density']*100), int(inp['budgetprop']*100),int(inp['itemcosts']), os[1:])
        else:
            filename = "data\{0}_{1}_{2}_{3}_{4}-{5}-{6}.txt".format(inp['gtype'], inp['numnodes'], int(inp['p']*100), int(inp['edge_density']*100), int(inp['budgetprop']*100), int(inp['itemcosts']), os[1:])

        file = open(filename, "r")
        
        for line in file:
            data = line.strip().split()
            if data[0] == 'p': 
                self.inp['numnodes'] = int(data[1])
                self.N = int(data[1])
                self.G = nx.empty_graph(self.inp['numnodes'])
            elif data[0] == 'e': self.G.add_edge(int(data[1]), int(data[2]))

            elif data[0] == 'cU': 
                self.cU = list(map(int, data[1:])) 
            elif data[0] == 'cL': 
                self.cL = list(map(int, data[1:]))
            elif data[0] == 'iL': 
                self.ciL = list(map(int, data[1:]))
            elif data[0] == 'bU': 
                self.bU = int(data[1])
            elif data[0] == 'bL': 
                self.bL = int(data[1])

        print (self.cU, self.cL, self.ciL, self.bU, self.bL)


