import cplex
import numpy as np
import matplotlib as plt
import networkx as nx
import sys
import math
import random
from cplex.exceptions import CplexSolverError
import cplex.callbacks as CPX_CB
from cplex.callbacks import UserCutCallback, LazyConstraintCallback, IncumbentCallback, BranchCallback
from datetime import datetime



class MyHeuristic(CPX_CB.HeuristicCallback):
    def __call__(self):
        dbug = self.dbug
        
        node_data = self.get_node_data()
        if self.get_node_ID() % 50 != 0: return
        incum_value = self.get_incumbent_objective_value()
        if dbug: print ("Entering heuristic callback", self.get_node_ID(), self.get_current_node_depth(), incum_value)

        if node_data is None:
            
            lp_yvals = np.array(self.get_values(self.yvars))
            
            subsidy_cost = 0
            
            yvals = [0 for i in range(self.N)]
            
            for i in range(self.N):
                if np.random.random() > lp_yvals[i]: continue #do a randomised rounding on the LP value
                if subsidy_cost + self.cU[i] > self.bU: continue
                yvals[i] = 1
                subsidy_cost += self.cU[i]
            
            cost = list(map(lambda x: self.cL[x] + self.ciL[x]*yvals[x], range(self.N)))#subsidised cost of items - list
            cost_new = cost.copy()
            x_inc = self.x_inc
            u_inc = self.u_inc
            
            self.incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=x_inc, val=cost)],
                                       senses = ["L"],
                                        rhs = [self.bL], names = ["newcon"])
            for i in range(self.N):
                cost_new[i] += self.bL+1 - cost[i]
                self.incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind= x_inc, val = cost_new)],
                                                     senses = ["G"],
                                                     rhs = [self.bL+1 - cost[i]], names = ["maxcon"+str(i)])
                cost_new[i] += self.bL+1 - cost[i]
            #todo remove this constraint after solve()
            if not dbug: self.incmodel.parameters.mip.display.set(0)
            if not dbug: self.incmodel.set_results_stream(None)
            #if dbug: self.incmodel.write("innerprobHeur.lp")
            self.incmodel.solve()
            if dbug: print ("Solved LL problem", self.incmodel.solution.get_status())
            #self.incmodel.linear_constraints.delete("newCon")
            #todo do early termination to exceed incumbent  value
            
            incmodelObj = self.incmodel.solution.get_objective_value()
            
            if incmodelObj < incum_value - 1e-4: return
            
            #get u values, x values from LL problem and y values from upper level problem
            uvals = []
            solution = self.incmodel.solution
            for i in range(self.N):
                uvals.append([])
                uvals[i] = list(map(lambda x: 1 if x > 1-1e-4 else 0, solution.get_values(self.u_inc[i])))
            
            xvals = list(map(lambda x: 1 if x > 1e-4 else 0, solution.get_values(self.x_inc)))
            zvals = list(map(lambda x: 1 if xvals[x] > 1 - 1e-4 and yvals[x] > 1 - 1e-4 else 0, range(self.N)))
            index = self.xvars + self.yvars + self.zvars
            vals = xvals + yvals + zvals
            for i in range(self.N):
                index += self.uvars[i]
                vals += list(uvals[i])
            
            self.incmodel.linear_constraints.delete("newcon")
            for i in range(self.N):
                self.incmodel.linear_constraints.delete("maxcon"+str(i))
            self.set_solution([index, vals], objective_value = incmodelObj)
        else:
            
            xvals, innerobj, incval, yvals, uvals = node_data
            
            if innerobj <= incum_value - 1e-4: return
            zvals = list(map(lambda x: 1 if xvals[x] > 1 - 1e-4 and yvals[x] > 1 - 1e-4 else 0, range(self.N)))
            #zvals = xvals*yvals
            #
            index = self.xvars + self.yvars + self.zvars
            vals = xvals + yvals + zvals
            
            for i in range(self.N):
                index += self.uvars[i]
                vals += list(uvals[i])
                
            
            self.set_solution([index, vals])

        

class MyBranch(BranchCallback):
    def __call__(self):
        dbug = self.dbug
        if dbug: print ("Entering branching callback", self.get_node_ID(), self.get_current_node_depth())
        node_data = self.get_node_data()
        if node_data is not None:
            xsol, innerobj, incval, yvals, uvals = node_data
        else: 
            return
            
        if dbug: print ("Making branches in branch callback...", xsol, innerobj, incval)

        rhs = self.bL
        ind = []
        val = []
        for x in range(self.N):
            if xsol[x] < 1e-4: continue
            ind += [self.yvars[x]]
            val += [self.ciL[x]*xsol[x]]
            rhs -= self.cL[x]*xsol[x]
           
        index = []
        value = []
        
        #values = [1.0 for i in range(self.N)]
        for x in range(self.N):
            index += self.uvars[x]
            #value += values
            value += self.uobj[x]
        
        #left branch
        #TODO: add the corresponding covercuts here for the left branch
        con = [(cplex.SparsePair(ind = ind, val = val), "G", rhs+1)] #budget based cut
        
        con += [(cplex.SparsePair(ind = index, val = value), "L", incval)] #UB based on rejected incumbent's value
        #Add Nogood cut
        con += [(cplex.SparsePair(ind = self.xvars, 
                                  val = list(map(lambda x: 1.0 if x > 1 - 1e-4 else -1.0, xsol))), "L", sum(xsol) - 1) ]
        #upperbound the objective by incval
        con += [(cplex.SparsePair(ind = index, val= value), "L", incval)]

        self.make_branch(incval, constraints = con)	
        
        #rightbranch
        con = [(cplex.SparsePair(ind = ind, val = val), "L", rhs)]
        
        con += [(cplex.SparsePair(ind = index, val= value), "L", innerobj)]
        self.make_branch(innerobj, constraints = con)	

class MyIncumbent(IncumbentCallback):

    #def __init__(self, env):
    #    super().__init__(env)

    def __call__(self):

        dbug = self.dbug
        if dbug: print ("Entering incumbent callback", self.get_node_ID(), self.get_current_node_depth())
        xvals = list(map(lambda x: 1 if x > 1-1e-4 else 0, self.get_values(self.xvars)))
        yvals = list(map(lambda x: 1 if x > 1-1e-4 else 0, self.get_values(self.yvars)))
        #yvals = np.array(self.get_values(self.yvars))
        
        #zvals = list(map(lambda x: 1 if x > 1-1e-4 else 0, self.get_values(self.zvars)))
        incval =  self.get_objective_value()
        incval1 = self.get_incumbent_objective_value()
        
        cost = list(map(lambda x: self.cL[x] + self.ciL[x]*yvals[x], range(self.N)))#subsidised cost of items - list
        #cost = [self.cL[i] + self.ciL[i]*yvals[i] for i in range(self.N)] #subsidised cost of items - list

        x_inc = self.x_inc
        u_inc = self.u_inc
        
        self.incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=x_inc, val=cost)],
                                             senses = ["L"],
                                             rhs = [self.bL],
                                             names = ["newcon"])
        #todo remove this constraint after solve()
        if not dbug: self.incmodel.parameters.mip.display.set(0)
        if not dbug: self.incmodel.set_results_stream(None)
        #if dbug: self.incmodel.write("incumlp.lp")
        self.incmodel.solve()
        if dbug: print ("Solved LL problem")
        #self.incmodel.linear_constraints.delete("newCon")
        #todo do early termination to exceed incumbent  value
        incmodelObj = self.incmodel.solution.get_objective_value()
        if dbug: print ("Incumbent rejection check", incmodelObj, incval, incval1 )
        if incmodelObj < incval - 1e-4:
            if incmodelObj > incval1 + 1e-4: #only if innerobj is better than current incumbent value, set noded data
                solution = self.incmodel.solution
                #usol = sum(solution.get_values(self.u_inc))
                #usolution, yvals to be passed to heuristics. Pack them first
                uvals = []
                for i in range(self.N):
                    uvals.append([])
                    uvals[i] = list(map(lambda x: 1 if x > 1-1e-4 else 0, solution.get_values(self.u_inc[i])))
                
                xsol = list(map(lambda x: 1 if x > 1e-4 else 0, solution.get_values(self.x_inc)))
                innerobj = incmodelObj
                self.set_node_data([xsol, innerobj, incval, yvals, uvals])
            #delete the parameterised budget constraint (parametrisation on the y-values set by upper level)
            self.incmodel.linear_constraints.delete("newcon")
            if dbug: print ("Incumbent rejected at callback....")
            #TODO: use the lower level solution to set incumbent (node data that can be used in heuristiccallback)
            self.reject()
        else:
            self.incmodel.linear_constraints.delete("newcon")
            self.accept()
        
class MyCut(UserCutCallback):
    def __call__(self):
        dbug = self.dbug 
        current_node_depth = self.get_current_node_depth()
        if dbug: print ("Entering user cut callback", self.get_node_ID(), current_node_depth)
        if current_node_depth % 10 != 0: return
        G = self.G
        N = self.N
        save_total = self.totcuts
        save_nodedepth = self.nodedepth
        startnode = self.startnode
        numcuts = math.ceil(math.log(current_node_depth) if current_node_depth != 0 else 0)*30
        if save_nodedepth == current_node_depth and save_total >= 200:
            return
        if save_nodedepth != current_node_depth:
            self.nodedepth = current_node_depth
            self.totcuts = 0
            
        #if startnode != 0: print ("Startnode in usercut", startnode,numcuts, current_node_depth)
        attr = {}
        rvals = self.get_values(self.rvars)
        for k,(i,j) in enumerate(G.edges):
            attr[(i,j)] = {'capacity':rvals[k]}
        nx.set_edge_attributes(G, attr)
        totCut = 0
        nodelist = list(range(startnode, N)) + list(range(startnode))
        for t in nodelist:
            uvals = self.get_values(self.uvars[t])
            if totCut > numcuts: break
            for s in range(t+1, N):
                if dbug: print("Running min cut for", s, t)
                cut = nx.minimum_cut(G, s,t)
                if cut[0] < uvals[s] - 1e-4:
                    if dbug: print("cut violation, adding cuts:", cut[0], s,t, uvals[s])
                    nodecutset1 = []
                    nodecutset2 = []
                    nodeval1 = []
                    nodeval2 = []
                    cutset = [self.uvars[s][t]]
                    vals = [-1.0]
                    for k, (i,j) in enumerate(G.edges):
                        if i in cut[1][0] and j in cut[1][1]:
                            cutset += [self.rvars[k]]
                            vals += [1.0]
                            nodecutset1 += [self.xvars[i]]
                            nodecutset2 += [self.xvars[j]]
                            nodeval1 += [1.0]
                            nodeval2 += [1.0]
                        if j in cut[1][0] and i in cut[1][1]:
                            cutset += [self.rvars[k]]
                            vals += [1.0] 
                            nodecutset1 += [self.xvars[j]]
                            nodecutset2 += [self.xvars[i]]
                            nodeval1 += [1.0]
                            nodeval2 += [1.0]
                    if dbug: print ("cutset found", cutset)       
                    self.add(cut=cplex.SparsePair(ind=cutset, val= vals), 
                             sense= "G",
                             rhs = 0)
                    nodeval1 += [1]
                    nodeval2 += [1]
                    for p in cut[1][0]: 
                        for q in cut[1][1]:
                            nodecutset1 += [self.uvars[p][q]]
                            nodecutset2 += [self.uvars[p][q]]
                            self.add(cut=cplex.SparsePair(ind=nodecutset1, val= nodeval1), 
                                     sense= "L",
                                     rhs = len(nodeval1)-1)
                            self.add(cut=cplex.SparsePair(ind=nodecutset2, val= nodeval2), 
                                     sense= "L",
                                     rhs = len(nodeval2)-1)
                            nodecutset1.pop()
                            nodecutset2.pop()
                    totCut += 1
        
        self.startnode = t
        self.totcuts += totCut
            
        if totCut > 0 and dbug: print("Added ", totCut, "cutset cuts")
                

class MyLazy(LazyConstraintCallback):
    """Lazy constraint callback to enforce maximal subsets.
    If used then the callback is invoked for every integer feasible
    solution CPLEX finds. 
    """

    # Callback constructor. Fields 'budgets and costs' are set externally after registering the callback.
    def __init__(self, env):
        super().__init__(env)

    def __call__(self):
        dbug = self.dbug 
        if dbug: print ("Entering lazy callback", self.get_node_ID(), self.get_current_node_depth())
        dbug = self.dbug
        xvals = np.array(self.get_values(self.xvars))
        yvals = np.array(self.get_values(self.yvars))
        zvals = np.array(self.get_values(self.zvars))
        #xvals = list(map(lambda x: 1 if x > 1-1e-4 else 0, self.get_values(self.xvars)))
        #yvals = list(map(lambda x: 1 if x > 1-1e-4 else 0, self.get_values(self.yvars)))
        #zvals = list(map(lambda x: 1 if x > 1-1e-4 else 0, self.get_values(self.zvars)))
        #cost = np.array(self.I.cL) + np.array(self.I.ciL)*yvals
        cost = list(map(lambda x: self.I.cL[x] + self.I.ciL[x]*yvals[x], range(self.I.N)))#subsidised cost of items - list
        #budget = self.I.bL + 1 - np.sum(cost*xvals) #leftover budget
        budget = self.I.bL + 1 - sum(map(lambda x: cost[x]*xvals[x], range(self.I.N))) #leftover budget
    
        #Heuristically generate sets and add cuts
        cutcount = 0
        for i in range(1,self.cutlimit+1):
            if i >= 3 and cutcount/i <= 0.2: break
            if dbug: print ("Entering lazy callback3", i)
            items = list(filter(lambda x: xvals[x] < 1e-4, range(self.I.N))) #items not selected - list
            if dbug: print ("Entering lazy callback4", i, len(items),np.round(xvals,1),yvals)
            random.shuffle(items) #shuffle the list
            if dbug: print ("Entering lazy callback5", i, items, cost, budget)
            cover = []
            costsum = 0
            covercost = 0
            
            for item in items:
                if costsum + cost[item] >= budget: #if the cost of including a non-selected item exceeds the left over budget
                    if dbug: print ("checking cover")
                    if len(cover) == 0: break #cover empty
                    ind, val = [], []
                    for j in range(self.I.N): #add the cut
                        if j in cover:
                            ind += [self.yvars[j]]
                            val += [self.I.ciL[j]]
                        else:
                            ind += [self.xvars[j], self.zvars[j]]
                            val += [self.I.cL[j],self.I.ciL[j]]                           
                    self.add(constraint=cplex.SparsePair(ind = ind, val = val), sense='G', rhs=self.I.bL+1 - covercost)
                    cutcount += 1
                    #print (cover, covercost, costsum, budget, xvals, self.I.bL, yvals, 
                    #        - self.I.bL - 1 + covercost)
                    break
                else:
                    cover += [item]
                    costsum += cost[item]
                    covercost += self.I.cL[item]
                    if dbug: print ("cover,....", cover, costsum)

        #to be implemented
        #if cutcount == 0: do exact separation

        if dbug: print ("Added ", cutcount, " Lazy constraints..........") 

