import cplex
import geninstance
import BiCNDP_callbacks
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






def biCNDP(**kwargs):
    # Getting the current date and time
    dt = datetime.now()
    heurFlag = 0
    # getting the timestamp
    ts = datetime.timestamp(dt)

    inp = {'numnodes':100, 'gtype':"er", 'p':0.3, 
           'edge_density':0.2, 'ws_param':0.05, 
           'itemcosts':50, 'budgetprop':0.6, 
           'dbug':True, 'test':True, 'ccut':True}
    for key,value in kwargs.items():
            inp[key] = value    
        
    random.seed(ts)
    DBUG = inp['dbug']
    I = genInput()
    I.initialise_data(numnodes=inp['numnodes'], test=inp['test'])
    n=I.N
    ccut = inp['ccut']

    cpx = cplex.Cplex()
    incmodel = cplex.Cplex()
    UB = n
    LB = 0.5
    
    
    
    #Add variables

    x = list(cpx.variables.add(obj = [0]*n,
                               lb = [0]*n,
                               ub = [1]*n,
                               types = ["B"]*n,
                               names = list(map(lambda x: "x[" +str(x)+"]", range(n)))))

    y = list(cpx.variables.add(obj = [0]*n,
                               lb = [0]*n,
                               ub = [1]*n,
                               types = ["B"]*n,
                               names = list(map(lambda x:"y[" + str(x)+"]", range(n)))))
    z = list(cpx.variables.add(obj = [0]*n,
                               lb = [0]*n,
                               ub = [1]*n,
                               types = ["B"]*n,
                               names = list(map(lambda x: "z[" + str(x)+"]", range(n)))))

    x_inc = list(incmodel.variables.add(obj = [0]*n,
                               lb = [0]*n,
                               ub = [1]*n,
                               types = ["B"]*n,
                               names = list(map(lambda x: "x[" + str(x)+"]", range(n)))))
    
    u = list(map(lambda x: [], range(n)))
    u_inc = list(map(lambda x: [], range(n)))
    usol = list(map(lambda x: 0, range(n)))
    usols = []
    uobj = []
    
    for i in range(n):
        usols += [usol]
        obj = [1.0 if j < i else 0 for j in range(n)]
        uobj.append(obj)
        u[i] = list(cpx.variables.add(obj = obj,
                               lb = [0]*n,
                               ub = [1]*n,
                               types = ["B"]*n,
                               names = list(map(lambda x: "u[" + str(x) + "][" + str(i)+"]", range(n)))))

        u_inc[i] = list(incmodel.variables.add(obj = obj,
                               lb = [0]*n,
                               ub = [1]*n,
                               types = ["B"]*n,
                               names = list(map(lambda x: "u[" + str(x) + "][" + str(i)+"]", range(n)))))

    m  = len(I.G.edges)
    r = list(cpx.variables.add(obj = [0]*m,
                               lb = [0]*m,
                               ub = [1]*m,
                               types = ["B"]*m,
                               names = list(map(lambda x: "r[" +str(x)+"]", I.G.edges))))
    
    #print ("Testing ", len([*x, *z]), len(I.cU), len([*I.cL, *I.ciL]))
    cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=y, val = list(I.cU))],
                              senses = ["L"],
                              rhs = [I.bU])

    cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[*x,*z], val = [*list(I.cL), *list(I.ciL)])],
                              senses = ["L"],
                              rhs = [I.bL])
    #add maximal sets cut
    copy_cL = [i for i in I.cL]
    for i in range(n):
        copy_cL[i] += I.bL+1 - I.cL[i]
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[*x,*z, y[i]], val = [*copy_cL, *list(I.ciL), I.ciL[i]])],
                              senses = ["G"],
                              rhs = [I.bL+1-I.cL[i]])
        copy_cL[i] -= I.bL+1-I.cL[i]
        
        

    for i in range(n):
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [z[i],y[i]], val = [1.0, -1.0])],
                                   senses = ["L"],
                                   rhs = [0.0])
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [z[i],x[i]], val = [1.0, -1.0])],
                                   senses = ["L"],
                                   rhs = [0.0])
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [z[i],x[i], y[i]], val = [-1.0, 1.0, 1.0])],
                                   senses = ["L"],
                                   rhs = [1.0])
                                   
    

    if inp['test']:cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [y[0]], val = [1])],
                                              senses = ["E"],
                                              rhs = [1.0])
    for k,(i,j) in enumerate(I.G.edges):

        [i,j]= [i,j] if i < j else [j,i]
                    
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u[i][j],x[i],x[j]], val=[1, 1, 1])],
                                   senses = ["G"],
                                    rhs = [1.0])
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[r[k],x[i]], val=[1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[r[k],x[j]], val=[1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u[i][j],x[i]], val=[1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])
        cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u[i][j],x[j]], val=[1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])
        incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u_inc[i][j],x_inc[i],x_inc[j]], val=[1, 1, 1])],
                                   senses = ["G"],
                                    rhs = [1.0])
        incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u_inc[i][j],x_inc[i]], val=[1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])
        incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u_inc[i][j],x_inc[j]], val=[1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])

        

    for i in range(n):
        for j in range(n):
            if i == j: continue
            cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u[i][j],u[j][i]], val=[1, - 1])],
                                   senses = ["E"],
                                    rhs = [0])
            incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u_inc[i][j],u_inc[j][i]], val=[1, - 1])],
                                   senses = ["E"],
                                    rhs = [0])
            for k in range(n):
                if j == k or k==i: continue
                cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u[i][j],u[j][k], u[k][i]], val=[1, 1, -1])],
                                   senses = ["L"],
                                    rhs = [1.0])
                cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u[i][j],u[j][k], u[k][i]], val=[1, -1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])
                cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u[i][j],u[j][k], u[k][i]], val=[-1, 1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])

                incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u_inc[i][j],u_inc[j][k], u_inc[k][i]], val=[1, 1, -1])],
                                   senses = ["L"],
                                    rhs = [1.0])
                incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u_inc[i][j],u_inc[j][k], u_inc[k][i]], val=[1, -1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])
                incmodel.linear_constraints.add(lin_expr = [cplex.SparsePair(ind=[u_inc[i][j],u_inc[j][k], u_inc[k][i]], val=[-1, 1, 1])],
                                   senses = ["L"],
                                    rhs = [1.0])

    xsol = list(map(lambda x: 0, range(I.N)))
    cpx.objective.set_sense(cpx.objective.sense.maximize)

    cpx.write("biCNDP.lp")
    #cpx.set_results_stream(resFile)
    #cpx.parameters.preprocessing.reduce.set(0)
    #cpx.parameters.preprocessing.presolve.set(0)
    cpx.parameters.mip.tolerances.mipgap.set(0.05)
    
    #cpx.parameters.mip.strategy.probe.set(2)
    #cpx.parameters.emphasis.mip.set(0)
    #cpx.parameters.mip.strategy.dive.set(2)
    #cpx.parameters.mip.strategy.bbinterval.set(5)
    #cpx.parameters.mip.display.set(4)
    cpx.parameters.timelimit.set(3600)
    #cpx.parameters.threads.set(1)
    #cpx.parameters.mip.strategy.backtrack.set(0.01)
    '''lazycb = cpx.register_callback(MyLazy)
    lazycb.I = I
    lazycb.dbug = DBUG
    lazycb.xvars = x
    lazycb.yvars = y
    lazycb.zvars = z
    lazycb.cutlimit = 30
    lazycb.xsol = xsol
    lazycb.usol = usols
    '''

    if ccut:
        usercb = cpx.register_callback(MyCut)
        usercb.rvars = r
        usercb.uvars = u
        usercb.xvars = x
        usercb.G = I.G.copy()
        usercb.N = I.N
        usercb.dbug = DBUG
        usercb.startnode = 0
        usercb.nodedepth = 0
        usercb.totcuts = 0
    
    
    inccb = cpx.register_callback(MyIncumbent)
    inccb.incmodel = incmodel
    inccb.dbug = DBUG
    inccb.x_inc = x_inc
    inccb.u_inc = u_inc
    inccb.ciL = I.ciL
    inccb.cL = I.cL
    inccb.bL = I.bL
    inccb.N = I.N
    inccb.xvars = x
    inccb.yvars = y
    inccb.xsol = xsol
    inccb.usol = usols
    inccb.innerobj = 0
    
    brcb = cpx.register_callback(MyBranch)
    brcb.dbug = DBUG
    brcb.yvars = y
    brcb.xvars = x
    brcb.uvars = u
    brcb.ciL = I.ciL
    brcb.cL = I.cL
    brcb.bL = I.bL
    brcb.N  = n
    brcb.uobj = uobj
    
    #lazycb.cU = I.cU
    #lazycb.cL = I.cL
    #lazycb.bU = I.bU
    #lazycb.bL = I.bL
    #lazycb.ciL = I.ciL

    #usercb = cpx.register_callback(MyCut)
    #usercb = x
    #usercb = u
    
    #cpx.register_callback(MyIncumbent)
    heurcb = cpx.register_callback(MyHeuristic)
    heurcb.xvars = x
    heurcb.dbug = DBUG
    heurcb.yvars = y
    heurcb.zvars = z
    heurcb.uvars = u
    heurcb.incmodel = incmodel
    heurcb.x_inc = x_inc
    heurcb.u_inc = u_inc
    heurcb.N = I.N
    heurcb.ciL = I.ciL
    heurcb.cL = I.cL
    
    heurcb.bL = I.bL
    heurcb.N = I.N
    heurcb.cU = I.cU
    heurcb.bU = I.bU
    #cpx.register_callback(MyCut)
    #cpx.register_callback(MyBranch)

    cpx.solve()
    solution = cpx.solution
    mip = solution.MIP
    progress = solution.progress
    xvals = solution.get_values(x)
    yvals = solution.get_values(y)
    zvals = solution.get_values(z)
    print ("Total xvars picked", xvals, yvals, zvals)
    obj = cpx.solution.get_objective_value()
    print ("objective:", obj, mip.get_mip_relative_gap())
