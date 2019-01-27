# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 19:37:05 2019

@author: Quentin
"""

import cplex
from cplex.callbacks import LazyConstraintCallback
import time
import numpy as np
import random

import files

# Fonctions utilitaires

def getVarName(var, m):
    return var + '_' + str(m[0]) + '_' + str(m[1])

def getNodeIndex(node):
    return unique.index(node)
        
def printSolution(objectif, solution):
    print "Objectif :", objectif
    arks = []
    for i in range(len(Mat)):
        if solution[i] != 0:
            arks.append((Mat[i][0], Mat[i][1]))
    chemin = "Chemin : "
    curr = s
    while curr != t:
        for ark in arks:
            if ark[0] == curr:
                chemin+= str(curr) + " - "
                curr = ark[1]
                break
    chemin += str(t)
    print chemin
    return chemin


# Lecture des données & pre-processing
def readDataFromFile(filename):
    with open(filename, 'r') as f:
        global n, s, t, S, d1, d2, p, ph, Mat, Dist, Delta, unique
        Mat = []
        n = int(f.readline().split("=")[1])
        Dist = [[0]*n for _ in range(n)]
        Delta = [[0.0]*n for _ in range(n)]
        s = int(f.readline().split("=")[1])
        t = int(f.readline().split("=")[1])
        S = int(f.readline().split("=")[1])
        d1 = float(f.readline().split("=")[1])
        d2 = float(f.readline().split("=")[1])
        #print n,s,t,S,d1,d2
        p = [int(x) for x in f.readline().split("[")[1].split(']')[0].split(',')]
        ph = [int(x) for x in f.readline().split("[")[1].split(']')[0].split(',')]
        f.readline()
        for line in f:
            x = line.split(";")[0].split(']')[0].split()
            Mat.append([int(x[0]),int(x[1]), int(x[2]), float(x[3])])
            Dist[int(x[0])-1][int(x[1])-1] = int(x[2])
            Delta[int(x[0])-1][int(x[1])-1] = float(x[3])
    # Pre-processing : retrait des noeuds doublons
    #nei = [[i+1 for i in range(n) if Dist[i] == Dist[j] and Delta[i] == Delta[j] and p[i] == p[j]] for j in range(n)]
    #dbl = np.unique(nei)
    dbl = []
    already = []
    for i in range(n):
        if i not in already:
            db = [i+1]
            for j in range(i+1,n):
                if Dist[i] == Dist[j] and Delta[i] == Delta[j] and p[i] == p[j]:
                    db.append(j+1)
                    already.append(j)
            dbl.append(db)
            
    if len(dbl) < n:
        unique = []
    
        for db in dbl:
            if s in db:
                unique.append(s)
            elif t in db:
                unique.append(t)
            else:
                ph_dbl = [ph[node-1] for node in db]
                if len(ph_dbl) > 0:
                    ind = np.argmin(ph_dbl)
                    unique.append(db[ind])
        dbl = np.unique([[i for i in unique if Dist[i-1] == Dist[j-1] and Delta[i-1] == Delta[j-1]] for j in unique])
        if len(dbl) < len(unique):
            unique = []
            for db in dbl:
                if s in db:
                    unique.append(s)
                elif t in db:
                    unique.append(t)
                else:
                    ind = np.argmin([p[node-1] for node in db])
                    if ph[db[ind]-1] == np.min([ph[node-1] for node in db]):
                        unique.append(db[ind])
                    else:
                        unique += db
    Dist = [[Dist[j-1][i-1] for i in unique] for j in unique]
    Delta = [[Delta[j-1][i-1] for i in unique] for j in unique]
    Mat = [m for m in Mat if m[0] in unique and m[1] in unique]
    p = [p[i-1] for i in unique]
    ph = [ph[i-1] for i in unique]
    n = len(unique)
    
class CuttingPlaneCallback(LazyConstraintCallback):
    def __call__(self):
        # Enregistrement des résultats du problème maître
        m_star = []
        sum_x_i = [0] * n
        x_star = []
        z = round(self.get_values('z'),5)
        
        for m in Mat:
            x_star.append(self.get_values(getVarName('x',m)))
            if self.get_values(getVarName('x',m)) > 0:
                m_star.append(m)
                sum_x_i[getNodeIndex(m[0])] +=1
        
        # Sous-problème lié à l'objectif
        obj_val, delta_1_star = worstObjective_heu(x_star)
        if obj_val > z:
            ind = ['z'] + [getVarName('x',m) for m in Mat]
            val = [1] + [-Mat[i][2] * (1+delta_1_star[i]) for i in range(len(Mat))]
            self.add(constraint=cplex.SparsePair(ind=[name2idx[name] for name in ind],\
                                                 val=val),\
                     rhs=0, sense="G")
        else:
            con_val, delta_2_star = worstWeight_heu(x_star)
            if con_val > S:
                p_2 = [p[getNodeIndex(i)] + delta_2_star[getNodeIndex(i)] * ph[getNodeIndex(i)] for i in unique]
                ind = [getVarName('x',m) for m in Mat]
                val = [p_2[getNodeIndex(m[0])] for m in Mat]
                self.add(constraint=cplex.SparsePair(ind=[name2idx[name] for name in ind], val=val),\
                         rhs=S - p_2[getNodeIndex(t)], sense="L")
            




def createVariables(cpl, z = False, dual = False):
    global name2idx
    # Ajout des variables X
    if z:
        cpl.variables.add(names=[getVarName('x',m) for m in Mat],\
                  obj=[0] * len(Mat),\
                  types=[cpl.variables.type.binary] * len(Mat))
        
        cpl.variables.add(names=["z"], obj = [1], types=[cpl.variables.type.continuous])
    else:
        cpl.variables.add(names=[getVarName('x',m) for m in Mat],\
                  obj=[m[2] for m in Mat],\
                  types=[cpl.variables.type.binary] * len(Mat))
    # Ajout des variables duales
    if dual:
        cpl.variables.add(names=["alpha"], obj = [d1], types=[cpl.variables.type.continuous])
    
        cpl.variables.add(names=[getVarName('beta',m) for m in Mat],\
                          obj=[m[3] for m in Mat], \
                          types=[cpl.variables.type.continuous] * len(Mat))
        
        cpl.variables.add(names=["gamma"], obj = [0], types=[cpl.variables.type.continuous])
        
        cpl.variables.add(names=["epsilon_"+str(i+1) for i in range(n)], obj = [0]*n,\
                             types=[cpl.variables.type.continuous] * n)
        
    name2idx = {n : j for j,n in enumerate(cpl.variables.get_names())}




def flowConstraints(cpl):
    # Contraintes de flot
    c = []
    values = []
    rhs = []
    senses = ""
    
    for i in range(n):
        c.append([])
        values.append([])
        rhs.append(1 if i == getNodeIndex(s) else -1 if i == getNodeIndex(t) else 0)
        senses += "E"
    
    for m in Mat:
        ## Contraintes de flot
        # Arcs entrants
        c[getNodeIndex(m[0])].append(getVarName('x',m))
        values[getNodeIndex(m[0])].append(1)
        # Arcs sortants
        c[getNodeIndex(m[1])].append(getVarName('x',m))
        values[getNodeIndex(m[1])].append(-1)

    cpl.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[name2idx[name] for name in c[i]], val=values[i]) for i in range(n)],\
                           rhs=rhs, senses=senses)



def dualConstraints(cpl):
    # Contraintes duales de l'objectif
    c_d = []
    
    # Contraintes duales de poids
    c_w = []
    val_w = []
    rhs_w = []
    senses_w = []
    
    for i in range(n):
        c_w.append(["gamma", "epsilon_"+str(i+1)])
        val_w.append([1,1])
        rhs_w.append(0 if i!=getNodeIndex(t) else ph[getNodeIndex(t)])
        senses_w += "G"
    c_w.append(["gamma"] + ["epsilon_"+str(i) for i in range(1,n+1)])
    val_w.append([d2] + [2] * n)
    rhs_w.append(S - p[getNodeIndex(t)])
    senses_w += "L"

    for m in Mat:
        ## Contraintes duales objectifs
        c_d.append(["alpha", getVarName("beta",m),getVarName("x",m)])
        
        ## Contraintes duales poids
        if(m[0] != t):
            c_w[getNodeIndex(m[0])].append(getVarName("x",m))
            val_w[getNodeIndex(m[0])].append(-ph[getNodeIndex(m[0])])
        
        c_w[n].append(getVarName("x",m))
        val_w[n].append(p[getNodeIndex(m[0])])
    
    # Contraintes duales de l'objectif    
    cpl.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[name2idx[name] for name in c_d[i]], val=[1,1,-Mat[i][2]]) for i in range(len(Mat))],\
                                rhs=[0 for _ in c_d], senses = ['G' for _ in c_d])
    # Contraintes duales de poids
    cpl.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[name2idx[name] for name in c_w[i]], val=val_w[i]) for i in range(n+1)],\
                            rhs=rhs_w, senses=senses_w)
    
    
def addObjectiveCut(cpl, delta_1, matr):
    ind = ['z'] + [getVarName('x',m) for m in matr]
    val = [1] + [-matr[i][2] * (1+delta_1[i]) for i in range(len(matr))]
    cpl.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[name2idx[name] for name in ind], val=val)],\
                                rhs=[0], senses="G")
    
    
def addWeightCut(cpl, delta_2):
    def getP(d):
        return p[getNodeIndex(d)] + delta_2[getNodeIndex(d)] * ph[getNodeIndex(d)]
    ind = [getVarName('x',m) for m in Mat]
    val = [getP(m[0])for m in Mat]
    cpl.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[name2idx[name] for name in ind], val=val)],\
                                rhs=[S - getP(t)], senses="L")

def worstObjective(matr):
    cpl_obj = cplex.Cplex()
    
    cpl_obj.set_log_stream(None)
    cpl_obj.set_results_stream(None)

    
    cpl_obj.objective.set_sense(cpl_obj.objective.sense.maximize)
    
    ub = [m[3] for m in matr]
    
    cpl_obj.variables.add(names=[getVarName('delta_1',m) for m in matr],\
                          obj=[m[2] for m in matr],\
                          ub = ub)
    
    cpl_obj.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[i for i in range(len(matr))], val=[1] * len(matr))],
                                    rhs = [d1], senses="L")
    cpl_obj.solve()
    
    obj_coef = []
    obj_val = 0
    
    for m in matr:
        coef = cpl_obj.solution.get_values(getVarName('delta_1',m))
        obj_coef.append(coef)
        obj_val += m[2] * (1+coef)

    return round(obj_val,5), obj_coef

def worstObjective2(x_star):
    m_star = []
    for i in range(len(Mat)):
        if x_star[i] > 0:
            m_star.append(Mat[i])
            
    obj_val, delta_1_star = worstObjective(m_star)
    delta_1_star_all = [0] * len(Mat)
    for i in range(len(delta_1_star)):
        delta_1_star_all[Mat.index(m_star[i])] = delta_1_star[i]
    return obj_val, delta_1_star_all
    
  
def worstObjective_heu(x_star):
    m_star = []
    for i in range(len(Mat)):
        if x_star[i] > 0:
            m_star.append(Mat[i])
    
    m_star.sort(key=lambda x: x[2], reverse=True)
    obj_val = sum([m[2] for m in m_star])
    delta_1_star = [0] * len(Mat)
    d = d1
    i = 0
    while d > 0 and i < len(m_star):
        delta = min([d,m_star[i][3]])
        d -= delta
        delta_1_star[Mat.index(m_star[i])] = delta
        obj_val += m_star[i][2] * delta
        i += 1
    
    return round(obj_val,5), delta_1_star


          

def worstWeight(y): # where y_i = sum x_ij
    cpl_con = cplex.Cplex()
    cpl_con.set_log_stream(None)
    cpl_con.set_results_stream(None)

    
    cpl_con.objective.set_sense(cpl_con.objective.sense.maximize)
    
    obj_d2 = [ph[i]*y[i] for i in range(n)]
    obj_d2[getNodeIndex(t)] += ph[getNodeIndex(t)] 
    
    cpl_con.variables.add(names = ["delta_2_"+str(i) for i in range(n)],\
                                   obj=obj_d2,\
                                   ub = [2] * n)
    cpl_con.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[i for i in range(n)], val = [1]*n)],\
                                             rhs = [d2], senses="L")
    cpl_con.solve()

    con_coef = []
    con_val = 0
    for i in range(n):
        coef = cpl_con.solution.get_values('delta_2_' + str(i))
        con_coef.append(coef)
        con_val += (p[i] + coef * ph[i]) * (y[i] if i != getNodeIndex(t) else 1)

    return con_val, con_coef
    

def worstWeight2(x_star):
    y = [0] * n
    for i in range(len(Mat)):
        if x_star[i] > 0:
            y[getNodeIndex(Mat[i][0])] += 1
    con_val, delta_2_star = worstWeight(y)
    return con_val, delta_2_star

def worstWeight_heu(x_star):
    delta_2_star = [0] * n
    y = [t]
    for i in range(len(Mat)):
        if x_star[i] > 0:
            y.append(Mat[i][0])
    y.sort(key=lambda x: ph[getNodeIndex(x)], reverse = True)
    con_val = sum(p[getNodeIndex(x)] for x in y)
    d = d2
    i = 0
    while d > 0 and i < len(y):
        delta = min([d, 2])
        d -= delta
        delta_2_star[getNodeIndex(y[i])] = delta
        con_val += delta * ph[getNodeIndex(y[i])]
        i += 1
    return con_val, delta_2_star


def cuttingPlanes():
    cpl = cplex.Cplex()    
    cpl.set_log_stream(None)
    cpl.set_results_stream(None)

    createVariables(cpl, z=True)
    
    flowConstraints(cpl)
    
    
    sum_dist = sum([sum(d) for d in Dist])
    addObjectiveCut(cpl, [min(d1*m[2]/sum_dist, m[3]) for m in Mat], Mat)
    
    sum_ph = sum(ph)
    addWeightCut(cpl, [min(ph[i]*d2/sum_ph, 2) for i in range(n)])
    
    is_opt = False
    while not is_opt:
        cpl.solve()
        
        # Enregistrement des résultats du problème maître
#        m_star = []
#        sum_x_i = [0] * n
        x_star = []
        z = round(cpl.solution.get_values('z'),5)
        
        for m in Mat:
            x_star.append(cpl.solution.get_values(getVarName('x',m)))
#            if cpl.solution.get_values(getVarName('x',m)) > 0:
#                m_star.append(m)
#                sum_x_i[getNodeIndex(m[0])] +=1

        # Sous-problème lié à l'objectif
        obj_val, delta_1_star = worstObjective_heu(x_star)
        if obj_val > z:
            addObjectiveCut(cpl, delta_1_star, Mat)
        else:
            is_opt = True
        
            # Sous-problème lié au poids
            con_val, delta_2_star = worstWeight_heu(x_star)
            if con_val > S:
                addWeightCut(cpl, delta_2_star)
                is_opt = False

    return cpl.solution.get_objective_value(), cpl.solution.get_values()


def branchAndCut():
    cpl = cplex.Cplex()
    
    cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    
    createVariables(cpl, z=True)
    
    flowConstraints(cpl)
    
    
    sum_dist = sum([sum(d) for d in Dist])
    addObjectiveCut(cpl, [min(d1*m[2]/sum_dist, m[3]) for m in Mat], Mat)
    
    sum_ph = sum(ph)
    addWeightCut(cpl, [min(ph[i]*d2/sum_ph, 2) for i in range(n)])
    
    cpl.register_callback(CuttingPlaneCallback)
    cpl.solve()

    return cpl.solution.get_objective_value(), cpl.solution.get_values()

def dualSolve():
    cpl = cplex.Cplex()
    
    cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    
    createVariables(cpl, dual=True)
    
    flowConstraints(cpl)
    
    dualConstraints(cpl)
    
    cpl.solve()
    
    return cpl.solution.get_objective_value(), cpl.solution.get_values()


def nonRobustSolve():
    cpl = cplex.Cplex()
    
    cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    
    createVariables(cpl, z=True)

    flowConstraints(cpl)
    
    addObjectiveCut(cpl, [0]*len(Mat), Mat)
    addWeightCut(cpl, [0]*n)

    cpl.solve()
    
    return cpl.solution.get_objective_value(), cpl.solution.get_values()


def heuristicSolve_____():
    ArkMap = [[-1]*n for _ in range(n)]
    for i in range(len(Mat)):
        ArkMap[getNodeIndex(Mat[i][0])][getNodeIndex(Mat[i][1])] = i

    P = []
    Pb = [x for x in unique]
    d = [float("inf")] * n
    pred = [-1] * n
    d[getNodeIndex(s)] = 0


    while len(P) < n:
        a = sorted([x for x in Pb],key=lambda x: d[getNodeIndex(x)])[0]
        P.append(a)
        Pb.remove(a)
        
        if d[getNodeIndex(a)] >= d[getNodeIndex(t)]:
            print "break", len(P), n
            break
        for b in Pb:
                if Dist[getNodeIndex(a)][getNodeIndex(b)] > 0:
                    #distance = d[getNodeIndex(a)] + Dist[getNodeIndex(a)][getNodeIndex(b)] * (1+Delta[getNodeIndex(a)][getNodeIndex(b)])
                    distance = d[getNodeIndex(a)] + p[getNodeIndex(b)] # + 2*ph[getNodeIndex(b)]
                    #wei = weight[getNodeIndex(a)] + p[getNodeIndex(b)]
                    if d[getNodeIndex(b)] > distance:# and wei <= S - wei_add:
                        pred[getNodeIndex(b)] = a
                        d[getNodeIndex(b)] = distance
                        #weight[getNodeIndex(b)] = wei


    curr = t
    path = []
    values = [0] * len(Mat)
    #print pred, d
    while curr != s:
        path = [curr] + path
        values[ArkMap[getNodeIndex(pred[getNodeIndex(curr)])][getNodeIndex(curr)]] = 1
        curr = pred[getNodeIndex(curr)]
    
    path = [s] + path
    inside = [x for x in path if x != s and x != t]
    
    con_val, delta_2_star = worstWeight_heu(values)
#    while con_val > S:
#        change = False
#        for remove in sorted(inside, key=lambda x: p[getNodeIndex(x)] + delta_2_star[getNodeIndex(x)] * ph[getNodeIndex(x)]):
#            before = path[path.index(remove)-1]
#            after = path[path.index(remove)+1]
#            for add in sorted([node for node in unique if node not in path and\
#                                                            ArkMap[getNodeIndex(before)][getNodeIndex(node)] >= 0 and\
#                                                            ArkMap[getNodeIndex(node)][getNodeIndex(after)] >= 0],
#                                key=lambda x: p[getNodeIndex(x)] + delta_2_star[getNodeIndex(x)] * ph[getNodeIndex(x)]):
##                if p[getNodeIndex(add)] + delta_2_star[getNodeIndex(remove)] * ph[getNodeIndex(add)] > p[getNodeIndex(remove)] + delta_2_star[getNodeIndex(remove)] * ph[getNodeIndex(remove)]:
##                    continue
#                print path, remove, add
#                values[ArkMap[getNodeIndex(before)][getNodeIndex(add)]] = 1
#                values[ArkMap[getNodeIndex(add)][getNodeIndex(after)]] = 1
#                values[ArkMap[getNodeIndex(before)][getNodeIndex(remove)]] = 0
#                values[ArkMap[getNodeIndex(remove)][getNodeIndex(after)]] = 0
#                
#                path[path.index(remove)] = add
#
#                inside[inside.index(remove)] = add
#                
#                change = True
#                break
#            if change:
#                break
#        if change:    
#            print change
#            con_val, delta_2_star = worstWeight_heu(values)
#            print con_val, S
#            continue
#        else:
#            print "BROKEN"
#            break
    if con_val > S:
        print "BROKEN"
        obj_val = 0
    else:
        obj_val, _ = worstObjective_heu(values)
    
    
    return obj_val, values

def heuristicSolve():
    ArkMap = [[-1]*n for _ in range(n)]
    for i in range(len(Mat)):
        ArkMap[getNodeIndex(Mat[i][0])][getNodeIndex(Mat[i][1])] = i

    banned = []
    saved = []
    saving = -1
    feasible = False
    
    while True:
        print "BAN", banned
        print "SAVED", saved
        P = []
        Pb = [x for x in unique if x not in banned]
        d = [float("inf")] * n
        pred = [-1] * n
        d[getNodeIndex(s)] = 0
    
    
        while len(P) < n:
            a = sorted([x for x in Pb],key=lambda x: d[getNodeIndex(x)])[0]
            P.append(a)
            Pb.remove(a)
            
            if d[getNodeIndex(a)] >= d[getNodeIndex(t)]:
                break
            for b in Pb:
                    if Dist[getNodeIndex(a)][getNodeIndex(b)] > 0:
                        distance = d[getNodeIndex(a)] + Dist[getNodeIndex(a)][getNodeIndex(b)] * (1+Delta[getNodeIndex(a)][getNodeIndex(b)])
                        #wei = weight[getNodeIndex(a)] + p[getNodeIndex(b)]
                        if d[getNodeIndex(b)] > distance:# and wei <= S - wei_add:
                            pred[getNodeIndex(b)] = a
                            d[getNodeIndex(b)] = distance
                            #weight[getNodeIndex(b)] = wei
        
        if d[getNodeIndex(t)] == float("inf"):
            saving = banned[0]
            saved.append(saving)
            banned = []
        else:
        
            curr = t
            path = []
            values = [0] * len(Mat)
            #print pred, d
            while curr != s:
                path = [curr] + path
                values[ArkMap[getNodeIndex(pred[getNodeIndex(curr)])][getNodeIndex(curr)]] = 1
                curr = pred[getNodeIndex(curr)]
            
            path = [s] + path
            inside = [x for x in path if x != s and x != t and x != saving]
            
            
            con_val, _ = worstWeight_heu(values)
            
            if con_val > S:
                if inside == []:
                    saving = banned[0]
                    saved.append(saving)
                    banned = []
                    
                else:
                    print path, inside
                    banned.append(sorted(inside, key=lambda x: ph[getNodeIndex(x)])[0])
            else:
                feasible = True
                break
        if feasible:
            break
    
    obj_val, _ = worstObjective_heu(values)
    
    
    return obj_val, values


#modes = ["CP", "dual", "b&c"]
modes = ["b&c"]

results = {mode : [] for mode in modes}


for f in files.files[:]:
    for mode in modes:
        start_time = time.clock()
        print f
        
        readDataFromFile(f)
        
        data_time = time.clock()
        if mode == "CP":
            objectif, solution = cuttingPlanes()
        if mode == "dual":
            objectif, solution = dualSolve()
        if mode == "b&c":
            objectif, solution = branchAndCut()
        if mode == "classic":
            objectif, solution = nonRobustSolve()
        if mode == "heuristic":
            objectif, solution = heuristicSolve_____()
        
        solve_time = time.clock() - data_time
        exec_time = time.clock() - start_time
        
        path = printSolution(objectif, solution)
        results[mode].append((f, objectif, path, exec_time, solve_time))
        
        print "time :", exec_time, "s"
        print "dont data :", data_time - start_time
        print "dont solving :", solve_time
        print "\n\n"
        break
        if objectif == 0:
            print "STOP"
            break
for i in range(len(results[modes[0]])):
    print results[modes[0]][i][0]
    for mode in modes:
        print results[mode][i][1], "(", results[mode][i][3], "s., solving", results[mode][i][4], ",s.)"








































