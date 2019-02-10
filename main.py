# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 19:37:05 2019

@author: Quentin PENA, Guillaume FAVIER
"""

import cplex
from cplex.callbacks import LazyConstraintCallback
import time
import numpy as np

import files

# Fonctions utilitaires

# Retourne le nom de la variable indexée par l'arc m
def getVarName(var, m):
    return var + '_' + str(m[0]) + '_' + str(m[1])

# Retourne l'index du noeud dans le graphe réduit
def getNodeIndex(node):
    return unique.index(node)

# Affiche le chemin solution et sa valeur
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
    ## Pre-processing : retrait des noeuds doublons
    # Construction des listes de noeuds équivalents
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
    
    # Choix du meilleur élément équivalent pour chaque liste de noeuds équivalents        
    if len(dbl) < n:
        unique = []
    
        for db in dbl:
            # On garde la source ou la destination en priorité
            if s in db:
                unique.append(s)
            elif t in db:
                unique.append(t)
            # Sinon, on prend le noeud de plus petit p chapeau
            else:
                ph_dbl = [ph[node-1] for node in db]
                if len(ph_dbl) > 0:
                    ind = np.argmin(ph_dbl)
                    unique.append(db[ind])
        # On construit les listes de noeuds jumeaux pour les noeuds que l'on a gardés
        dbl = np.unique([[i for i in unique if Dist[i-1] == Dist[j-1] and Delta[i-1] == Delta[j-1]] for j in unique])
        # On choisit le meilleur élément pour chaque liste
        if len(dbl) < len(unique):
            unique = []
            for db in dbl:
                # On garde la source et la destination en priorité
                if s in db:
                    unique.append(s)
                elif t in db:
                    unique.append(t)
                # Si l'élément de plus petit p et aussi celui de plus petit p chapeau
                # on ne garde que lui, sinon on les garde tous
                else:
                    ind = np.argmin([p[node-1] for node in db])
                    if ph[db[ind]-1] == np.min([ph[node-1] for node in db]):
                        unique.append(db[ind])
                    else:
                        unique += db
    # On met à jour le graphe en ne gardant que les noeuds choisis
    Dist = [[Dist[j-1][i-1] for i in unique] for j in unique]
    Delta = [[Delta[j-1][i-1] for i in unique] for j in unique]
    Mat = [m for m in Mat if m[0] in unique and m[1] in unique]
    p = [p[i-1] for i in unique]
    ph = [ph[i-1] for i in unique]
    n = len(unique)


# Callback pour le branch & cut
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
            # On ajoute la coupe
            ind = ['z'] + [getVarName('x',m) for m in Mat]
            val = [1] + [-Mat[i][2] * (1+delta_1_star[i]) for i in range(len(Mat))]
            self.add(constraint=cplex.SparsePair(ind=[name2idx[name] for name in ind],\
                                                 val=val),\
                     rhs=0, sense="G")
        else:
            # Sous-problème lié au poids
            con_val, delta_2_star = worstWeight_heu(x_star)
            if con_val > S:
                # On ajoute la coupe
                p_2 = [p[getNodeIndex(i)] + delta_2_star[getNodeIndex(i)] * ph[getNodeIndex(i)] for i in unique]
                ind = [getVarName('x',m) for m in Mat]
                val = [p_2[getNodeIndex(m[0])] for m in Mat]
                self.add(constraint=cplex.SparsePair(ind=[name2idx[name] for name in ind], val=val),\
                         rhs=S - p_2[getNodeIndex(t)], sense="L")
            

# Création des variables pour le programme linéaire
# Si z = True, on ajoute la variable z et elle remplace les x_ij dans l'objectif
# Si dual = True, on ajoute les variables duales
def createVariables(cpl, z = False, dual = False):
    global name2idx # dictionnaire nom de variables => indice (pour accélerer le chargement des contraintes)
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



# Ajout des contraintes de flot dans le programme linéaire
def flowConstraints(cpl):
    # Contraintes de flot
    c = [] # Liste de listes de variables
    values = [] # Liste de listes des coefficients associés aux variables
    rhs = [] # Liste de listes de membres droits
    senses = "" # Chaîne qui représente le sens de l'inégalité
    
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


# Ajout des contraintes duales dans le programme linéaire
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
    

# Ajoute la coupe liée à l'objectif à partir des delta_1 calculés
def addObjectiveCut(cpl, delta_1):
    ind = ['z'] + [getVarName('x',m) for m in Mat]
    val = [1] + [-Mat[i][2] * (1+delta_1[i]) for i in range(len(Mat))]
    cpl.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[name2idx[name] for name in ind], val=val)],\
                                rhs=[0], senses="G")
    
# Ajoute la coupe liée au poids à partir des delta_2 calculés
def addWeightCut(cpl, delta_2):
    def getP(d):
        return p[getNodeIndex(d)] + delta_2[getNodeIndex(d)] * ph[getNodeIndex(d)]
    ind = [getVarName('x',m) for m in Mat]
    val = [getP(m[0])for m in Mat]
    cpl.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=[name2idx[name] for name in ind], val=val)],\
                                rhs=[S - getP(t)], senses="L")

# Résolution de la valeur robuste du chemin par un programme linéaire
# Note : cette fonction n'est pas utilisée
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

# Calcul de la valeur robuste du chemin associé à la solution x* par un PL
# Note: cette fonction n'est plus utilisée
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
    
# Calcul de la valeur robuste du chemin associé à la solution x* par un algo exact
def worstObjective_heu(x_star):
    # On récupère les arcs sélectionnés par x*
    m_star = []
    for i in range(len(Mat)):
        if x_star[i] > 0:
            m_star.append(Mat[i])
    # On les trie par distances décroissantes
    m_star.sort(key=lambda x: x[2], reverse=True)
    obj_val = sum([m[2] for m in m_star])
    delta_1_star = [0] * len(Mat)
    d = d1
    i = 0
    # Tant qu'on a pas distribué d1 ou que tous les arcs ont leur D_ij saturé
    # On ajoute le maximum possible sur les arcs de plus grande distance
    while d > 0 and i < len(m_star):
        delta = min([d,m_star[i][3]])
        d -= delta
        delta_1_star[Mat.index(m_star[i])] = delta
        obj_val += m_star[i][2] * delta
        i += 1
    
    return round(obj_val,5), delta_1_star

        
# Résolution du poids robuste du chemin par un programme linéaire
# Note: cette fonction n'est plus utilisée
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
    
# Calcul du poids robuste du chemin associé à la solution x* par un PL
# Note: cette fonction n'est plus utilisée
def worstWeight2(x_star):
    y = [0] * n
    for i in range(len(Mat)):
        if x_star[i] > 0:
            y[getNodeIndex(Mat[i][0])] += 1
    con_val, delta_2_star = worstWeight(y)
    return con_val, delta_2_star

# Calcul de la valeur robuste du chemin associé à la solution x* par un algo exact
def worstWeight_heu(x_star):
    delta_2_star = [0] * n
    # On récupère les noeuds traversés par le chemin associé à x*
    y = [t]
    for i in range(len(Mat)):
        if x_star[i] > 0:
            y.append(Mat[i][0])
    # On trie les noeuds par p chapeau décroissants
    y.sort(key=lambda x: ph[getNodeIndex(x)], reverse = True)
    con_val = sum(p[getNodeIndex(x)] for x in y)
    d = d2
    i = 0
    # On distribue d2 sur tous les noeuds au maximum dans l'ordre des p chapeau décroissants
    while d > 0 and i < len(y):
        delta = min([d, 2])
        d -= delta
        delta_2_star[getNodeIndex(y[i])] = delta
        con_val += delta * ph[getNodeIndex(y[i])]
        i += 1
    return con_val, delta_2_star

# Résolution du problème par plans coupants
def cuttingPlanes():
    cpl = cplex.Cplex()    
    cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    # Création des variables
    createVariables(cpl, z=True)
    # Ajout des contraintes de flot
    flowConstraints(cpl)
    
    # Ajout des premières coupes U0 et U1
    sum_dist = sum([sum(d) for d in Dist])
    addObjectiveCut(cpl, [min(d1*m[2]/sum_dist, m[3]) for m in Mat])
    
    sum_ph = sum(ph)
    addWeightCut(cpl, [min(ph[i]*d2/sum_ph, 2) for i in range(n)])
    
    is_opt = False
    # On boucle tant que l'on trouve une coupe liée à l'objectif ou au poids
    while not is_opt:
        # On résout le problème maître
        cpl.solve()
        
        # Enregistrement des résultats du problème maître
        x_star = []
        z = round(cpl.solution.get_values('z'),5)
        for m in Mat:
            x_star.append(cpl.solution.get_values(getVarName('x',m)))

        # Sous-problème lié à l'objectif
        obj_val, delta_1_star = worstObjective_heu(x_star)
        if obj_val > z:
            # On ajoute une coupe
            addObjectiveCut(cpl, delta_1_star)
        else:
            is_opt = True
        
            # Sous-problème lié au poids
            con_val, delta_2_star = worstWeight_heu(x_star)
            if con_val > S:
                # On ajoute une coupe
                addWeightCut(cpl, delta_2_star)
                is_opt = False

    return cpl.solution.get_objective_value(), cpl.solution.get_values()

# Résolution du problème par branch & cut
def branchAndCut():
    cpl = cplex.Cplex()
    
    cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    # On ajoute les variables
    createVariables(cpl, z=True)
    # On ajoute les contraintes de flot
    flowConstraints(cpl)
    
    # On ajoute les coupes initiales U0 et U1
    sum_dist = sum([sum(d) for d in Dist])
    addObjectiveCut(cpl, [min(d1*m[2]/sum_dist, m[3]) for m in Mat])
    sum_ph = sum(ph)
    addWeightCut(cpl, [min(ph[i]*d2/sum_ph, 2) for i in range(n)])
    
    # On ajoute le callback
    cpl.register_callback(CuttingPlaneCallback)
    # On résout le problème
    cpl.solve()

    return cpl.solution.get_objective_value(), cpl.solution.get_values()

# Résolution du problème dual
def dualSolve():
    cpl = cplex.Cplex()
    
    cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    # On ajoute les variables, dont les duales
    createVariables(cpl, dual=True)
    # On ajoute les contraintes de flot
    flowConstraints(cpl)
    # On ajoute les contraintes duales
    dualConstraints(cpl)
    
    # On résout le problème
    cpl.solve()
    
    return cpl.solution.get_objective_value(), cpl.solution.get_values()

# Résolution du problème dans sa version non robuste
def nonRobustSolve():
    cpl = cplex.Cplex()
    
    cpl.set_log_stream(None)
    cpl.set_results_stream(None)
    # On ajoute les variables
    createVariables(cpl, z=True)
    # On ajoute les contraintes de flot
    flowConstraints(cpl)
    
    # On ajoute la contrainte de poids
    addWeightCut(cpl, [0]*n)
    
    # On résout le problème
    cpl.solve()
    
    return cpl.solution.get_objective_value(), cpl.solution.get_values()

# Résolution heuristique du problème        
def heuristicSolve():
    # Retourne les n plus grands éléments de la concaténation des deux tableaux
    # où n est la taille du premier tableau
    def getMax(array1, array2):
        concat = np.concatenate((array1, array2))
        sort = sorted(concat, reverse=True)
        return np.array(sort[:len(array1)])
    
    # Calcule la distance dans le pire des cas pour les d_ij et les D_ij associés donnés
    def getWorstDist(dist_array, D_array):
        dist = 0.0
        d = d1
        sort_ind = sorted(range(len(dist_array)), key=dist_array.__getitem__, reverse=True)
        for ind in sort_ind:
            delta = min(D_array[ind], d)
            dist += dist_array[ind] * (1+delta)
            d -= delta
        return dist
        
    # Création du matrice nxn où ArkMap_ij est l'indice de l'arc ij dans le tableau Mat
    ArkMap = [[-1]*n for _ in range(n)]
    for i in range(len(Mat)):
        ArkMap[getNodeIndex(Mat[i][0])][getNodeIndex(Mat[i][1])] = i
    
    # Constantes
    ph_v = np.array([2.0] * int(d2/2) + ([1.0] if d2%2 == 1 else []))

    # Déclaration des tableaux
    P = [] # noeuds sélectionnés
    Pb = [x for x in unique] # noeuds non sélectionnés
       
    worst_w = [float("inf")] * n # poids dans le pire des cas du plus court chemin du noeud vers t
    weight = [0.0] * n # poids non robuste de ce chemin
    ph_sup = [np.array([0] * len(ph_v))]*n # les p chapeau non nuls dans ce chemin
    pred = [-1] * n # prédecesseur du noeud


    # Initialisation
    ph_sup[getNodeIndex(t)] = getMax(ph_sup[getNodeIndex(t)], [ph[getNodeIndex(t)]])
    weight[getNodeIndex(t)] += p[getNodeIndex(t)]
    worst_w[getNodeIndex(t)] = weight[getNodeIndex(t)] + sum(ph_v*ph_sup[getNodeIndex(t)])

    # Dijkstra sur les poids des noeuds
    while len(P) < n:
        # On sélectionne le noeud non sélectionnés de plus petit worst_w
        a = sorted([x for x in Pb],key=lambda x: worst_w[getNodeIndex(x)])[0]
        P.append(a)
        Pb.remove(a)

        for b in Pb:
            # Pour tous les noeuds non sélectionnés qui sont reliés au noeud sélectionné
            if Dist[getNodeIndex(b)][getNodeIndex(a)] > 0:
                # On calcule le poids d'un chemin du noeud b vers t passant par a
                n_weight = weight[getNodeIndex(a)] + p[getNodeIndex(b)]
                n_ph_sup = getMax(ph_sup[getNodeIndex(a)], [ph[getNodeIndex(b)]]) 
                n_worst = n_weight + sum(n_ph_sup*ph_v)
                # Si ce poids est meilleur que le poids actuel d'un chemin de b vers t
                # on met à jour les valeurs
                if n_worst < worst_w[getNodeIndex(b)]:
                    worst_w[getNodeIndex(b)] = n_worst
                    ph_sup[getNodeIndex(b)] = n_ph_sup
                    weight[getNodeIndex(b)] = n_weight
                    pred[getNodeIndex(b)] = a
                        
        
    # Déclaration des tableaux pour le second Dijkstra
    dist = [float('inf')] * n # plus courte distance dans le pire des cas d'un chemin de s vers ce noeud
    path_d = [np.array([])] * n # les d_ij des arcs traversés par ce chemin
    path_D = [np.array([])] * n # les D_ij associés
    weight_d = [0.0] *n # le poids non robuste du chemin
    ph_inf = [np.array([0] * len(ph_v))]*n # les p chapeau associés à des coefficients non nuls pour ce chemin
    pred_d = [-1] * n # le prédecesseur du noeud
    

    P = [] # noeuds sélectionnés
    Pb = [x for x in unique] # noeuds non sélectionnés
    
    # Initialisation
    dist[getNodeIndex(s)] = 0.0
    weight_d[getNodeIndex(s)] += p[getNodeIndex(s)]
    ph_inf[getNodeIndex(s)] = getMax(ph_inf[getNodeIndex(s)], [ph[getNodeIndex(s)]])
    
    while len(P) < n:
        # On sélectionne le noeud non encore sélectionné de plus petite distance
        a = sorted([x for x in Pb],key=lambda x: dist[getNodeIndex(x)])[0]
        P.append(a)
        Pb.remove(a)
        
        if a == t:
            break
        
        for b in Pb:
            # Pour tous les noeuds non sélectionnés qui sont reliés au noeud sélectionné
            if Dist[getNodeIndex(a)][getNodeIndex(b)] > 0: #and weight_d[getNodeIndex(a)] + worst_w[getNodeIndex(b)] <= S:
                # On calcule la distance du chemin de s vers b passant par a
                n_path_d = np.append(path_d[getNodeIndex(a)], Dist[getNodeIndex(a)][getNodeIndex(b)])
                n_path_D = np.append(path_D[getNodeIndex(a)], Delta[getNodeIndex(a)][getNodeIndex(b)])
                n_dist = getWorstDist(n_path_d, n_path_D)
                # On calcule le poids minimal du chemin de s vers t passant par a puis b
                # à partir des résultats du premier Dijkstra
                worst_ph = getMax(ph_inf[getNodeIndex(a)], ph_sup[getNodeIndex(b)])
                worst_weight = weight_d[getNodeIndex(a)] + weight[getNodeIndex(b)] + sum(ph_v*worst_ph)
                # Si la contrainte de poids est satisfaite et que la distance est strictement meilleure
                if worst_weight <= S and dist[getNodeIndex(b)] > n_dist:
                    # On met à jour les valeurs
                    dist[getNodeIndex(b)] = n_dist
                    path_d[getNodeIndex(b)] = n_path_d
                    path_D[getNodeIndex(b)] = n_path_D
                    weight_d[getNodeIndex(b)] = weight_d[getNodeIndex(a)] + p[getNodeIndex(b)]
                    ph_inf[getNodeIndex(b)] = getMax(ph_inf[getNodeIndex(a)], [ph[getNodeIndex(b)]])
                    pred_d[getNodeIndex(b)] = a
    
    # On recrée le chemin solution en remontant les prédecesseurs à partir de t
    curr = t
    path = []
    values = [0] * len(Mat)
    try:
        while curr != s:
            #print curr
            path = [curr] + path
            values[ArkMap[getNodeIndex(pred_d[getNodeIndex(curr)])][getNodeIndex(curr)]] = 1
            curr = pred_d[getNodeIndex(curr)]
        path = [s] + path
            
        obj_val, _ = worstObjective_heu(values)
        if obj_val != round(dist[getNodeIndex(t)],5):
            obj_val = -1
        else:
            con_val, delta_2_star = worstWeight_heu(values)
            if con_val > S:
                obj_val = 0
        return obj_val, values
    except:
        return 0,[]




###

#modes = ["CP", "dual", "b&c", "classic", "heuristic"]
modes = ["dual"] # Choix des modes de résolution

results = {mode : [] for mode in modes} # stockage des résultats

# Boucle sur les fichiers de données
for f in files.files[:]:
    # Boucle sur les modes de résolution
    for mode in modes:
        start_time = time.clock()
        print f
        
        # Lecture des données et pré-traitement
        readDataFromFile(f)
        
        data_time = time.clock()
        
        # Résolution du problème
        if mode == "CP":
            objectif, solution = cuttingPlanes()
        if mode == "dual":
            objectif, solution = dualSolve()
        if mode == "b&c":
            objectif, solution = branchAndCut()
        if mode == "classic":
            objectif, solution = nonRobustSolve()
        if mode == "heuristic":
            objectif, solution = heuristicSolve()
        
        print "\n"
        
        solve_time = time.clock() - data_time
        exec_time = time.clock() - start_time
        
        # Enregistrement des résultats
        # f est le fichier
        # objectif la valeur du chemin solution
        # path le chemin solution
        # exec_time le temps de résolution total
        # solve_time le temps de résolution sans la lecture et le pré-traitement des données
        if objectif > 0:
            path = printSolution(objectif, solution)
        results[mode].append((f, objectif, path, exec_time, solve_time))
        
        print "time :", exec_time, "s"
        print "dont data :", data_time - start_time
        print "dont solving :", solve_time
        print "\n\n"
        
# Affichage des résultats
for i in range(len(results[modes[0]])):
    for mode in modes:
        print results[mode][i][1], "(", results[mode][i][3], "s., solving", results[mode][i][4], ",s.)"
        







































