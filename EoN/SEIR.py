#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 14:53:47 2020

@author: Justin Sheen

This code is a partial remix of Joel C. Miller's code for an SEIR epidemic. I
include the simulation of social distancing to the mix. I also include the idea
of expanding your quarantine circle after some amount of time.

Joel C. Miller's original code can be found here: https://epidemicsonnetworks.readthedocs.io/en/latest/functions/EoN.Gillespie_simple_contagion.html#EoN.Gillespie_simple_contagion
"""

import EoN
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
import random
import scipy
import numpy as np

"""
Initial set up of of the original network as well as the social distancing 
network.
The original network is made up of two graphs: one for local connections, and
the other for long distance connections.
The social distancing network is built out of the original network.
"""
random.seed(0)
philly_covid_direc = "/Users/Justin/"

# Local connections -----------------------------------------------------------
N = 1000
local_raw = np.genfromtxt(philly_covid_direc + "Philly_Covid/example_network.csv", delimiter=',')
rows, cols = np.where(local_raw == 1)
edges = zip(rows.tolist(), cols.tolist())
local = nx.Graph()
local.add_edges_from(edges)
for node_num in list(range(1000)):
    local.add_node(node_num)

# Long distance connections ---------------------------------------------------
ld = nx.fast_gnp_random_graph(N, 0.02) # about 10000 edges out of possible 500000 possible edges

# Combine the local connections graph and long distance connections graph -----
O = nx.Graph()
O.add_nodes_from(local)
O.add_edges_from(local.edges())
O.add_edges_from(ld.edges())
if ((O.number_of_nodes() != N) | 
    O.number_of_edges() >= (local.number_of_edges() + ld.number_of_edges())):
    raise NameError("O is not correct.")
node_attribute_dict = {node: 0.5+random.random() for node in O.nodes()}
edge_attribute_dict = {edge: 0.5+random.random() for edge in O.edges()}
nx.set_node_attributes(O, values=node_attribute_dict, name='expose2infect_weight')
nx.set_edge_attributes(O, values=edge_attribute_dict, name='transmission_weight')

# Create the SD graph out of the O graph --------------------------------------
SD = nx.Graph()
SD.add_nodes_from(local)
SD.add_edges_from(local.edges())
edges_to_delete = random.sample(list(range(ld.number_of_edges())), round(ld.number_of_edges() / 2))
SD_ld_edges = list(ld.edges)
for index in sorted(edges_to_delete, reverse=True):
    del SD_ld_edges[index]
SD.add_edges_from(SD_ld_edges)
nx.set_node_attributes(SD, values=node_attribute_dict, name='expose2infect_weight')
nx.set_edge_attributes(SD, values=edge_attribute_dict, name='transmission_weight')
if ((SD.number_of_nodes() != N) |
    SD.number_of_edges() >= (local.number_of_edges() + ld.number_of_edges() - round(ld.number_of_edges() / 2))):
    raise NameError("SD is not correct.")
    
# Method to expand your quarantine circle -------------------------------------
def expand_local(local):
    ccs = list(nx.connected_component_subgraphs(local))
    ccs_to_expand = random.sample(list(range(len(ccs))), round(len(ccs)/ 2))
    #ccs_to_expand = random.sample(list(range(len(ccs))), round(len(ccs)))
    pairs = list()
    while (len(ccs_to_expand) > 0):
        first_island_dex = random.randrange(0, len(ccs_to_expand))
        first_island = ccs_to_expand.pop(first_island_dex)
        
        second_island_dex = random.randrange(0, len(ccs_to_expand))
        second_island = ccs_to_expand.pop(second_island_dex)
        
        pairs.append((first_island, second_island))
    
    new_local = nx.Graph()
    for pair in pairs:
        island_one = pair[0]
        island_two = pair[1]
        n_nodes = ccs[island_one].number_of_nodes() + ccs[island_two].number_of_nodes()
        
        temp = nx.Graph()
        temp.add_nodes_from(ccs[island_one])
        temp.add_nodes_from(ccs[island_two])
        
        node_names = dict()
        for node_dex in list(range(n_nodes)):
            node_names[node_dex] = list(temp.nodes())[node_dex]
        
        fused = nx.complete_graph(n_nodes)
        fused = nx.relabel_nodes(fused, node_names)
        
        new_local.add_nodes_from(fused.nodes())
        new_local.add_edges_from(fused.edges())
    
    ccs_not_to_expand = np.setdiff1d(list(range(len(ccs))), ccs_to_expand)
    for cc_num in ccs_not_to_expand:
        new_local.add_nodes_from(ccs[cc_num])
        new_local.add_edges_from(ccs[cc_num].edges())
    
    return new_local
expanded_local = expand_local(local.copy())

# Create expanded quarantine graph out of SD graph ----------------------------
# Assign edge attributes to new edges of expanded_local graph -----------------
expanded_edge_attribute_dict = {edge: 0.5+random.random() for edge in expanded_local.edges()}

# Update dictionary to included attributes for new edges ----------------------
edge_attribute_dict.update(expanded_edge_attribute_dict)

# Add the new edges to copy of SD graph ---------------------------------------
expanded_SD = SD.copy()
expanded_SD.add_edges_from(expanded_local.edges())

# Add node and edge attributes to expanded_SD graph ---------------------------
nx.set_node_attributes(expanded_SD, values=node_attribute_dict, name='expose2infect_weight')
nx.set_edge_attributes(expanded_SD, values=edge_attribute_dict, name='transmission_weight')

"""
Define simple ODE mass-action SEIR model with social distancing

def ode_model(z, t, beta, beta_two, sigma, gamma, SD_day):
    
    # Reference https://www.idmod.org/docs/hiv/model-seir.html
    
    S, E, I, R = z
    N = S + E + I + R
    if (t < SD_day):
        dSdt = -beta * S * I / N
        dEdt = beta * S * I / N- sigma * E
        dIdt = sigma * E - gamma * I
        dRdt = gamma * I
    else:
        dSdt = -beta_two * S / N * I
        dEdt = beta_two * S * I / N - sigma * E
        dIdt = sigma * E - gamma * I
        dRdt = gamma * I

    return [dSdt, dEdt, dIdt, dRdt]

def ode_solver(t, initial_conditions, params, SD_day):
    initE, initI, initR, initN = initial_conditions
    beta, beta_two, sigma, gamma = params
    initS = initN - (initE + initI + initR)
    res = scipy.integrate.odeint(ode_model, [initS, initE, initI, initR], t, args=(beta, beta_two, sigma, gamma, SD_day))
    return res

soln = ode_solver(t=list(range(0, 200, 1)),
                        initial_conditions=[0, 5, 0, N], 
                        params=[0.5, 0.25, 0.2, 0.167],
                        SD_day=50)

# Plot the SEIR ODE model dynamics --------------------------------------------
# plt.plot(list(range(0, 200, 1)), soln[:,0])
# plt.plot(list(range(0, 200, 1)), soln[:,1])
# plt.plot(list(range(0, 200, 1)), soln[:,2])
# plt.plot(list(range(0, 200, 1)), soln[:,3])

# Get ODE solution --------------------------------------------------------
soln = ode_solver(t=list(range(0, 200, 1)),
                    initial_conditions=[0, 5, 0, N],
                    params=[0.5, 0.25, 0.2, 0.167],
                    SD_day=SD_days[SD_day_dex])
soln_I = soln[:, 2]
"""

"""
Set the spontaneous parameters and transmission parameters for the SEIR simulation
"""

H = nx.DiGraph()
H.add_node('S')
H.add_edge('E', 'I', rate = 3* 1/5, weight_label='expose2infect_weight')
H.add_edge('I', 'R', rate = 3* 1/8)

J = nx.DiGraph()
J.add_edge(('I', 'S'), ('I', 'E'), rate = (2.5 / (O.number_of_edges() / N) / 5), weight_label='transmission_weight')
IC = defaultdict(lambda: 'S')
for node in range(3):
    IC[node] = 'I'

return_statuses = ('S', 'E', 'I', 'R')

"""
Run the SEIR epidemic on the O graph until the SD day, then switch to the SD 
graph, then switch to the expanding quarantine circle graph
"""
cumul_sims_SD = list()
cumul_sims_EQ = list()
cumul_sims_RL = list()
SD_days = list(range(2, 92, 2))
EQ_implemented = 21
for SD_day in SD_days:
    for iteration in list(range(5)):
        # First, run on the original graph ------------------------------------
        full_O = EoN.Gillespie_simple_contagion(O, H, J, IC, return_statuses, tmax = SD_day, return_full_data=True)
        t_O = full_O.t()
        S_O = full_O.S()
        E_O = full_O.summary()[1]['E']
        I_O = full_O.I()
        R_O = full_O.R()
        nodes_O_final = full_O.get_statuses(list(range(N)), t_O[-1])
        
        # Next, run on the SD graph -------------------------------------------
        SD_IC = defaultdict(lambda: 'S')
        for node in range(N):
            SD_IC[node] = nodes_O_final[node]
        full_SD = EoN.Gillespie_simple_contagion(SD, H, J, SD_IC, return_statuses, tmax = EQ_implemented, return_full_data=True)
        t_SD = full_SD.t()
        S_SD = full_SD.S()
        E_SD = full_SD.summary()[1]['E']
        I_SD = full_SD.I()
        R_SD = full_SD.R()
        nodes_SD_final = full_SD.get_statuses(list(range(N)), t_SD[-1])
        
        # Next, finish sim on the expanded quarantine graph, SD graph, or
        # relaxing lockdown ---------------------------------------------------
        EQ_IC = defaultdict(lambda: 'S')
        for node in range(N):
            EQ_IC[node] = nodes_SD_final[node]
        full_EQ = EoN.Gillespie_simple_contagion(expanded_SD, H, J, EQ_IC, return_statuses, tmax = float('Inf'), return_full_data=True)
        t_EQ = full_EQ.t()
        S_EQ = full_EQ.S()
        E_EQ = full_EQ.summary()[1]['E']
        I_EQ = full_EQ.I()
        R_EQ = full_EQ.R()
        
        FSD_IC = defaultdict(lambda: 'S')
        for node in range(N):
            FSD_IC[node] = nodes_SD_final[node]
        full_FSD = EoN.Gillespie_simple_contagion(SD, H, J, FSD_IC, return_statuses, tmax = float('Inf'), return_full_data=True)
        t_FSD = full_FSD.t()
        S_FSD = full_FSD.S()
        E_FSD = full_FSD.summary()[1]['E']
        I_FSD = full_FSD.I()
        R_FSD = full_FSD.R()
        
        RL_IC = defaultdict(lambda: 'S')
        for node in range(N):
            RL_IC[node] = nodes_SD_final[node]
        full_RL = EoN.Gillespie_simple_contagion(O, H, J, RL_IC, return_statuses, tmax = float('Inf'), return_full_data=True)
        t_RL = full_RL.t()
        S_RL = full_RL.S()
        E_RL = full_RL.summary()[1]['E']
        I_RL = full_RL.I()
        R_RL = full_RL.R()
        
        # Combine the time series of SD ---------------------------------------
        t = np.concatenate((t_O, (t_SD + t_O[-1])), axis=None)
        S = np.concatenate((S_O, S_SD), axis=None)
        E = np.concatenate((E_O, E_SD), axis=None)
        I = np.concatenate((I_O, I_SD), axis=None)
        R = np.concatenate((R_O, R_SD), axis=None)
        t = np.concatenate((t, (t_FSD + t_SD[-1] + t_O[-1])), axis = None)
        S = np.concatenate((S, S_FSD), axis=None)
        E = np.concatenate((E, E_FSD), axis=None)
        I = np.concatenate((I, I_FSD), axis=None)
        R = np.concatenate((R, R_FSD), axis=None)
    
        # Save SD results in list of lists ------------------------------------
        to_add = list()
        to_add.append(t)
        to_add.append(S)
        to_add.append(E)
        to_add.append(I)
        to_add.append(R)
        cumul_sims_SD.append(to_add)
        
        # Combine the time series of EQ ---------------------------------------
        t = np.concatenate((t_O, (t_SD + t_O[-1])), axis=None)
        S = np.concatenate((S_O, S_SD), axis=None)
        E = np.concatenate((E_O, E_SD), axis=None)
        I = np.concatenate((I_O, I_SD), axis=None)
        R = np.concatenate((R_O, R_SD), axis=None)
        t = np.concatenate((t, (t_EQ + t_SD[-1] + t_O[-1])), axis = None)
        S = np.concatenate((S, S_EQ), axis=None)
        E = np.concatenate((E, E_EQ), axis=None)
        I = np.concatenate((I, I_EQ), axis=None)
        R = np.concatenate((R, R_EQ), axis=None)
    
        # Save EQ results in list of lists ------------------------------------
        to_add = list()
        to_add.append(t)
        to_add.append(S)
        to_add.append(E)
        to_add.append(I)
        to_add.append(R)
        cumul_sims_EQ.append(to_add)
        
        # Combine the time series of RL ---------------------------------------
        t = np.concatenate((t_O, (t_SD + t_O[-1])), axis=None)
        S = np.concatenate((S_O, S_SD), axis=None)
        E = np.concatenate((E_O, E_SD), axis=None)
        I = np.concatenate((I_O, I_SD), axis=None)
        R = np.concatenate((R_O, R_SD), axis=None)
        t = np.concatenate((t, (t_RL + t_SD[-1] + t_O[-1])), axis = None)
        S = np.concatenate((S, S_RL), axis=None)
        E = np.concatenate((E, E_RL), axis=None)
        I = np.concatenate((I, I_RL), axis=None)
        R = np.concatenate((R, R_RL), axis=None)
    
        # Save RL results in list of lists ------------------------------------
        to_add = list()
        to_add.append(t)
        to_add.append(S)
        to_add.append(E)
        to_add.append(I)
        to_add.append(R)
        cumul_sims_RL.append(to_add)
        
"""
Plot comparison of the ODE model solution to the stochastic simulations of the 
network model
"""
SD_day_dex = 0
for sim_num in range(len(cumul_sims_SD)):
    # Plot the stochastic simulations together --------------------------------
    plt.xlim(0, 200)
    plt.ylim(0, 0.4)
    plt.ylabel('I/N')
    plt.xlabel('t')
    plt.title('SD: ' + str(SD_days[SD_day_dex]))
    #plt.plot(list(range(0, 200, 1)), soln_I, label="mass-action", color="red")
    plt.plot(cumul_sims_SD[sim_num][0], cumul_sims_SD[sim_num][3] / N, label='network', color="black")
    
    if (np.where(cumul_sims_EQ[sim_num][0] >= (SD_days[SD_day_dex] + EQ_implemented))[0].size != 0):
        first_EQ = np.where(cumul_sims_EQ[sim_num][0] >= (SD_days[SD_day_dex] + EQ_implemented))[0][0]
        plt.plot(cumul_sims_EQ[sim_num][0][first_EQ:-1], cumul_sims_EQ[sim_num][3][first_EQ:-1] / N, label='network', color="red")
    
    if (np.where(cumul_sims_RL[sim_num][0] >= (SD_days[SD_day_dex] + EQ_implemented))[0].size != 0):
        first_RL = np.where(cumul_sims_RL[sim_num][0] >= (SD_days[SD_day_dex] + EQ_implemented))[0][0]
        plt.plot(cumul_sims_RL[sim_num][0][first_RL:-1], cumul_sims_RL[sim_num][3][first_RL:-1] / N, label='network', color="green")
    
    if (np.mod(sim_num + 1, 5) == 0):
        plt.axvline(x = SD_days[SD_day_dex])
        plt.axvline(x = (SD_days[SD_day_dex] + EQ_implemented), linestyle='--')
        file_name = '/Users/Justin/Philly_Covid/EoN/res/SEIR_' + str(SD_days[SD_day_dex]) + '.png'
        plt.savefig(file_name)
        plt.clf()
        SD_day_dex += 1