# -*- coding: utf-8 -*-

# Functions to compute the maximum productivity of a
# flexible Net integrating the E. coli metabolic network,
# citramalate synthesis and bioreactor variables.

from __future__ import division, print_function
import numpy as np
import cobra
from cobra import Reaction, Metabolite
from fnyzer import FNFactory, cobra2fn
# from CHOBiorFN_medium import EcoMinimizeMediumim
import os
import pandas as pd
from pandas_ods_reader import read_ods
from cobra import Metabolite, Reaction, Gene
import time
import argparse

# Constants
glucose_mmass = 180.1577 # Glucose molar mass (g/mol)
gln_mmass = 146.14
def_dilution_rt = 0.0166 
def_biomass_ini = 3.18727272727273 
def_biomass_fin = 3.18727272727273 
def_dataset = 'HP'
def_solver = 'cplex_direct'


# Parse parameters
parser = argparse.ArgumentParser()
parser.add_argument("-D", "--dilution_rate", default=def_dilution_rt, help="Dilution rate (h-1)")
parser.add_argument("-X_ini", "--biomass_ini", default=def_biomass_ini, help="Biomass initial (gdcw L-1)")
parser.add_argument("-X_fin", "--biomass_fin", default=def_biomass_fin, help="Biomass final (gdcw L-1)")
parser.add_argument("-data", "--dataset", default=def_dataset, help="Dataset used")
parser.add_argument("-solver", "--solver", default=def_solver, help="Solver used")
args = parser.parse_args()

D = float(args.dilution_rate)
X0 = float(args.biomass_ini)
Xf = float(args.biomass_fin)
data = args.dataset
solver = args.solver


if data == 'HP':
    lim_exch_phe = 0.0065*1.1 
    lim_exch_glc = 0.271*1.1
    lim_exch_gln = 0.0028*1.1
    lim_exch_asn = 0.075*1.1 
    lim_exch_ser = 0.04*1.1
    lim_exch_his = 0.0065*1.1
    lim_exch_thr = 0.015*1.1
    lim_exch_arg = 0.011*1.1
    lim_exch_tyr = 0.0053*1.1
    lim_exch_val = 0.015*1.1
    lim_exch_met = 0.0048*1.1
    lim_exch_trp = 0.0032*1.1
    lim_exch_ile = 0.012*1.1 
    lim_exch_leu = 0.02*1.1
    lim_exch_lys = 0.0126*1.1
    lim_exch_pro = 0.0116*1.1 
    lim_exch_asp = 0.025*1.1
    lim_exch_cys = 0.005*1.1
	
elif data == 'Late_Exp':
    lim_exch_phe = 0.0034*1.1 
    lim_exch_glc = 0.182*1.1
    lim_exch_gln = 0.0384*1.1
    lim_exch_asn = 0.022*1.1 
    lim_exch_ser = 0.019*1.1
    lim_exch_his = 0.00214*1.1 
    lim_exch_thr = 0.0059*1.1 
    lim_exch_arg = 0.00664*1.1
    lim_exch_tyr = 0.0045*1.1
    lim_exch_val = 0.00684*1.1 
    lim_exch_met = 0.00261*1.1 
    lim_exch_trp = 0.0032*1.1
    lim_exch_ile = 0.00458*1.1 
    lim_exch_leu = 0.00797*1.1 
    lim_exch_lys = 0.00714*1.1 
    lim_exch_pro = 0.00536*1.1 
    lim_exch_asp = 0.0035*1.1 
    lim_exch_cys = 0.00227*1.1 
    
if solver == 'cplex_direct':
    solver = 'cplex_direct'
    
elif solver == 'gurobi':
    solver = 'gurobi'



def comProductivity(fnet, abmw, c, glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu, D, X0, Xf):
    # Computes the antibody maximal production for a given dilution rate, biomass and the concentrations of the metabolites in the medium

    solution_list = [] 
   
    nut = []
    aa = ['glc', 'gln', 'phe', 'arg', 'asn', 'asp', 'cys', 'his', 'ile', 'leu', 'lys', 'met', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'glu']
    for i in aa:
        j = i.capitalize()
        nut.append(j)
   
    st = time.time()    

    genCHOBiorFN(fnet, D, glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu, Xmin = X0, Xmax = Xf)

    try:
		        netobj = FNFactory(fnet) # Build net object
		        netobj.optimize()
		        solution = netobj.objval
		        biomass =  netobj.places['X'].avm
		        print("  Growth with Xf in [", f"{Xf:.2f}", "] gdcw L-1:")
		        # for k in nut:
		            # print(k + " in tank: ", netobj.places[k].avm, 'mM')
		        
		        print("Solution:", solution, 'mM h-1') 
		        et = time.time()
		        elapsed_time = et - st
		        # print('Execution time:', elapsed_time, 'seconds')
    except:
    
        print('Did you install a solver? If not, please see Solvers section in the README')
        pass        
		                  
    return solution


def loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = solver):
    #### Parameters
    # filename: file with the SBML model
    # name: name of the FN
    # solver: solver to be used
   
    CHOcobramodel = cobra.io.read_sbml_model(filename+'.xml')
    
    """
    
    Add the antibodies to the model
    
    """
    
    antibody = read_ods("added_ab_reactions.ods")
    model_name = "iCHOv1.xml"


    ab_names = ["antiCD20", "iggM1", "iggM2", "iggM3", "mAb", "mAb2", "herceptin_subunit", "Igg_final_subunit"]
    for ab in ab_names:
        reaction = cobra.Reaction(ab)
        reaction.name = "%s synthesis" % ab
        reaction.subsystem = "PROTEIN PRODUCTION"
        reaction.lower_bound = 0
        reaction.upper_bound = 1000
        
        ab_met = cobra.Metabolite(id = ab,  # create a metabolite for each antibody and add it to the model
            name= ab,
            compartment='c')
        
        CHOcobramodel.add_metabolites(ab_met)

        r_dict = {}
        for i in range(len(antibody)):
                
            met = CHOcobramodel.metabolites.get_by_id(antibody.loc[i]["name"][2:])  # create metabolites for the reaction of each antibody based on the added_ab_reactions (3)dup.ods
            r_dict[met] = round(antibody.loc[i]["%s_coef" % ab] , 3)
            ab_dict = {ab_met: 1}  # adding the antibody to the products
                
                
        # print(ab_dict)
        r_dict.update(ab_dict)
        # print(r_dict)
        reaction.add_metabolites(r_dict)
        CHOcobramodel.add_reactions([reaction])
        
        # Exchange reaction for every antibody except for herceptin subunit (we need 2 subunits for the whole antibody)   
        
        if ab != 'Igg_final_subunit':   
        
            reaccEX = Reaction('EX_'+ ab)
            reaccEX.name = 'Exchange reaction to allow ' + ab + ' to leave the system'
            reaccEX.lower_bound = 0.0
            reaccEX.upper_bound = 1000.0                         
            reaccEX.add_metabolites({CHOcobramodel.metabolites.get_by_id(ab): -1.0})                         
            CHOcobramodel.add_reactions([reaccEX]) 
    
    # Herceptin formation reaction combines two herceptin subunits to produce the whole herceptin antibody
    
    reaction = Reaction('Igg_final_form')
    reaction.name = " IgG formation" 
    reaction.subsystem = "PROTEIN PRODUCTION"
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    
    
    igg_met = Metabolite(id = 'IgG',
            name= 'IgG',
            compartment='c')
    igg_subunit_met = CHOcobramodel.metabolites.get_by_id('Igg_final_subunit')
    
    reaction.add_metabolites({igg_met: 1, 
                              igg_subunit_met: -2})
    CHOcobramodel.add_reactions([reaction])    

    # Create the herceptin exchange reaction so the antibody is allowed to leave the system and it reaches steady state
    
    reaccEX = Reaction('EX_'+ 'IgG')
    reaccEX.name = 'Exchange reaction to allow ' + 'IgG' + ' to leave the system'
    reaccEX.lower_bound = 0.0
    reaccEX.upper_bound = 1000.0                               
    reaccEX.add_metabolites({CHOcobramodel.metabolites.get_by_id('IgG'): -1.0})                          
    CHOcobramodel.add_reactions([reaccEX]) 
    
    cobra.io.write_sbml_model(CHOcobramodel, model_name.replace(".xml", "_mapping_last.xml"))
    
    """
    
    Modify the exchange reactions of the model. 
    This is necessary beacause, otherwise, the iCHOv1_final model is too constrained and infeasible solutions are obtained.
    Only the upper bound for the backward reaction of the exchange reaction is considered. This adjusts the gluccose entrance flux to 0,1979 mmol_per_gDW_per_hr
    which corresponds with the measurements.
    
    """
    glc_reac = 'EX_glc__D_e'
    gln_reac = 'EX_gln__L_e'
    phe_reac = 'EX_phe__L_e'
    exch = list(CHOcobramodel.exchanges)
    for i in exch: 
        if i.lower_bound != 0.0:
                i.lower_bound = -1000
                # print(i.id, i.upper_bound, i.lower_bound)
        
    CHOcobramodel.reactions.get_by_id('EX_glu__L_e').lower_bound = -1000  # Glutamate exchange in the model is set to be only the forward reaction, we need it in the backward sense of the reaction
    
    fnet = cobra2fn(CHOcobramodel) # Build Flexible Net
    fnet['name'] =  name
    fnet['solver'] = solver
    
    return fnet


def genCHOBiorFN(fnet, D, glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu, Xmin, Xmax):
    # Adds to fnet the bioreactor dynamics
    #### Parameters, variables and units
    # D (h-1) Dilution rate
    # glcmM (mM) Glucose concentration in medium in mM
    # Xmin (gdcw L-1) Minimum density of cells in tank
    # Xmax (gdcw L-1) Maximum density of cells in tank
    # X (gdcw L-1) Density of cells in the tank 
    # G (mM) Concentration of glucose in the tank 
    # C (mM) Concentration of citramalate in the tank 
    # Fluxes are expressed in concentration of the reactant or product per hour
    
    
    aa = ['glc', 'gln', 'phe', 'arg', 'asn', 'asp', 'cys', 'his', 'ile', 'leu', 'lys', 'met', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'glu']
    aa_conc = [glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu]
    
    ### Bioreactor(tank) concentrations
    fnet['places']['X'] = 1.9512
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
    fnet['places']['A'] = 0.0


    fnet['shandlers'] = {}
            
        
    
    for i, k in zip(aa, aa_conc):
        if k != 0:
            j = i.capitalize()
            # print(j)
            fnet['places'][j] = 0.0
        
            ### Bioreactor reactions (transitions and handlers)
            
            # Glucose feed
            fnet['trans']['t' + i + 'in'] = {'l0': D*k, 'a0': 0}
            fnet['vhandlers']['v' + i + 'in'] = [{'a':('v' + i + 'in', j), 'v':('t' + i + 'in','v' + i + 'in')},
                                    'a == v']

            # Glucose uptake (from the tank into the cell)
            fnet['trans']['t' + i + 't'] = {'l0': 0, 'a0': 0}
            fnet['vhandlers']['v' + i + 't'] = [{'a':(j,'v' + i + 't'), 'v':('t' + i + 't','v' + i + 't')},
                                    'a == v']
            
            if i == 'glc':
                nut_ex = 'EX_glc__D_e' # Reaction ID of glucose exchange in metabolic network
            else:
                nut_ex = 'EX_' + i + '__L_e'
            
            nut_in = 't_'+nut_ex+'_b'  # b is for backward reaction and f is for forward reaction
            fnet['trans'][nut_in]['l0'] = 0 # Glucose uptake determined by its intensity handler
            nut_out = 't_'+nut_ex+'_f'
            fnet['trans'][nut_out]['l0'] = 0 # No glucose out of the cell allowed
            
            if i == 'phe':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_phe)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'glc':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_glc)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'gln':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_gln)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'asn':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asn)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'ser':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ser)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
           
            elif i == 'his':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_his)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'thr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_thr)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'arg':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_arg)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'tyr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_tyr)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'val':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_val)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'met':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_met)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'trp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_trp)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'ile':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ile)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'leu':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,   i + 'g <= ' + str(lim_exch_leu)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'lys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_lys)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'pro':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_pro)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'asp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asp)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            elif i == 'cys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_cys)] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]
            
            
            
            else: 
            
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax)] #, 'ug <= 0.1979'] # This ug is the upper bound for glucose exchange [mmol_per_gDW_per_hr]

            # Glucose out of the tank (to effluent)
            fnet['trans']['t' + i + 'out'] = {'l0': 0, 'a0': 0}
            fnet['vhandlers']['v' + i + 'out'] = [{'a':(j,'v' + i + 'out'), 'v':('t' +  i + 'out','v' + i + 'out')},
                                    'a == v']
            fnet['shandlers']['s' + i + 'out'] = [{'g':(j,'s' + i + 'out'), 'r':('s' + i + 'out','t' + i + 'out')},
                                    'r == g*'+str(D)]


# ----------------------------
    



    # Cell growth
    fnet['trans']['txt'] = {'l0': 0, 'a0': 0}  
    fnet['vhandlers']['vxt'] = [{'a':('vxt','X'), 'v':('txt','vxt')},
                                'a == v']
    biomass_reaction = 'BIOMASS_cho' # Reaction ID of biomass in metabolic network
    tBiomass = 't_'+biomass_reaction+ '_f' # transition for biomass
    fnet['trans'][tBiomass]['l0'] = 0 # Biomass production (i.e. growth rate) determined by its intensity handler. It should be equal to D in steady state.
    fnet['shandlers']['hr'] = [{'r':('hr',tBiomass), 'rt':('hr','txt')},
                                'r*'+str(Xmin)+'<= rt', 'rt <= r*'+str(Xmax)]

    # Cells out of the tank (to effluent)
    fnet['trans']['txout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vxout'] = [{'a':('X','vxout'), 'v':('txout','vxout')},
                                  'a == v']
    fnet['shandlers']['sxout'] = [{'x':('X','sxout'), 'r':('sxout','txout')},
                                   'r == x*'+str(D)]
    
    
    # Antibody production
    
    # Antibody production (from cell to tank)
    
    fnet['trans']['tat'] = {'l0': 0, 'a0': 0} 
    fnet['vhandlers']['vat'] = [{'a':('vat','A'), 'v':('tat','vat')},
                                'a == v']
    tExAb = 't_EX_IgG_f'
    fnet['trans'][tExAb]['l0'] = 0 # Citramalate production determined by its intensity handler.
    fnet['shandlers']['ha'] = [{'z':('ha',tExAb), 'at':('ha','tat')},
                                'z*'+str(Xmin)+'<= at', 'at <= z*'+str(Xmax)]

    # Antibody out of the tank (to effluent)
    fnet['trans']['taout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vaout'] = [{'a':('A','vaout'), 'v':('taout','vaout')},
                                  'a == v']
    fnet['shandlers']['saout'] = [{'a':('A','saout'), 'r':('saout','taout')},
                                   'r == a*'+str(D)]    

    fnet['obj'] = {'f': "avl['tat']", 'sense': 'max'}
    fnet['extrans'] = 'all'
    
    # --------------------
    
    cap_aa = []
    for i, k in zip(aa, aa_conc):
        if k != 0:
            j = i.capitalize()
            cap_aa.append(j)
            # print(j)
        
    fnet['actavplaces'] = ['X', 'A'] + cap_aa
    fnet['actplaces'] = ['X', 'A'] + cap_aa
    fnet['options'] = {
            'antype': 'cst',
            'savenet': False,            
            'printres': False,
            'printmodel': False,
            'writevars': {
                'avm': ['X', 'A'] + cap_aa,
                'avl':'all'},
            'plotres': False,
            'writexls': True,
            }

#### Uncomment next lines

fnet = loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = solver)

solution = comProductivity(fnet = fnet, abmw = 145881, c = 300e-12, glcmM = 35.15619796, gln = 7.872825297, phe = 1.275823895, arg = 1.946111508, asn = 5.122945878, asp = 1.354860924, cys = 0.379, his = 1.216191003, ile = 2.944437965, leu = 3.964891601, lys = 2.380196919, met = 1.010750472, pro = 4.621036929, ser = 4.956440488, thr = 2.774941185, trp = 0.926515767, tyr = 0.978250505, val = 2.946237191, glu = 1.84340765248991, D = D, X0 = X0, Xf = Xf) 

# EcoMinimizeMediumim(fnet=fnet, D=D, Xmin = X0, Xmax = Xf, ab_fixed = solution)
