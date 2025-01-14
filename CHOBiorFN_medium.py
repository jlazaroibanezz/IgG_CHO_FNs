# from CHOBiorFN_all_aa_nutrients_HP_pru1_def import loadCHOmodel

import numpy as np
import cobra
from cobra import Reaction, Metabolite
from fnyzer import FNFactory, cobra2fn
import os
import pandas as pd
from pandas_ods_reader import read_ods
from cobra import Metabolite, Reaction, Gene
import argparse
from CHOBiorFN import fnet, solution


# Constants
glucose_mmass = 180.1577 # Glucose molar mass (g/mol)
gln_mmass = 146.14
def_dilution_rt = 0.0166 
def_biomass_ini = 3.18727272727273 
def_biomass_fin = 3.18727272727273 
def_dataset = 'HP'

# Solver
solver = 'cplex_direct' # glpk

# Parse parameters
parser = argparse.ArgumentParser()
parser.add_argument("-D", "--dilution_rate", default=def_dilution_rt, help="Dilution rate (h-1)")
parser.add_argument("-X_ini", "--biomass_ini", default=def_biomass_ini, help="Biomass initial (gdcw L-1)")
parser.add_argument("-X_fin", "--biomass_fin", default=def_biomass_fin, help="Biomass final (gdcw L-1)")
parser.add_argument("-data", "--dataset", default=def_dataset, help="Dataset used")
args = parser.parse_args()

D = float(args.dilution_rate)
X0 = float(args.biomass_ini)
Xf = float(args.biomass_fin)
data = args.dataset


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


def EcoMinimizeMediumim(fnet, D, Xmin, Xmax, ab_fixed):
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
    
    # aa_conc = [glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu]
    aa_price = [0.92, 2.89, 1.615, 1.615, 1.535, 2.79, 1.615, 1.445, 1.445, 1.615, 1.7, 1.615, 1.615, 1.615, 1.445, 1.615, 0.558, 1.615, 2.21]
    
    ### Bioreactor(tank) concentrations
    fnet['places']['X'] = 1.9512
    #X 
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
    # fnet['places']['G'] = 0.0
    fnet['places']['A'] = 0.0
    fnet['places']['W'] = 0.0
    fnet['shandlers'] = {}
            
        
    
    for i in aa:
        # if k != 0:
            j = i.capitalize()
            # print(j)
            # print(type(i))
            fnet['places'][j] = 0.0
        
            ### Bioreactor reactions (transitions and handlers)
            

            # Glucose feed
            fnet['trans']['t' + i + 'in'] = {'l0': 0, 'a0': 0}
            fnet['vhandlers']['v' + i + 'in'] = [{'a':('v' + i + 'in', j), 'v':('t' + i + 'in','v' + i + 'in')},
                                    'a == v']

            # Glucose uptake (from the tank into the cell)
            fnet['trans']['t' + i + 't'] = {'l0': 0, 'a0': 0}
            fnet['vhandlers']['v' + i + 't'] = [{'a':(j,'v' + i + 't'), 'v':('t' + i + 't','v' + i + 't')},
                                    'a == v']
            
            fnet['shandlers']['s' + i + 'in'] = [{'n':(j, 's' + i + 'in'), 'o':('s' + i + 'in', 't' + i + 'in')},
                                   'o == n*'+str(D)]
            
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
    # glc_reac = 'EX_glc__D_e'
    tBiomass = 't_'+biomass_reaction+ '_f' # transition for biomass
    # tglcb = 't_'+glc_reac+ '_b'
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
                                'z*'+str(Xmin)+'<= at', 'at <= z*'+str(Xmax)] # , 'z >= 0.00087']

    # Antibody out of the tank (to effluent)
    fnet['trans']['taout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vaout'] = [{'a':('A','vaout'), 'v':('taout','vaout')},
                                  'a == v']
    fnet['shandlers']['saout'] = [{'a':('A','saout'), 'r':('saout','taout')},
                                   'r == a*'+str(D)]    

   

    # Cells out of the tank (to effluent)
    fnet['trans']['txout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vxout'] = [{'a':('X','vxout'), 'v':('txout','vxout')},
                                  'a == v']
    fnet['shandlers']['sxout'] = [{'x':('X','sxout'), 'r':('sxout','txout')},
                                   'r == x*'+str(D)]

    fnet['obj'] = {'f': "0.92*avl['tglcin'] + 2.89*avl['tglnin'] + 1.615*avl['tphein'] + 1.615*avl['targin'] + 1.535*avl['tasnin'] + 2.79*avl['taspin'] + 1.615*avl['tcysin'] + 1.445*avl['thisin'] + 1.445*avl['tilein'] + 1.615*avl['tleuin'] + 1.7*avl['tlysin'] + 1.615*avl['tmetin'] + 1.615*avl['tproin'] + 1.615*avl['tserin'] + 1.445*avl['tthrin'] + 1.615*avl['ttrpin'] + 0.558*avl['ttyrin'] + 1.615*avl['tvalin'] + 2.21*avl['tgluin']", 'sense': 'min'}
    
    
    fnet['extracons'] = ["avl['tat']== " + str(ab_fixed)]
    
    # --------------------
    
    cap_aa = []
    for i in aa:
            j = i.capitalize()
            cap_aa.append(j)
            
            
    netobj = FNFactory(fnet)
    netobj.optimize()
    ecominmediumim = netobj.objval
    for i in aa:
        j = i.capitalize()
        print(j + " enters in tank: ", netobj.trans['t' + str(i) + 'in'].avl/D, 'mM')
    print('The weighted sum for the economic minimization of the medium is ', ecominmediumim)
        
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
    return ecominmediumim            
       
# fnet = loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = solver)

# solution = comProductivity(fnet = fnet, abmw = 145881, c = 300e-12, glcmM = 35.15619796, gln = 7.872825297, phe = 1.275823895, arg = 1.946111508, asn = 5.122945878, asp = 1.354860924, cys = 0.379, his = 1.216191003, ile = 2.944437965, leu = 3.964891601, lys = 2.380196919, met = 1.010750472, pro = 4.621036929, ser = 4.956440488, thr = 2.774941185, trp = 0.926515767, tyr = 0.978250505, val = 2.946237191, glu = 1.84340765248991, D = D, X0 = X0, Xf = Xf) 

EcoMinimizeMediumim(fnet=fnet, D=D, Xmin = X0, Xmax = Xf, ab_fixed = solution)
