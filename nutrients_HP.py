# -*- coding: utf-8 -*-

# Functions to compute the maximum productivity of a
# flexible Net integrating the E. coli metabolic network,
# IgG synthesis and bioreactor variables.

from __future__ import division, print_function
import numpy as np
import cobra
from cobra import Reaction, Metabolite
from fnyzer import FNFactory, cobra2fn
import os
import pandas as pd
from pandas_ods_reader import read_ods
from cobra import Metabolite, Reaction, Gene

def loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = "gurobi"):
    #### Parameters
    # filename: file with the SBML model
    # name: name of the FN
    # solver: solver to be used
   
    CHOcobramodel = cobra.io.read_sbml_model(filename+'.xml')
    
    """
    
    Add the antibodies to the model
    
    """
    
    antibody = read_ods("added_ab_reactions (3)dup.ods") # Antibody reaction
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
            # r_dict[met] = antibody.loc[i]["%s_coef" % ab]
            ab_dict = {ab_met: 1}  # adding the antibody to the products
                
                
        print(ab_dict)
        r_dict.update(ab_dict)
        print(r_dict)
        reaction.add_metabolites(r_dict)
        CHOcobramodel.add_reactions([reaction])
        
        # Exchange reaction for every antibody (we need 2 subunits for the whole antibody)   
        
        if ab != 'Igg_final_subunit':   
        
            reaccEX = Reaction('EX_'+ ab)
            reaccEX.name = 'Exchange reaction to allow ' + ab + ' to leave the system'
            reaccEX.lower_bound = 0.0
            reaccEX.upper_bound = 1000.0                       
            reaccEX.add_metabolites({CHOcobramodel.metabolites.get_by_id(ab): -1.0})                        
            CHOcobramodel.add_reactions([reaccEX]) 

    
    # IgG formation reaction combines two IgG subunits to produce the whole herceptin antibody
    
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

    # Create the IgG exchange reaction so the antibody is allowed to leave the system and it reaches steady state
    
    reaccEX = Reaction('EX_'+ 'IgG')
    reaccEX.name = 'Exchange reaction to allow ' + 'IgG' + ' to leave the system'
    reaccEX.lower_bound = 0.0
    reaccEX.upper_bound = 1000.0                               
    reaccEX.add_metabolites({CHOcobramodel.metabolites.get_by_id('IgG'): -1.0})                          
    CHOcobramodel.add_reactions([reaccEX]) 
    
    cobra.io.write_sbml_model(CHOcobramodel, model_name.replace(".xml", "_mapping_last.xml"))
    
    """
    
    Modify the exchange reactions of the model to unconstrained
    
    """
    glc_reac = 'EX_glc__D_e'
    gln_reac = 'EX_gln__L_e'
    phe_reac = 'EX_phe__L_e'
    exch = list(CHOcobramodel.exchanges)
    for i in exch: 
        if i.lower_bound != 0.0:
                i.lower_bound = -1000
                print(i.id, i.upper_bound, i.lower_bound)
        
    CHOcobramodel.reactions.get_by_id('EX_glu__L_e').lower_bound = -1000  # Glutamate exchange in the model is set to be only the forward reaction, we need it in the backward sense of the reaction
    
    
    fnet = cobra2fn(CHOcobramodel) # Build Flexible Net
    fnet['name'] =  name
    fnet['solver'] = solver
    
    return fnet


def genCHOBiorFNmax(fnet, D, glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu, Xmin, Xmax, ab_fixed, aminosol): 
    # Adds to fnet the bioreactor dynamics
    #### Parameters, variables and units
    # D (h-1) Dilution rate
    # glcmM (mM) Glucose concentration in medium in mM
    # Xmin (gdcw L-1) Minimum density of cells in tank
    # Xmax (gdcw L-1) Maximum density of cells in tank
    # X (gdcw L-1) Density of cells in the tank 
    # G (mM) Concentration of glucose in the tank 
    # C (mM) Concentration of IgG in the tank 
    # Fluxes are expressed in concentration of the reactant or product per hour
    
    
    aa = ['glc', 'gln', 'phe', 'arg', 'asn', 'asp', 'cys', 'his', 'ile', 'leu', 'lys', 'met', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'glu']
    
    aa_conc = [glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu]


    # Adding the maximum uptake values
    lim_exch_phe = 0.0065*1.1  # 10% of uncertainty so it is not so restrictive
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

    
    ### Bioreactor(tank) concentrations
    fnet['places']['X'] = 1.9512
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
    fnet['places']['A'] = 0.0
    fnet['shandlers'] = {}
            
        
    
    for i, k in zip(aa, aa_conc):
        if k != 0:
            j = i.capitalize()
            print(j)
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
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_phe)] 
            
            elif i == 'glc':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_glc)] 
            
            elif i == 'gln':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_gln)] 
            
            elif i == 'asn':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asn)] 
            
            elif i == 'ser':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ser)] 
           
            elif i == 'his':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_his)] 
            
            elif i == 'thr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_thr)] 
            
            elif i == 'arg':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_arg)] 
            
            elif i == 'tyr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_tyr)] 
            
            elif i == 'val':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_val)] 
            
            elif i == 'met':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_met)] 
            
            elif i == 'trp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_trp)] 
            
            elif i == 'ile':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ile)] 
            
            elif i == 'leu':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,   i + 'g <= ' + str(lim_exch_leu)] 
            
            elif i == 'lys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_lys)] 
            
            elif i == 'pro':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_pro)] 
                
            elif i == 'asp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asp)] 
            
            elif i == 'cys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_cys)] 
            
            
            else: 
            
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax)] #, 'ug <= 0.1979'] 
                
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
    
    # Antibody production (from cell to tank)
    
    fnet['trans']['tat'] = {'l0': 0, 'a0': 0} 
    fnet['vhandlers']['vat'] = [{'a':('vat','A'), 'v':('tat','vat')},
                                'a == v']
    tExAb = 't_EX_IgG_f'
    fnet['trans'][tExAb]['l0'] = 0 # IgG production determined by its intensity handler.
    fnet['shandlers']['ha'] = [{'z':('ha',tExAb), 'at':('ha','tat')},
                                'z*'+str(Xmin)+'<= at', 'at <= z*'+str(Xmax)] # , 'z >= 0.00087']

    # Antibody out of the tank (to effluent)
    fnet['trans']['taout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vaout'] = [{'a':('A','vaout'), 'v':('taout','vaout')},
                                  'a == v']
    fnet['shandlers']['saout'] = [{'a':('A','saout'), 'r':('saout','taout')},
                                   'r == a*'+str(D)]    

    
    if aminosol == 'glc':
        fnet['obj'] = {'f': "avl['t_EX_glc__D_e_b']", 'sense': 'max'}
    
    else:
        fnet['obj'] = {'f': "avl['t_EX_" + aminosol +"__L_e_b']", 'sense': 'max'}
        

    fnet['extrans'] = 'all'
    
    # Constraints if we want to maximize the concentration of any nutrient in the tank. Antibody flux is fixed to its maximum computed before
    

    fnet['extracons'] = ["avl['tat']== " + str(ab_fixed)]
    
    cap_aa = []
    for i, k in zip(aa, aa_conc):
        if k != 0:
            j = i.capitalize()
            cap_aa.append(j)
        
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
    
    netobj = FNFactory(fnet)
    netobj.optimize()
    solutionmax = netobj.objval
    print(solutionmax)
    return solutionmax
    

def genCHOBiorFNmin(fnet, D, glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu, Xmin, Xmax, ab_fixed, aminosol):   #############
    # Adds to fnet the bioreactor dynamics
    #### Parameters, variables and units
    # D (h-1) Dilution rate
    # glcmM (mM) Glucose concentration in medium in mM
    # Xmin (gdcw L-1) Minimum density of cells in tank
    # Xmax (gdcw L-1) Maximum density of cells in tank
    # X (gdcw L-1) Density of cells in the tank 
    # G (mM) Concentration of glucose in the tank 
    # C (mM) Concentration of IgG in the tank 
    # Fluxes are expressed in concentration of the reactant or product per hour
    
    
    aa = ['glc', 'gln', 'phe', 'arg', 'asn', 'asp', 'cys', 'his', 'ile', 'leu', 'lys', 'met', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'glu']
    
    aa_conc = [glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu]


    # Adding the maximum uptake values
    lim_exch_phe = 0.0065*1.1  # 10% of uncertainty so it is not so restrictive
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
    
    
    ### Bioreactor(tank) concentrations
    fnet['places']['X'] = 1.9512
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
    fnet['places']['A'] = 0.0
    fnet['shandlers'] = {}
            
        
    
    for i, k in zip(aa, aa_conc):
        if k != 0:
            j = i.capitalize()
            print(j)
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
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_phe)] 
            
            elif i == 'glc':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_glc)] 
            
            elif i == 'gln':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_gln)] 
            
            elif i == 'asn':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asn)] 
            
            elif i == 'ser':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ser)] 
           
            elif i == 'his':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_his)] 
            
            elif i == 'thr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_thr)] 
            
            elif i == 'arg':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_arg)] 
            
            elif i == 'tyr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_tyr)] 
            
            elif i == 'val':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_val)] 
            
            elif i == 'met':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_met)] 
            
            elif i == 'trp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_trp)] 
            
            elif i == 'ile':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ile)] 
            
            elif i == 'leu':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,   i + 'g <= ' + str(lim_exch_leu)] 
            
            elif i == 'lys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_lys)] 
            
            elif i == 'pro':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_pro)] 
            
            elif i == 'asp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asp)] 
            
            elif i == 'cys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_cys)] 
                
            else: 
            
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax)] #, 'ug <= 0.1979']
                
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
   
    
    # Antibody production (from cell to tank)
    
    fnet['trans']['tat'] = {'l0': 0, 'a0': 0} 
    fnet['vhandlers']['vat'] = [{'a':('vat','A'), 'v':('tat','vat')},
                                'a == v']
    tExAb = 't_EX_IgG_f'
    fnet['trans'][tExAb]['l0'] = 0 # IgG production determined by its intensity handler.
    fnet['shandlers']['ha'] = [{'z':('ha',tExAb), 'at':('ha','tat')},
                                'z*'+str(Xmin)+'<= at', 'at <= z*'+str(Xmax)] # , 'z >= 0.00087']

    # Antibody out of the tank (to effluent)
    fnet['trans']['taout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vaout'] = [{'a':('A','vaout'), 'v':('taout','vaout')},
                                  'a == v']
    fnet['shandlers']['saout'] = [{'a':('A','saout'), 'r':('saout','taout')},
                                   'r == a*'+str(D)]    

    if aminosol == 'glc':
        fnet['obj'] = {'f': "avl['t_EX_glc__D_e_b']", 'sense': 'min'}
    
    else:
        fnet['obj'] = {'f': "avl['t_EX_" + aminosol +"__L_e_b']", 'sense': 'min'}
        
    fnet['extrans'] = 'all'
    
    # Constraints if we want to maximize the concentration of any nutrient in the tank. Antibody flux is fixed to its maximum computed before
    

    fnet['extracons'] = ["avl['tat']== " + str(ab_fixed)]

    
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
    netobj = FNFactory(fnet)
    netobj.optimize()
    solutionmin = netobj.objval
    print(solutionmin)
    return solutionmin


def MinimizeMediumim(fnet, D, Xmin, Xmax, ab_fixed, tau):  
    # Adds to fnet the bioreactor dynamics
    #### Parameters, variables and units
    # D (h-1) Dilution rate
    # glcmM (mM) Glucose concentration in medium in mM
    # Xmin (gdcw L-1) Minimum density of cells in tank
    # Xmax (gdcw L-1) Maximum density of cells in tank
    # X (gdcw L-1) Density of cells in the tank 
    # G (mM) Concentration of glucose in the tank 
    # C (mM) Concentration of IgG in the tank 
    # Fluxes are expressed in concentration of the reactant or product per hour
    
    
    aa = ['glc', 'gln', 'phe', 'arg', 'asn', 'asp', 'cys', 'his', 'ile', 'leu', 'lys', 'met', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'glu']
    
    aa_price = [0.92, 2.89, 1.615, 1.615, 1.535, 2.79, 1.615, 1.445, 1.445, 1.615, 1.7, 1.615, 1.615, 1.615, 1.445, 1.615, 0.558, 1.615, 2.21]

    # Adding the maximum uptake values
    lim_exch_phe = 0.0065*1.1  # 10% of uncertainty so it is not so restrictive
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
    
    ### Bioreactor(tank) concentrations
    fnet['places']['X'] = 1.9512
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
    fnet['places']['A'] = 0.0
    fnet['places']['W'] = 0.0
    fnet['shandlers'] = {}
            
        
    
    for i in aa:

            j = i.capitalize()
            print(j)
            fnet['places'][j] = 0.0
            

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
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_phe)] 
            
            elif i == 'glc':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_glc)] 
            
            elif i == 'gln':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_gln)] 
            
            elif i == 'asn':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asn)] 
            
            elif i == 'ser':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ser)] 
           
            elif i == 'his':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_his)] 
            
            elif i == 'thr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_thr)] 
            
            elif i == 'arg':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_arg)] 
            
            elif i == 'tyr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_tyr)] 
            
            elif i == 'val':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_val)] 
            
            elif i == 'met':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_met)] 
            
            elif i == 'trp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_trp)] 
            
            elif i == 'ile':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ile)] 
            
            elif i == 'leu':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,   i + 'g <= ' + str(lim_exch_leu)] 
            
            elif i == 'lys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_lys)] 
            
            elif i == 'pro':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_pro)] 
            
            elif i == 'asp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asp)] 
            
            elif i == 'cys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_cys)] 
            
            
            else: 
            
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax)] #, 'ug <= 0.1979'] 
                
            # Glucose out of the tank (to effluent)
            fnet['trans']['t' + i + 'out'] = {'l0': 0, 'a0': 0}
            fnet['vhandlers']['v' + i + 'out'] = [{'a':(j,'v' + i + 'out'), 'v':('t' +  i + 'out','v' + i + 'out')},
                                    'a == v']
            fnet['shandlers']['s' + i + 'out'] = [{'g':(j,'s' + i + 'out'), 'r':('s' + i + 'out','t' + i + 'out')},
                                    'r == g*'+str(D)]


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
    
    # Antibody production (from cell to tank)
    
    fnet['trans']['tat'] = {'l0': 0, 'a0': 0} 
    fnet['vhandlers']['vat'] = [{'a':('vat','A'), 'v':('tat','vat')},
                                'a == v']
    tExAb = 't_EX_IgG_f'
    fnet['trans'][tExAb]['l0'] = 0 # IgG production determined by its intensity handler.
    fnet['shandlers']['ha'] = [{'z':('ha',tExAb), 'at':('ha','tat')},
                                'z*'+str(Xmin)+'<= at', 'at <= z*'+str(Xmax)] # , 'z >= 0.00087']

    # Antibody out of the tank (to effluent)
    fnet['trans']['taout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vaout'] = [{'a':('A','vaout'), 'v':('taout','vaout')},
                                  'a == v']
    fnet['shandlers']['saout'] = [{'a':('A','saout'), 'r':('saout','taout')},
                                   'r == a*'+str(D)]    

   
    fnet['obj'] = {'f': "avl['tglcin'] + avl['tglnin'] + avl['tphein'] + avl['targin'] + avl['tasnin'] + avl['taspin'] + avl['tcysin'] + avl['thisin'] + avl['tilein'] + avl['tleuin'] + avl['tlysin'] + avl['tmetin'] + avl['tproin'] + avl['tserin'] + avl['tthrin'] + avl['ttrpin'] + avl['ttyrin'] + avl['tvalin'] + avl['tgluin']", 'sense': 'min'}
    
   
    fnet['extracons'] = ["avl['tat']== " + str(ab_fixed)]
    
    cap_aa = []
    for i in aa:
        # if k != 0:
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
            
def EcoMinimizeMediumim(fnet, D, Xmin, Xmax, ab_fixed, tau): 
    # Adds to fnet the bioreactor dynamics
    #### Parameters, variables and units
    # D (h-1) Dilution rate
    # glcmM (mM) Glucose concentration in medium in mM
    # Xmin (gdcw L-1) Minimum density of cells in tank
    # Xmax (gdcw L-1) Maximum density of cells in tank
    # X (gdcw L-1) Density of cells in the tank 
    # G (mM) Concentration of glucose in the tank 
    # C (mM) Concentration of IgG in the tank 
    # Fluxes are expressed in concentration of the reactant or product per hour
    
    
    aa = ['glc', 'gln', 'phe', 'arg', 'asn', 'asp', 'cys', 'his', 'ile', 'leu', 'lys', 'met', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'glu']
    
    aa_price = [0.92, 2.89, 1.615, 1.615, 1.535, 2.79, 1.615, 1.445, 1.445, 1.615, 1.7, 1.615, 1.615, 1.615, 1.445, 1.615, 0.558, 1.615, 2.21]

    
    # Adding the maximum uptake values
    lim_exch_phe = 0.0065*1.1  # 10% of uncertainty so it is not so restrictive
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
    
    
    ### Bioreactor(tank) concentrations
    fnet['places']['X'] = 1.9512
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
    fnet['places']['A'] = 0.0
    fnet['places']['W'] = 0.0
    fnet['shandlers'] = {}
            
        
    
    for i in aa:
            j = i.capitalize()
            print(j)
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
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_phe)] 
            
            elif i == 'glc':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_glc)] 
            
            elif i == 'gln':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_gln)] 
            
            elif i == 'asn':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asn)] 
            
            elif i == 'ser':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ser)] 
           
            elif i == 'his':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_his)] 
            
            elif i == 'thr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_thr)] 
            
            elif i == 'arg':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_arg)] 
            
            elif i == 'tyr':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) , i + 'g <= ' + str(lim_exch_tyr)] 
            
            elif i == 'val':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_val)] 
            
            elif i == 'met':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_met)] 
            
            elif i == 'trp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_trp)] 
            
            elif i == 'ile':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_ile)] 
            
            elif i == 'leu':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,   i + 'g <= ' + str(lim_exch_leu)] 
            
            elif i == 'lys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_lys)] 
            
            elif i == 'pro':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_pro)] 
            
            elif i == 'asp':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_asp)] 
            
            elif i == 'cys':
                
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax) ,  i + 'g <= ' + str(lim_exch_cys)] 
            
            
            else: 
            
                fnet['shandlers']['h' + i] = [{i + 'g':('h' + i, nut_in), i + 't':('h' + i,'t' + i + 't')}, i + 'g*'+str(Xmin)+'<=' + i + 't', i + 't <=' + i + 'g*'+str(Xmax)] #, 'ug <= 0.1979'] 
                
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
    
    # Antibody production (from cell to tank)
    
    fnet['trans']['tat'] = {'l0': 0, 'a0': 0} 
    fnet['vhandlers']['vat'] = [{'a':('vat','A'), 'v':('tat','vat')},
                                'a == v']
    tExAb = 't_EX_IgG_f'
    fnet['trans'][tExAb]['l0'] = 0 # IgG production determined by its intensity handler.
    fnet['shandlers']['ha'] = [{'z':('ha',tExAb), 'at':('ha','tat')},
                                'z*'+str(Xmin)+'<= at', 'at <= z*'+str(Xmax)] # , 'z >= 0.00087']

    # Antibody out of the tank (to effluent)
    fnet['trans']['taout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vaout'] = [{'a':('A','vaout'), 'v':('taout','vaout')},
                                  'a == v']
    fnet['shandlers']['saout'] = [{'a':('A','saout'), 'r':('saout','taout')},
                                   'r == a*'+str(D)]    

    
    fnet['obj'] = {'f': "0.92*avl['tglcin'] + 2.89*avl['tglnin'] + 1.615*avl['tphein'] + 1.615*avl['targin'] + 1.535*avl['tasnin'] + 2.79*avl['taspin'] + 1.615*avl['tcysin'] + 1.445*avl['thisin'] + 1.445*avl['tilein'] + 1.615*avl['tleuin'] + 1.7*avl['tlysin'] + 1.615*avl['tmetin'] + 1.615*avl['tproin'] + 1.615*avl['tserin'] + 1.445*avl['tthrin'] + 1.615*avl['ttrpin'] + 0.558*avl['ttyrin'] + 1.615*avl['tvalin'] + 2.21*avl['tgluin']", 'sense': 'min'}
    
    
    fnet['extracons'] = ["avl['tat']== " + str(ab_fixed)]

    
    cap_aa = []
    for i in aa:
        # if k != 0:
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

