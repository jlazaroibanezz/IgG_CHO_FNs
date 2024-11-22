# -*- coding: utf-8 -*-

# Functions to compute the fluxes of a
# flexible Net integrating a CHO cells' metabolic network,
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
from nutrients_HP_git import genCHOBiorFNmax, genCHOBiorFNmin, MinimizeMediumim, EcoMinimizeMediumim

# Constants
glucose_mmass = 180.1577 # Glucose molar mass (g/mol)
tau = 0.0022 # toxicity (not included in this work)



def comProductivity(fnet, abmw, c, glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu, D = 0.025, X0 = 0.5, Xf = 5.0, XNsamps = 4, Xint = 0.0002):
    
    '''
    Step 1. Computes the antibody maximal production for a given dilution rate, biomass and the concentrations of the metabolites in the medium
    Step 2. Stores the antibody maximal production
    Step 3. Fixes the antibody maximal production and maximizes and minimizes the uptake of each aminoacid. This last step involves the two functions introduced in nutrients.py
    '''
    
    bio_list = []
    anti_list = []
    glc_list = []
    gln_list = []
    phe_list = []
    solution_list = [] 
    limiting_met = []
    thrmax = []
    thrmin = []
    lysmax = []
    lysmin = []

    
    nut = []
    aa = ['glc', 'gln', 'phe', 'arg', 'asn', 'asp', 'cys', 'his', 'ile', 'leu', 'lys', 'met', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'glu']
    for i in aa:
        j = i.capitalize()
        nut.append(j)
    
    ori_sol_list_index = -1

    
    for X in np.linspace(X0, Xf, XNsamps):   # Sampling the X (biomass)
        fnet = loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = "cplex_direct")  # load the model
        genCHOBiorFN(fnet, D, glcmM, gln, phe, arg, asn, asp, cys, his, ile, leu, lys, met, pro, ser, thr, trp, tyr, val, glu, Xmin = X, Xmax = X)  # Generate FN from model
        netobj = FNFactory(fnet) # Build net object
        ori_sol_list = []  # list containing the antibody solutions obtained in Step 1
        aminolistmax = []  # list containing the maximum uptake for each aminoacid
        aminolistmin = []  # list containing the minimum uptake of each aminoacid
        dict_solutions = {}  # dictionary containg the values in aminolistmax, aminolistmin and its difference
	
	# Start simulation 
	
        try:
                
                netobj.optimize()
                solution = netobj.objval
                biomass =  netobj.places['X'].avm
                glc_tank = netobj.places['Glc'].avm
                glc_out = netobj.trans['tglcout'].avl
                glc_in = netobj.trans['tglcin'].avl
                glc_incell = netobj.trans['tglct'].avl
                print("  Growth with X in [", f"{X:.2f}",",", f"{X+Xint:.2f}", "] gdcw L-1:")
              
                for k in nut:
                    print(k + " in tank: ", netobj.places[k].avm, 'mM')
                
                print(netobj.trans['tat'].avl)
                print("Solution:", solution, 'mM h-1') 
                conv_sol = solution*24*abmw*10**9*c/biomass  # conversion units
                print("converted solution", conv_sol)
                print(netobj.trans['t_herceptin_subunit_f'].avl)
                print('Glucose in tank: ', glc_tank, 'mM') 
                print('Glucose leaving through the effluent:', glc_out, 'mM h-1') # glucose in tank * Dilution rate, e.g =  12.168086498662001 * 0,026
                print('Glucose entering the tank:', glc_in, 'mM h-1')  # glucose in supply medium * Dilution rate, e.g =  35,076 * 0,026
                print('Glucose entering the cell:', glc_incell, 'mM h-1') # biomass * uptake rate (here is a constraint), e.g =  1,43 * 0,418
                print('Antibody concentration: ', netobj.trans['taout'].avl/D, 'mM') 
                print('Antibody concentration 2: ', netobj.places['A'].avm, 'mM') 
                print('Antibody concentration in g/l: ', netobj.places['A'].avm*abmw/1000, 'g/l')
                print('Antibody leaving through effluent: ', netobj.trans['taout'].avl, 'mM h-1')
                print('Glutamine in tank: ', netobj.places['Gln'].avm, 'mM')
                print('Glutamine leaving through the effluent:', netobj.trans['tglnout'].avl, 'mM h-1')
                print('Glutamine entering the tank:', netobj.trans['tglnin'].avl, 'mM h-1')
                print('Glutamine entering the cell:', netobj.trans['tglnt'].avl, 'mM h-1')
                print('Phenylalanine in tank: ', netobj.places['Phe'].avm, 'mM')
                print('Phenylalanine leaving through the effluent:', netobj.trans['tpheout'].avl, 'mM h-1')
                print('Phenylalanine entering the tank:', netobj.trans['tphein'].avl/D, 'mM h-1')
                print('Phenylalanine entering the cell:', netobj.trans['tphet'].avl, 'mM h-1')
                print('Antibody production intracellular flux: ', netobj.trans['t_mAb2_f'].avl, 'mmol_per_gDW_per_hr')
                print('Biomass or growth production intracellular flux (equal to dilution rate): ', netobj.trans['t_BIOMASS_cho_f'].avl, 'h-1')
                print('growth production x cell number: ', netobj.trans['txout'].avl, 'h-1')
                print('Glucose exchange flux: ', netobj.trans['t_EX_glc__D_e_b'].avl, 'mmol_per_gDW_per_hr')
                print('Glutamine exchange flux: ', netobj.trans['t_EX_gln__L_e_b'].avl, 'mmol_per_gDW_per_hr')
                print('Phenylalanine exchange flux: ', netobj.trans['t_EX_phe__L_e_b'].avl, 'mmol_per_gDW_per_hr')
                
                print('Biomass: ', biomass)
                bio_list.append(biomass)
                ori_sol_list.append(solution)
                anti_list.append(conv_sol)
                glc_list.append(glc_tank)
                gln_list.append(netobj.places['Gln'].avm)
                phe_list.append(netobj.places['Phe'].avm)
                solution_list.append(ori_sol_list)
                
                # Medium minimization
                
                try: 
                
			
                    fnet = loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = "cplex_direct")  # change to glpk or cplex
                    
                    
                    EcoMinimizeMediumim(fnet, D=D, Xmin = X, Xmax = X, ab_fixed = ori_sol_list[0], tau = tau)  # Function defined in nutrients_HP_git.py
                    
                    
                    netobj = FNFactory(fnet)

                    netobj.optimize()
                    
                    ecominmediumim = netobj.objval
                    for i in aa:
                    	j = i.capitalize()
                    	print(j + " enters in tank: ", netobj.trans['t' + str(i) + 'in'].avl/D, 'mM')
                    	
                    print(ecominmediumim)
                    
                    fnet = loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = "cplex_direct")  # change to glpk or cplex
                    MinimizeMediumim(fnet, D=D, Xmin = X, Xmax = X, ab_fixed = ori_sol_list[0], tau = tau)
                    netobj = FNFactory(fnet)
                    netobj.optimize()
                    minmediumim = netobj.objval
                    for i in aa:
                    	j = i.capitalize()
                    	print(j + " enters in tank: ", netobj.trans['t' + str(i) + 'in'].avl/D, 'mM')
                    	
                    print('This is the minimum sum of concentrations in the medium --- ' , minmediumim , ' mM')
                    
                    
                    # Finding the limiting metabolite
                    
                    for a in aa:        

                        
                            solutionmax = genCHOBiorFNmax(fnet, D=D, glcmM = 35.15619796, gln = 7.872825297, phe = 1.275823895, arg = 1.946111508, asn = 5.122945878, asp = 1.354860924, cys = 0.379, his = 1.216191003, ile = 2.944437965, leu = 3.964891601, lys = 2.380196919, met = 1.010750472, pro = 4.621036929, ser = 4.956440488, thr = 2.774941185, trp = 0.926515767, tyr = 0.978250505, val = 2.946237191, glu = 1.84340765248991, Xmin = X, Xmax = X, ab_fixed = ori_sol_list[0], aminosol = a) # importante poner la X en Xmax y Xmin
                            
                            aminolistmax.append(solutionmax)
                                    
                            solutionmin = genCHOBiorFNmin(fnet, D=D, glcmM = 35.15619796, gln = 7.872825297, phe = 1.275823895, arg = 1.946111508, asn = 5.122945878, asp = 1.354860924, cys = 0.379, his = 1.216191003, ile = 2.944437965, leu = 3.964891601, lys = 2.380196919, met = 1.010750472, pro = 4.621036929, ser = 4.956440488, thr = 2.774941185, trp = 0.926515767, tyr = 0.978250505, val = 2.946237191, glu = 1.84340765248991, Xmin = X, Xmax = X, ab_fixed = ori_sol_list[0], aminosol = a) # importante poner la X en Xmax y Xmin
                                
                            aminolistmin.append(solutionmin)
                            
                    for a, x, y in zip(aa, aminolistmax, aminolistmin):
            
                        dict_solutions[a] = x, y, x-y
        
                    print(dict_solutions)

                    dic_conc = {}
                    dif_conc = []

                    for key in dict_solutions:
                        dif_conc.append(dict_solutions[key][2])

                    for a, x in zip(aa, dif_conc):
                            dic_conc[a] = x
                            
                    print('The limiting metabolite is: ', min(dic_conc, key = dic_conc.get))  # gets the minimum value in dictionary that contains only the differences between maximum and minimum uptakes
                    
                    limiting_met.append(min(dic_conc, key = dic_conc.get)) # adds the metabolite whose difference between maximum and minimum uptake is the minimum for all metabolites to the list limiting_met. For each X and D
                    
                    for key in dict_solutions:
                        if key == 'thr':
                            thrmax.append(dict_solutions[key][0])
                            thrmin.append(dict_solutions[key][1])
                        elif key == 'lys':
                            lysmax.append(dict_solutions[key][0])
                            lysmin.append(dict_solutions[key][1])
                    
                    print(lysmax, lysmin, thrmax, thrmin)
                
                    
                
                except:
                    
                    pass
                    
                        
            
                 
        except:
                print("  Problem NOT feasible") # with X in [", f"{X:.2f}",",", f"{X+Xint:.2f}", "] gdcw L-1")
                bio_list.append(X)
                anti_list.append('Inf')
                glc_list.append('Inf')
                gln_list.append('Inf')
                phe_list.append('Inf')
                ori_sol_list.append('Inf')
                solution_list.append(ori_sol_list)
                limiting_met.append('No metabolite')
                thrmax.append('no data')
                thrmin.append('no data')
                lysmax.append('no data')
                lysmin.append('no data')
        

        
    return bio_list, anti_list, glc_list, gln_list, phe_list, solution_list, limiting_met, thrmax, thrmin, lysmax, lysmin #, dead_reac


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

    
       

    # Defining the objective
    
    fnet['obj'] = {'f': "avl['tat']", 'sense': 'max'}
    fnet['extrans'] = 'all'
    
    
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


fnet = loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = "cplex")

comProductivity(fnet, abmw = 145881, c = 300e-12, glcmM = 35.15619796, gln = 7.872825297, phe = 1.275823895, arg = 1.946111508, asn = 5.122945878, asp = 1.354860924, cys = 0.379, his = 1.216191003, ile = 2.944437965, leu = 3.964891601, lys = 2.380196919, met = 1.010750472, pro = 4.621036929, ser = 4.956440488, thr = 2.774941185, trp = 0.926515767, tyr = 0.978250505, val = 2.946237191, glu = 1.84340765248991, D = 0.0166, X0 = 3.18727272727273, Xf = 3.18727272727273, Xint = 0, XNsamps = 1) 


