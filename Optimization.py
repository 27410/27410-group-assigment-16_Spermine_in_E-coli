# %%
from cobra.io import read_sbml_model
import logging 
import time
import pandas
import os 
import prettytable
from Heterologuous_pathway_addition import model_assertion

os.system("clear")

level = logging.INFO	
fmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(level =level, format=fmt)


def main():

    model = read_sbml_model("GSM/iAF1260b.xml")
    model_assertion("GSM/iAF1260b.xml")

    BIOMASS_FVA = pandas.read_csv("FVA_BIOMASS.csv")
    Q_OBJ_FVA = pandas.read_csv("FVA_Q_OBJECTIVE.csv")

    BIOMASS_FVA = {x[0]:(min(x[1],x[2]),max(x[1],x[2])) for x in BIOMASS_FVA.to_numpy()}
    Q_OBJ_FVA = {x[0]:(min(x[1],x[2]),max(x[1],x[2])) for x in Q_OBJ_FVA.to_numpy()}

    reactions = list(BIOMASS_FVA.keys())

    MUST_U = list()
    MUST_L = list()
    MUST_X =list()

    Overlapping = list()

    delta_fluxes = dict()

    for values_BIOMASS,values_Q_OBJECTIVE in zip(BIOMASS_FVA.items(), Q_OBJ_FVA.items()):
        
        assert values_BIOMASS[0] == values_Q_OBJECTIVE[0]
        reaction =  values_BIOMASS[0]

        bounds_BIOMASS = min(values_BIOMASS[1]) , max(values_BIOMASS[1])
        bounds_Q_OBJECTIVE = min(values_Q_OBJECTIVE[1])  , max(values_Q_OBJECTIVE[1])

        if bounds_Q_OBJECTIVE == bounds_BIOMASS:
            reactions.remove(reaction)
            continue
        
        if bounds_Q_OBJECTIVE == (0,0):
            MUST_X.append(reaction)
            continue

        if bounds_Q_OBJECTIVE[1] < bounds_BIOMASS[0]:
            
            delta_fluxes[reaction] = bounds_Q_OBJECTIVE[0] - bounds_BIOMASS[1]
            MUST_L.append(reaction)
            op = "Down-reg"

        elif bounds_Q_OBJECTIVE[0] > bounds_BIOMASS[1]:

            delta_fluxes[reaction]  = bounds_Q_OBJECTIVE[0] - bounds_BIOMASS[1]
            MUST_U.append(reaction)
            op = "Up-reg"
        
        else:
            Overlapping.append(reaction)
            op = f"???"

        #print(f"{reaction} -> {op}, {bounds_BIOMASS=} VS {bounds_Q_OBJECTIVE=}")

    # ------------------------- RoI = reaction of interest ---------------------------
   
    """print(f"{MUST_U=}")
    print(f"{MUST_L=}")
    print(f"{MUST_X=}")"""
    RoI = MUST_U+MUST_L+MUST_X
    # ---------------------------------------------------------------------------------
    RoI.remove("SPRMS")

    max_SPRMS_flux = Q_OBJ_FVA["SPRMS"][1]

    model.objective = "BIOMASS_Ec_iAF1260_core_59p81M"
    Vn0 = model.optimize().fluxes["SPRMS"]
    Vn1 = Vn0 + 1

    Strategy = list()
    while Vn0 < Vn1:
        ranking = dict.fromkeys(RoI,0)
        
        Vn0 = Vn1
        for r in RoI:
            with model:
                model.reactions.get_by_id(r).bounds = Q_OBJ_FVA[r]
                new_V = model.optimize().fluxes["SPRMS"]
                delta_V = new_V - Vn0
                if delta_V > 0:
                    ranking[r] = delta_V
        
        ranking =  dict(sorted(ranking.items(), key=lambda item:item[1],reverse = True))
    
        Best_reaction = list(ranking.keys())[0]
        Vn1 = ranking[Best_reaction]
        RoI.remove(Best_reaction)
        model.reactions.get_by_id(Best_reaction).bounds = Q_OBJ_FVA[Best_reaction]

        Equivalent_reactions = [Best_reaction]
        print(ranking)
        
        for reaction,delta_V in ranking.items():
            if abs(delta_V-Vn1)/Vn1 <= 0.03 and reaction != Best_reaction:
                model.reactions.get_by_id(reaction).bounds = Q_OBJ_FVA[reaction]
                Equivalent_reactions.append(reaction)
                RoI.remove(reaction)
        
        Strategy.append(Equivalent_reactions)


    for i,elet in enumerate(Strategy):
        for i,reaction in enumerate(elet):
            if reaction in MUST_U:
                elet[i] = f"MUST_U {reaction}"
            if reaction in MUST_L:
                elet[i] = f"MUST_L {reaction}"
            if reaction in MUST_X:
                elet[i] = f"MUST_X {reaction}"

        print(f"{i} {elet}")

    print(f"{round(Vn1/max_SPRMS_flux*100)}% of max SPRMS flux")

    raise ValueError




    """

    model.objective = model.reactions.SPRMS
    quadratic_objective = model.problem.Objective(
    1 * model.reactions.SPRMS.flux_expression + 3.7* model.reactions.get_by_id(biomass_eq).flux_expression,
    direction='max')

    model.objective = quadratic_objective
    for r1 in MUST_U:

        with model:
            model.reactions.get_by_id(r1).knock_out()
            fluxes = model.optimize().fluxes
            for r2 in MUST_U:
                if r1 == r2:
                    continue 

                if fluxes[r2] < Q_OBJ_FVA[r2][0] or fluxes[r2] > Q_OBJ_FVA[r2][1]:
                    print(f"{r1} > > {r2}")
    """






    complementation_groups = list()
    for r1,delta1 in delta_fluxes.items():

        if delta1 == 0: continue

        for r2,delta2 in delta_fluxes.items():
            
            if r1 == r2:
                continue
            if abs(delta1-delta2)/abs(delta1) <= 0.01:
                print(f"{delta1=} {delta2=} {abs(delta1-delta2)/delta1=}")
                if sum([1 for x in complementation_groups if r1 in x]) == 1:

                    for group in complementation_groups:
                        if r1 in group and r2 not in group:
                            group.append(r2)
                            break
                elif sum([1 for x in complementation_groups if r2 in x]) == 1:

                    for group in complementation_groups:
                        if r2 in group and r1 not in group:
                            group.append(r1)
                            break
                else:
                    complementation_groups.append([r1,r2])

    print(complementation_groups)

    for i1,r1 in enumerate(reactions):
        r1_bounds_WT = [abs(x) for x in BIOMASS_FVA.get(r1)]
        r1_bounds_MUT = [abs(x) for x in Q_OBJ_FVA.get(r1)]
        
        for i2,r2 in enumerate(reactions):

            if r1 == r2:
                continue 

            if (r1,r2) in MUST_UU+MUST_UL+MUST_LL or (r2,r1) in MUST_UU+MUST_UL+MUST_LL:
                continue
            """
            if r1 in MUST_U and r2 in MUST_U:
                MUST_UU.append((r1,r2))
                continue
            if r1 in MUST_U and r2 in MUST_L:
                MUST_UL.append((r1,r2))
                continue
            if r1 in MUST_L and r2 in MUST_U:
                MUST_UL.append((r2,r1))
                continue
            if r1 in MUST_L and r2 in MUST_L:
                MUST_LL.append((r1,r2))
                continue
            """

            r2_bounds_WT = [abs(x) for x in BIOMASS_FVA.get(r2)]
            r2_bounds_MUT = [abs(x) for x in Q_OBJ_FVA.get(r2)]
      
            if min(r1_bounds_MUT)+min(r2_bounds_MUT) > max(r1_bounds_WT)+max(r2_bounds_WT):        
                print(f"UU {r1_bounds_WT=} {r1_bounds_MUT=}")
                print(f"UU {r2_bounds_WT=} {r2_bounds_MUT=}")
                MUST_UU.append((r1,r2))

            elif max(r1_bounds_MUT)+max(r2_bounds_MUT) < min(r1_bounds_WT)+min(r2_bounds_WT):        
                print(f"LL {r1_bounds_WT=} {r1_bounds_MUT=}")
                print(f"LL {r2_bounds_WT=} {r2_bounds_MUT=}")
                MUST_LL.append((r1,r2))

            else:
                LB_WT = min(r1_bounds_WT) - min(r2_bounds_WT)
                UB_WT = max(r1_bounds_WT) - max(r2_bounds_WT)

                LB_MUT = min(r1_bounds_MUT) - min(r2_bounds_MUT)
                UB_MUT = max(r1_bounds_MUT) - max(r2_bounds_MUT)
                
                if UB_WT < LB_MUT:
                    print(f"UL {r1_bounds_WT=} {r1_bounds_MUT=}")
                    print(f"UL {r2_bounds_WT=} {r2_bounds_MUT=}")
                    MUST_UL.append((r1,r2))
                   
                elif UB_MUT < LB_WT:
                    print(f"LU {r1_bounds_WT=} {r1_bounds_MUT=}")
                    print(f"LU {r2_bounds_WT=} {r2_bounds_MUT=}")
                    MUST_UL.append((r2,r1))
                   
    table = prettytable.PrettyTable()
    table.field_names = ["MUST UU","MUST UL","MUST LL"]
    
    _MUST_UU = MUST_UU + [(" "," ") for _ in range(max([len(MUST_UU),len(MUST_UL),len(MUST_LL)])-len(MUST_UU))]
    _MUST_UL = MUST_UL + [(" "," ") for _ in range(max([len(MUST_UU),len(MUST_UL),len(MUST_LL)])- len(MUST_UL))]
    _MUST_LL = MUST_LL + [(" "," ") for _ in range(max([len(MUST_UU),len(MUST_UL),len(MUST_LL)])-len(MUST_LL))]

    for elet in zip(_MUST_UU,_MUST_UL,_MUST_LL):
        table.add_row([f"{elet[0][0]},{elet[0][1]}",f"{elet[1][0]},{elet[1][1]}",f"{elet[2][0]},{elet[2][1]}"])

    print(table)

    assert len([x for x in MUST_UL if x in MUST_UU + MUST_LL]) == 0
    assert len([x for x in MUST_UU if x in MUST_LL + MUST_UL ]) == 0 
    assert len([x for x in MUST_LL if x in MUST_UU + MUST_UL ]) == 0 

    import networkx as nx
    import matplotlib.pyplot as plt

    G = nx.Graph()

    ind_r1 = [x[0] for x in MUST_UU] + [x[1] for x in MUST_UU]
    ind_r2 = [x[0] for x in MUST_UL] + [x[1] for x in MUST_UL]
    ind_r3 = [x[0] for x in MUST_LL] + [x[1] for x in MUST_LL]

    reactions_to_modify = dict.fromkeys(ind_r1+ind_r2+ind_r3)
    print([x for x in reactions_to_modify if x not in MUST_U+MUST_L ])

    color_map = list()
   
    for a,b in MUST_UU+MUST_UL+MUST_LL:
        
        if (a,b) == (0,0):
            continue

        if a not in G.nodes():
            G.add_node(a)
            if a in MUST_U:
               color_map.append("green")
            elif a in MUST_L:
               color_map.append("red")
            else:
                color_map.append("yellow")
        
        if b not in G.nodes():
            G.add_node(b)
            if b in MUST_U:
                color_map.append("green")
            elif b in MUST_L:
               color_map.append("red")
            else:
                color_map.append("yellow")

        G.add_edge(a,b)

    nx.draw(G, with_labels=True, font_weight='bold',node_color = color_map)        
    
    MUST_sets =  [x[0] for x in MUST_UU+MUST_UL+MUST_LL] + [x[1] for x in MUST_UU+MUST_UL+MUST_LL]
    Occurences  = dict.fromkeys(MUST_sets,0)

    for react in MUST_sets:
        Occurences[react] +=  1
    
    print(Occurences)
    plt.show()


if __name__ == "__main__":
        main()


# %%
