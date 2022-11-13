# %%
from cobra.io import read_sbml_model
import logging 
import time
import pandas
import os 
import prettytable

os.system("clear")

level = logging.INFO	
fmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(level =level, format=fmt)


def main():

    BIOMASS_FVA = pandas.read_csv("FVA_BIOMASS.csv")
    Q_OBJ_FVA = pandas.read_csv("FVA_Q_OBJECTIVE.csv")


    BIOMASS_FVA = {x[0]:(x[1],x[2]) for x in BIOMASS_FVA.to_numpy()}
    Q_OBJ_FVA = {x[0]:(x[1],x[2]) for x in Q_OBJ_FVA.to_numpy()}

    reactions = list(BIOMASS_FVA.keys())

    MUST_U = list()
    MUST_L = list()
    Overlapping = list()

    for values_BIOMASS,values_Q_OBJECTIVE in zip(BIOMASS_FVA.items(), Q_OBJ_FVA.items()):
        
        assert values_BIOMASS[0] == values_Q_OBJECTIVE[0]
        reaction =  values_BIOMASS[0]

        bounds_BIOMASS = values_BIOMASS[1][0] , values_BIOMASS[1][1]
        bounds_Q_OBJECTIVE = values_Q_OBJECTIVE[1][0]  , values_Q_OBJECTIVE[1][1]

        if bounds_Q_OBJECTIVE == bounds_BIOMASS:
            reactions.remove(reaction)
            continue

        elif bounds_Q_OBJECTIVE[1] < bounds_BIOMASS[0]:
            MUST_L.append(reaction)
            op = "Down-reg"

        elif bounds_Q_OBJECTIVE[0] > bounds_BIOMASS[1]:
            MUST_U.append(reaction)
            op = "Up-reg"

        else:
            Overlapping.append(reaction)
            op = f"???"

        logging.debug(f"{reaction} -> {op}, {bounds_BIOMASS=} VS {bounds_Q_OBJECTIVE=}")

    print(MUST_U)
    print(MUST_L)
    MUST_UU = list()
    MUST_UL = list()
    MUST_LL = list()

    for i1,r1 in enumerate(reactions):
        r1_bounds_WT = BIOMASS_FVA.get(r1)
        r1_bounds_MUT = Q_OBJ_FVA.get(r1)
        
        for i2,r2 in enumerate(reactions):

            if r1 == r2:
                continue 

            if (r1,r2) in MUST_UU+MUST_UL+MUST_LL or (r2,r1) in MUST_UU+MUST_UL+MUST_LL:
                continue

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

            logging.debug(f"{r1}/{r2} ???")
            r2_bounds_WT = BIOMASS_FVA.get(r2)
            r2_bounds_MUT = Q_OBJ_FVA.get(r2)
      
            if min(r1_bounds_MUT)+min(r2_bounds_MUT) > max(r1_bounds_WT)+max(r2_bounds_WT):        
                MUST_UU.append((r1,r2))

            elif max(r1_bounds_MUT)+max(r2_bounds_MUT) < min(r1_bounds_WT)+min(r2_bounds_WT):        
                MUST_LL.append((r1,r2))

            else:
                LB_WT = min(r1_bounds_WT) - min(r2_bounds_WT)
                UB_WT = max(r1_bounds_WT) - max(r2_bounds_WT)

                LB_MUT = min(r1_bounds_MUT) - min(r2_bounds_MUT)
                UB_MUT = max(r1_bounds_MUT) - max(r2_bounds_MUT)

                if UB_WT < LB_MUT:
                    MUST_UL.append((r1,r2))
                   
                elif UB_MUT < LB_WT:
                    MUST_UL.append((r2,r1))
                   
    table = prettytable.PrettyTable()
    table.field_names = ["MUST UU","MUST UL","MUST LL"]
    
    _MUST_UU = MUST_UU + [(" "," ") for _ in range(max([len(MUST_UU),len(MUST_UL),len(MUST_LL)])-len(MUST_UU))]
    _MUST_UL = MUST_UL + [(" "," ") for _ in range(max([len(MUST_UU),len(MUST_UL),len(MUST_LL)])- len(MUST_UL))]
    _MUST_LL = MUST_LL + [(" "," ") for _ in range(max([len(MUST_UU),len(MUST_UL),len(MUST_LL)])-len(MUST_LL))]

    for elet in zip(_MUST_UU,_MUST_UL,_MUST_LL):
        print("jkbh")
        table.add_row([f"{elet[0][0]},{elet[0][1]}",f"{elet[1][0]},{elet[1][1]}",f"{elet[2][0]},{elet[2][1]}"])

    print(table)

    assert len([x for x in MUST_UL if x in MUST_UU ]) == 0
    assert len([x for x in MUST_UL if x in MUST_LL ]) == 0
    assert len([x for x in MUST_UU if x in MUST_LL ]) == 0 

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
    plt.show()


if __name__ == "__main__":
        main()


# %%
