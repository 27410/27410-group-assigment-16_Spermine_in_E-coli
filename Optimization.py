# %%
from cobra.io import read_sbml_model
import logging 
import time
import pandas
import os 
import prettytable
from Heterologuous_pathway_addition import model_assertion
import itertools

level = logging.INFO	
fmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(level =level, format=fmt)


def main():

    model = read_sbml_model("GSM/iAF1260b.xml")

    model_assertion("GSM/iAF1260b.xml")

    BIOMASS_FVA = pandas.read_csv("full_FVA_BIOMASS.csv")
    Q_OBJ_FVA = pandas.read_csv("full_FVA_Q_OBJECTIVE.csv")

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

        elif bounds_Q_OBJECTIVE[1] < bounds_BIOMASS[0]:
            
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
    RoI = MUST_U + MUST_L + MUST_X
    # ---------------------------------------------------------------------------------
    RoI.remove("BIOMASS_Ec_iAF1260_core_59p81M")
    """
    a = ['UNK3', 'DKMPPD3', 'SPRMS', 'MTAN', 'MTRI', 'SPMS', 'MTRK', 'DM_spm_c', 'MDRPD', 'ADMDC']
    RoI = [x for x in RoI if x not in a]
    """
    max_SPRMS_flux = Q_OBJ_FVA["SPRMS"][1]

    model.objective = "BIOMASS_Ec_iAF1260_core_59p81M"
    Vn0 = model.optimize().fluxes["SPRMS"]
    Vn1 = Vn0 + 1
    Strategy = list()
    
    while Vn0 < Vn1:
        ranking = dict()
        Vn0 = model.optimize().fluxes["SPRMS"]

        for r in RoI:
            with model:
                model.reactions.get_by_id(r).bounds = Q_OBJ_FVA[r]
                new_V = model.optimize().fluxes["SPRMS"]
                if new_V > Vn0:
                    ranking[r] = new_V
                    print(f"{new_V-Vn0=}")
        
        ranking =  dict(sorted(ranking.items(), key=lambda item:item[1],reverse=True))
        print(ranking)
        try:
            Best_reaction = list(ranking.keys())[0]
        except IndexError:
            break

        model.reactions.get_by_id(Best_reaction).bounds = Q_OBJ_FVA[Best_reaction]
        Vn1 = model.optimize().fluxes["SPRMS"]
        RoI.remove(Best_reaction)

        Equivalent_reactions = [Best_reaction]
        
        for reaction,Vn in ranking.items():
            if abs(Vn-Vn1)/Vn1 <= 0.001 and reaction not in Equivalent_reactions:
                model.reactions.get_by_id(reaction).bounds = Q_OBJ_FVA[reaction]
                #print(f"{Vn-Vn1=}")
                Equivalent_reactions.append(reaction)
                RoI.remove(reaction)

        Strategy.append(Equivalent_reactions)
        print(f"{Vn1-Vn0=}")
        

    #os.system("clear")

    output_list = list()
    for i,elet in enumerate(Strategy):
        for i2,reaction in enumerate(elet):
            if reaction in MUST_U:
                elet[i2] = f"Up-reg {model.reactions.get_by_id(reaction).name} ({model.reactions.get_by_id(reaction).id})\n{model.reactions.get_by_id(reaction).reaction}\n"
            elif reaction in MUST_L:
                elet[i2] = f"Down-reg {model.reactions.get_by_id(reaction).name} ({model.reactions.get_by_id(reaction).id})\n{model.reactions.get_by_id(reaction).reaction}\n"
            elif reaction in MUST_X:
                elet[i2] = f"Knock-out {model.reactions.get_by_id(reaction).name} ({model.reactions.get_by_id(reaction).id})\n{model.reactions.get_by_id(reaction).reaction}\n"
            else:
                raise ValueError(f"{reaction} not part of any MUST sets")
        
        output_list.append(elet[:10])
        print(f"{i} {elet[:10]} ")

    print(f"{new_V/max_SPRMS_flux*100} % of max SPRMS flux")

    with open("Optimization startegy.txt","w") as f:
        f.write(f"{new_V/max_SPRMS_flux*100} % of max SPRMS flux for any reaction from rank 1 + any reaction from rank 2\n")
        for i,reaction_list in enumerate(output_list):
            f.write(f"\n#-------------------------------- Rank {i+1} reactions --------------------------------\n")
            f.writelines("\n".join(reaction_list))
        
        f.write("\n...")

if __name__ == "__main__":
        main()


# %%
