# %%
from cobra.io import read_sbml_model, save_json_model
import logging 
import time
import pandas
from Heterologuous_pathway_addition import model_assertion

level = logging.INFO	
fmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(level =level, format=fmt)


def main():

    model = read_sbml_model("GSM/iAF1260b.xml")
    save_json_model(model,f"GSM/iAF1260b_WT.json")

    model_assertion("GSM/iAF1260b.xml")

    #------------- Extract FVA data from output files -------------
    BIOMASS_FVA = pandas.read_csv("FVA_BIOMASS.csv")
    Q_OBJ_FVA = pandas.read_csv("FVA_Q_OBJECTIVE.csv")

    #------------- Format data in a more convenient way -------------
    BIOMASS_FVA = {x[0]:(min(x[1],x[2]),max(x[1],x[2])) for x in BIOMASS_FVA.to_numpy()}
    Q_OBJ_FVA = {x[0]:(min(x[1],x[2]),max(x[1],x[2])) for x in Q_OBJ_FVA.to_numpy()}

    reactions = list(BIOMASS_FVA.keys())

    MUST_U = list()
    MUST_L = list()
    MUST_X =list()

    Overlapping = list()

    delta_fluxes = dict()

    #------------- Establishing MUST sets -------------

    for values_BIOMASS,values_Q_OBJECTIVE in zip(BIOMASS_FVA.items(), Q_OBJ_FVA.items()):
        
        assert values_BIOMASS[0] == values_Q_OBJECTIVE[0] # assess reaction being analyzed is the same in WT and OS 
        reaction =  values_BIOMASS[0]
        
        bounds_BIOMASS = min(values_BIOMASS[1]) , max(values_BIOMASS[1])
        bounds_Q_OBJECTIVE = min(values_Q_OBJECTIVE[1])  , max(values_Q_OBJECTIVE[1])

        if bounds_Q_OBJECTIVE == bounds_BIOMASS: # No changes
            reactions.remove(reaction)
            continue
        
        if bounds_Q_OBJECTIVE == (0,0): # Knock-outs -> MUST_X
            MUST_X.append(reaction)

        elif bounds_Q_OBJECTIVE[1] < bounds_BIOMASS[0]: # Up-regulation -> MUST_U
            
            delta_fluxes[reaction] = bounds_Q_OBJECTIVE[0] - bounds_BIOMASS[1]
            MUST_L.append(reaction)

        elif bounds_Q_OBJECTIVE[0] > bounds_BIOMASS[1]: # Down-regulation -> MUST_L

            delta_fluxes[reaction]  = bounds_Q_OBJECTIVE[0] - bounds_BIOMASS[1]
            MUST_U.append(reaction)
        
        else:   
            Overlapping.append(reaction)

    # ------------------------- RoI = reaction of interest ---------------------------
   
    RoI = MUST_U + MUST_L + MUST_X
    RoI.remove("BIOMASS_Ec_iAF1260_core_59p81M") # Removing Biomass equation because doesn't represent 1 reaction 

    max_SPRMS_flux = Q_OBJ_FVA["SPRMS"][1] 

    model.objective = "BIOMASS_Ec_iAF1260_core_59p81M"
    Vn0 = model.optimize().fluxes["SPRMS"] # V_SPRMS(n) with n being the number of iterations of the loop
    Vn1 = Vn0 + 1   # V_SPRMS(n+1) 
    Strategy = list()
    
    while Vn0 < Vn1: # True as long as SPRMS flux keeps increasing 
        ranking = dict()
        Vn0 = model.optimize().fluxes["SPRMS"]

        for r in RoI:
            with model: # Temporary model 
                model.reactions.get_by_id(r).bounds = Q_OBJ_FVA[r] # Reaction bounds set as in OS 
                new_V = model.optimize().fluxes["SPRMS"] # SPRMS flux with modified Reaction 
                if new_V > Vn0: # If SPRMS flux increased, add reaction to ranking list 
                    ranking[r] = new_V
                    print(f"{new_V-Vn0=}")
        
        ranking =  dict(sorted(ranking.items(), key=lambda item:item[1],reverse=True)) # Ranking the reaction modifications by their associated increase in SPRMS flux

        try:
            Best_reaction = list(ranking.keys())[0] 
        except IndexError:
            break # If no element if ranking = no modifications increased the SPRMS flux 

        model.reactions.get_by_id(Best_reaction).bounds = Q_OBJ_FVA[Best_reaction] # adding Best reaction OS bounds as a new constraint
        Vn1 = model.optimize().fluxes["SPRMS"]
        RoI.remove(Best_reaction)

        Equivalent_reactions = [Best_reaction]
        
        # ------------------ Searching for the Reaction in a 3% range of the best increase ------------------
        for reaction,Vn in ranking.items():
            if abs(Vn-Vn1)/Vn1 <= 0.03 and reaction not in Equivalent_reactions:
                model.reactions.get_by_id(reaction).bounds = Q_OBJ_FVA[reaction]  # adding equivalent reaction OS bounds as a new constraint
                Equivalent_reactions.append(reaction)
                RoI.remove(reaction)

        Strategy.append(Equivalent_reactions)
        print(f"{Vn1-Vn0=}")
        

    # ------------------ Associating reaction modifications with their nature: Up-reg / Down-reg / Knock-out ------------------
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
        
        output_list.append(elet[:10]) # [:10] is only for display purposes and can be deleted for full results (long array) 
        print(f"{i} {elet[:10]} ")

    print(f"{new_V/max_SPRMS_flux*100} % of max SPRMS flux") # How close to optimal flux does the Strategy place the new strain model

    # ------------------ Printing results to output file "Optimization startegy.txt" ------------------

    with open("Optimization startegy.txt","w") as f:
        f.write(f"{new_V/max_SPRMS_flux*100} % of max SPRMS flux for any reaction from rank 1 + any reaction from rank 2\n")
        for i,reaction_list in enumerate(output_list):
            f.write(f"\n#-------------------------------- Rank {i+1} reactions --------------------------------\n")
            f.writelines("\n".join(reaction_list))
        
        f.write("\n...")

    save_json_model(model,f"GSM/iAF1260b_OS.json") # Save new strain model

if __name__ == "__main__":
        main()


# %%
