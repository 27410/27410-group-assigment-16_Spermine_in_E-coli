# %%
from cobra.io import read_sbml_model
import logging 
import time
import random
import os 

def main():
    level = logging.INFO	
    fmt = '[%(levelname)s] %(asctime)s - %(message)s'
    logging.basicConfig(level =level, format=fmt)

    from cobra import Reaction, Metabolite

    model = read_sbml_model("GSM/iAF1260b.xml")

    #-------------- Adding spm_c and spm_e to model's metabolites --------------

    Spermine_c = Metabolite(id="spm_c",name="spermine",compartment="c")
    Spermine_e = Metabolite(id="spm_e",name="spermine",compartment="e")

    model.add_metabolites([Spermine_c,Spermine_e])

    #-------------- Adding EX_sprm_e to model's reactions --------------
    EX_sprm_e = Reaction("EX_sprm_e")
    EX_sprm_e.add_metabolites({
        model.metabolites.spm_e:-1,
    })
    model.add_reaction(EX_sprm_e)
    logging.debug(f"EX_sprm_e added")

    #-------------- Adding SPRMS reaction --------------

    SPRMS = Reaction("SPRMS")
    SPRMS.add_metabolites({
        model.metabolites.ametam_c:-1,
        model.metabolites.spmd_c: -1,

        model.metabolites.get_by_id("5mta_c"): 1,
        model.metabolites.h_c:1,
        model.metabolites.spm_c:1
    })
    model.add_reaction(SPRMS)

    model.add_boundary(model.metabolites.spm_c, type = "demand")

    biomass_eq = "BIOMASS_Ec_iAF1260_core_59p81M"


    max_growth_rate = model.optimize().objective_value
    model.objective = model.reactions.SPRMS
    max_spm_rate = model.optimize().objective_value

    ratio = max_spm_rate/max_growth_rate

    quadratic_objective = model.problem.Objective(
    1 * model.reactions.SPRMS.flux_expression + ratio * model.reactions.get_by_id(biomass_eq).flux_expression,
    direction='max')

    from cobra.flux_analysis import flux_variability_analysis

    index = random.randint(0,len(model.reactions)) - 100
    a = [model.reactions.SPRMS,model.reactions.GSPMDS,model.reactions.GTHS] + model.reactions[index:index+20]

    model.objective = "BIOMASS_Ec_iAF1260_core_59p81M"
    start = time.perf_counter()
    FVA1 = flux_variability_analysis(model,a,loopless=True)
    FVA1.to_csv(f"FVA_BIOMASS.csv")

    print(model.objective)
    stop = time.perf_counter()
    print(f"FV1 Done in {stop-start}s")
    
    model.objective = quadratic_objective
    start = time.perf_counter()
    FVA2 = flux_variability_analysis(model,a,loopless=True)
    FVA2.to_csv(f"FVA_Q_OBJECTIVE.csv")
    print(model.objective)

    stop = time.perf_counter()
    print(f"FVA2 Done in {stop-start}s")

if __name__ == "__main__":
        main()

