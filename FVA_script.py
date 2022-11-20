# %%
from cobra.io import read_sbml_model
import logging 
import time
import random
import os 
from Heterologuous_pathway_addition import model_assertion


def main():
    level = logging.INFO	
    fmt = '[%(levelname)s] %(asctime)s - %(message)s'
    logging.basicConfig(level =level, format=fmt)

    model_assertion("GSM/iAF1260b.xml")
    model = read_sbml_model("GSM/iAF1260b.xml")

    biomass_eq = "BIOMASS_Ec_iAF1260_core_59p81M"

    max_growth_rate = model.optimize().objective_value
    model.objective = model.reactions.SPRMS
    max_spm_rate = model.optimize().objective_value

    ratio = max_spm_rate/max_growth_rate
    assert ratio > 0

    quadratic_objective = model.problem.Objective(
    1 * model.reactions.SPRMS.flux_expression + ratio * model.reactions.get_by_id(biomass_eq).flux_expression,
    direction='max')

    from cobra.flux_analysis import flux_variability_analysis

    index = random.randint(0,len(model.reactions)) - 100
    a = [model.reactions.SPRMS,model.reactions.SPMS,model.reactions.GSPMDS,model.reactions.GTHS,model.reactions.ORNDC] + model.reactions[index:index+50]

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

