# %%
from cobra.io import read_sbml_model
import logging 
import time
import random
import os 
from Heterologuous_pathway_addition import model_assertion
from cobra.util import get_solver_name

def main():

    level = logging.INFO	
    fmt = '[%(levelname)s] %(asctime)s - %(message)s'
    logging.basicConfig(level =level, format=fmt)

    model_assertion("GSM/iAF1260b.xml")
    model = read_sbml_model("GSM/iAF1260b.xml")

    print(f"Solver used: {get_solver_name()}")

    biomass_eq = "BIOMASS_Ec_iAF1260_core_59p81M"

    # ------------- Determining weighting ratio (alpha in the report ) -------------
    max_growth_rate = model.optimize().objective_value
    model.objective = model.reactions.SPRMS
    max_spm_rate = model.optimize().objective_value

    ratio = max_spm_rate/max_growth_rate
    assert ratio > 0

    #------------- Setting the quadratic objective -------------
    quadratic_objective = model.problem.Objective(
    1 * model.reactions.SPRMS.flux_expression + ratio * model.reactions.get_by_id(biomass_eq).flux_expression,
    direction='max')

    from cobra.flux_analysis import flux_variability_analysis

    # ------------- FVA of Wild Type model (WT in the report) -------------
    #                     /!\ Need high computing power

    model.objective = "BIOMASS_Ec_iAF1260_core_59p81M"
    start = time.perf_counter()
    FVA1 = flux_variability_analysis(model,model.reactions,loopless=True)
    FVA1.to_csv(f"FVA_BIOMASS.csv")

    print(model.objective)
    stop = time.perf_counter()
    print(f"FV1 Done in {stop-start}s")
    
    # ------------- FVA of Optimal strain (OS in the report) -------------
    #                     /!\ Need high computing power

    model.objective = quadratic_objective
    start = time.perf_counter()
    FVA2 = flux_variability_analysis(model,model.reactions,loopless=True)
    FVA2.to_csv(f"FVA_Q_OBJECTIVE.csv")
    print(model.objective)

    stop = time.perf_counter()
    print(f"FVA2 Done in {stop-start}s")

if __name__ == "__main__":
        main()

