#%%
from cobra.io import read_sbml_model
import logging 
import time

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

model.objective = quadratic_objective
solution = model.optimize(objective_sense=None)

from cobra.flux_analysis import flux_variability_analysis
start = time.perf_counter()
flux_variability_analysis(model,[SPRMS,model.reactions.SPMS,model.reactions.GSPMDS],loopless=True)
stop = time.perf_counter()
print(f"Done in {stop-start}s")
