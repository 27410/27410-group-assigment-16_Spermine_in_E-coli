#%%
from cobra.io import read_sbml_model
from cobra.flux_analysis import production_envelope
import matplotlib as plt 
import math

# --------------- Extracting data from the model and setting up quadratic objective ---------------
model = read_sbml_model("GSM/iAF1260b.xml")

max_growth_rate = model.optimize().objective_value
model.objective = model.reactions.SPRMS
max_spm_rate = model.optimize().objective_value

ratio = max_spm_rate/max_growth_rate
assert ratio > 0

biomass_eq = "BIOMASS_Ec_iAF1260_core_59p81M"

quadratic_objective = model.problem.Objective(
1 * model.reactions.SPRMS.flux_expression + ratio * model.reactions.get_by_id(biomass_eq).flux_expression,
direction='max')

model.objective = quadratic_objective

# --------------- Plotting phenotypic phase plane for 4 different carbon sources depending on O2 uptake ---------------
def plot():
    global model

    import matplotlib.pyplot as plt
    carbon_exchange = {"EX_glc__D_e":"D-Glucose","EX_gal_e":"Galactose","EX_sucr_e":"Sucrose","EX_lcts_e":"Lactose"}

    row = col = int(math.sqrt(len(carbon_exchange)))
    fig, axes = plt.subplots(nrows=row, ncols=col,figsize=(10, 10))
    fig.tight_layout(pad=3.0)

    i = j = 0

    for C_src in carbon_exchange.keys():
        prod_env = production_envelope(
        model, ["EX_o2_e"], objective="SPRMS", carbon_sources=C_src)
        
        prod_env.plot.area(
          x='EX_o2_e', y='carbon_yield_maximum',title=carbon_exchange[C_src],ax=axes[i,j])
        
        j+=1
        if j == col:
            j= 0
            i +=1

# --------------- Get needed medium composition + yield of Spermine in different units ---------------
def get_summary():
    global model
    solution = model.optimize()
    print(model.summary())
    print(f"Max SPRMS flux = {solution.fluxes['SPRMS']} mmol / [gDW h]")
    print(f"= {solution.fluxes['SPRMS']/abs(solution.fluxes['EX_glc__D_e'])} mmol Spermine / mmol of D-Glucose")

if __name__ == "__main__":
    plot()
    get_summary()

# %%
