from cobra.io import read_sbml_model
import glob 

import requests
import logging 
import csv 

level = logging.INFO	
fmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(level =level, format=fmt)

# ------------------- Extracting name of the models containing spmd metabolite from "list_spmd_models.txt" -------------------
spmd_models = list()
with open("list_spmd_models.txt","r") as f:
    for line in f.readlines():
        if line[0] == "i":
            spmd_models.append(line.split("\n")[0])
        else:
            logging.debug(f"erased {line}")

spmd_modelsf = list(dict.fromkeys(spmd_models))
logging.debug(f"{len(spmd_models)-len(spmd_modelsf)} doublons")


# ------------------- Preparing metabolites needed for SPRMS in case models don't contain them -------------------
from cobra import Reaction, Metabolite

S_Adenosylmethioninamine = Metabolite(id="ametam_c",name="S-Adenosylmethioninamine",compartment="c")
_5Methylthioadenosine = Metabolite(id="5mta_c",name="5-Methylthioadenosine",compartment="c")
Spermine_c = Metabolite(id="spm_c",name="spermine",compartment="c")
Spermine_e = Metabolite(id="spm_e",name="spermine",compartment="e")

Heterologuous_metabolites = {
    "ametam_c": S_Adenosylmethioninamine,
    "5mta_c": _5Methylthioadenosine,
    "spm_c": Spermine_c,
    "spm_e":Spermine_e
}

not_prod = 0
no_dev = 0
model_data = dict()

for XMLmodel in spmd_modelsf[:]:

# ------------------- if the model isn't in the "/GSM" folder, download it from BiGG -------------------
    if f"GSM/{XMLmodel}.xml" not in glob.glob("GSM/*.xml"):
        url = f"http://bigg.ucsd.edu/static/models/{XMLmodel}.xml"
        r = requests.get(url, allow_redirects=True)
        open(f"GSM/{XMLmodel}.xml", 'wb').write(r.content)
        logging.debug(f"{XMLmodel} download is done")

    model = read_sbml_model(f"GSM/{XMLmodel}.xml")
    model_metabolites_id = [x.id for x in model.metabolites]

    if "spmd_c" not in model_metabolites_id: # Shouldn't append if BiGG model list is reliable 
        logging.debug(f"removed {model}")
        spmd_modelsf.remove(model)
        continue

    # ------------------- Tests if the model contains needed Metabolites, if not adds them -------------------
    for id in Heterologuous_metabolites.keys():
        if id not in model_metabolites_id:
            model.add_metabolites([Heterologuous_metabolites[id]])
            logging.debug(f"Added {id} to {model}")
    
    # ------------------- Tests if the model contains needed Reactions, if not adds them -------------------
    try:
        r = model.reactions.EX_sprm_e
    except AttributeError:
        EX_sprm_e = Reaction("EX_sprm_e")
        EX_sprm_e.add_metabolites({
            model.metabolites.spm_e:-1,
        })
        EX_sprm_e.bounds = -1000,1000
        model.add_reaction(EX_sprm_e)
        logging.debug(f"{XMLmodel} EX_sprm_e added")

    try:
        r = model.reactions.SPRMS
    except AttributeError:
        SPRMS = Reaction("SPRMS")
        SPRMS.add_metabolites({
            model.metabolites.ametam_c:-1,
            model.metabolites.spmd_c: -1,

            model.metabolites.get_by_id("5mta_c"): 1,
            model.metabolites.h_c:1,
            model.metabolites.spm_c:1
        })
        SPRMS.bounds = -1000,1000
        model.add_reaction(SPRMS)
        logging.debug(f"{XMLmodel} SPRMS added")
        model.add_boundary(model.metabolites.spm_c, type = "demand")


    # ------------------- Finding the Biomass equation name (needed for quadratic objective) -------------------
    biomass_eq = "null"
    for reaction in model.reactions:
        if "biomass" in reaction.name.lower():
            biomass_eq = reaction.id

    if biomass_eq == "null": # Might happen in some model 
        logging.debug(f"No biomass equation found: {XMLmodel}, {model.objective=}")
        continue

    try: 
        r = model.reactions.SPMS
    except AttributeError: # if Spermidine synthase reaction is not present, stop this model analysis here 
        logging.debug(f"{XMLmodel} no SPMS")
        not_prod += 1
        continue

    # ------------- Determining weighting ratio (alpha in the report ) -------------
    max_growth_rate = model.optimize().objective_value
    model.objective = model.reactions.SPRMS
    max_spm_rate = model.optimize().objective_value

    try:
        ratio = max_spm_rate/max_growth_rate
    except ZeroDivisionError: # Happens if WT model has no growth rate = bad model
        logging.debug(f"{XMLmodel} no growth -> out")
        continue
    
    #------------- Setting the quadratic objective -------------
    quadratic_objective = model.problem.Objective(
    1 * model.reactions.SPRMS.flux_expression + ratio * model.reactions.get_by_id(biomass_eq).flux_expression,
    direction='max')
    model.objective = quadratic_objective
    print(f"{XMLmodel} {model.objective}")
    
    #------------- Performs FBA -------------
    solution = model.optimize(objective_sense=None)

    flux1 = solution.fluxes["SPMS"]
    flux2 = solution.fluxes["SPRMS"]
    growth_rate = solution.fluxes[biomass_eq]
    logging.info(f"{XMLmodel} SPMS -> {flux1}, SPRMS -> {flux2}, {growth_rate=}")
   
    if growth_rate > 0: # Verify that growth and Spermine production are compatible
        model_data[XMLmodel] = (flux2,growth_rate,0)
    else:
        no_dev += 1

with open("spmd_models.txt","w") as f:
    f.writelines([f"{x}\n" for x in spmd_modelsf])

logging.info(f"{not_prod/len(spmd_modelsf)*100} % not producer of spmd_c")
logging.info(f"{no_dev/len(spmd_modelsf)*100} % no LP with growth > 0")

model_data = dict(sorted(model_data.items(), key=lambda item: item[1][0], reverse = True))	

#mmol / [gDW h] (concentration per gram dry weight of cells and hour).

# ------------- Store ranked results in output file "GSM_choice.csv" ------------- 
with open("GSM_choice.csv","w") as file:
    writer = csv.writer(file)
    writer.writerow(["GSM","SPRMS flux","growth rate","Memote score"])
    for key,value in model_data.items():
        writer.writerow([key,value[0],value[1],value[2]])
    
