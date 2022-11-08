#%%
from cobra.io import save_json_model, read_sbml_model
from escher import Builder,list_available_maps

model =  read_sbml_model("GSM/iAF1260.xml")
save_json_model(model,"GSM/iAF1260.json")

print(list_available_maps())
builder = Builder(
    model_json='GSM/iAF1260.json'
)
builder