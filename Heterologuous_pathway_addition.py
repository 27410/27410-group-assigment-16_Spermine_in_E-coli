from cobra.io import read_sbml_model, write_sbml_model
import logging 

level = logging.INFO	
fmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(level =level, format=fmt)


def main(xml_model = "GSM/iAF1260b.xml"):

    from cobra import Reaction, Metabolite

    model = read_sbml_model(xml_model)

    #-------------- Adding spm_c and spm_e to model's metabolites --------------

    Spermine_c = Metabolite(id="spm_c",name="spermine",compartment="c")
    Spermine_e = Metabolite(id="spm_e",name="spermine",compartment="e")

    for metabolite in [Spermine_c,Spermine_e]:
        if metabolite not in model.metabolites:
            model.add_metabolites([metabolite])
            logging.info(f"{metabolite} metabolite added")

    #-------------- Adding SPRMS + demand reaction --------------

    SPRMS = Reaction("SPRMS")
    SPRMS.add_metabolites({
        model.metabolites.ametam_c:-1,
        model.metabolites.spmd_c: -1,

        model.metabolites.get_by_id("5mta_c"): 1,
        model.metabolites.h_c:1,
        model.metabolites.spm_c:1
    })
        
    if SPRMS not in model.reactions:
        model.add_reaction(SPRMS)
        model.add_boundary(model.metabolites.spm_c, type = "demand")
        logging.info(f"{SPRMS.id} reaction added")

    write_sbml_model(model,xml_model)

def model_assertion(xml_model = "GSM/iAF1260b.xml"):
    
    #------------- Security - Verify that requierements for Spermine production are fullfiled -------------

    model = read_sbml_model(xml_model)

    assert "spm_e" in [x.id for x in model.metabolites]
    assert "spm_c" in [x.id for x in model.metabolites]
    assert "DM_spm_c" in [x.id for x in model.reactions]
    assert "SPRMS" in [x.id for x in model.reactions]


if __name__ == "__main__":
    main()
    model_assertion()

