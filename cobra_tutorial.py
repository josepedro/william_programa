# Tutorial cobra
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
from os.path import join
from time import time
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

import cobra
import cobra.test
import os
import pandas

def capitulo_1():
    file = open("resultados_capitulo_1.txt","w") 
    model = cobra.test.create_test_model("textbook")
    file.write(str(len(model.reactions)))
    file.write("\n")
    file.write(str(len(model.metabolites)))
    file.write("\n")
    file.write(str(len(model.genes)))
    file.write("\n")
    model
    model.reactions[29]
    model.metabolites.get_by_id("atp_c")
    model.reactions.EX_glc__D_e.bounds
    pgi = model.reactions.get_by_id("PGI")
    pgi
    file.write(pgi.name)
    file.write("\n")
    file.write(pgi.reaction)
    file.write("\n")
    file.write(str(pgi.lower_bound) + "< pgi <" + str(pgi.upper_bound))
    file.write("\n")
    pgi.check_mass_balance()
    pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
    pgi.reaction
    pgi.check_mass_balance()
    pgi.subtract_metabolites({model.metabolites.get_by_id("h_c"): -1})
    file.write(pgi.reaction)
    file.write("\n")
    pgi.reaction = "g6p_c --> f6p_c + h_c + green_eggs + ham"
    pgi.reaction
    pgi.reaction = "g6p_c <=> f6p_c"
    pgi.reaction
    atp = model.metabolites.get_by_id("atp_c")
    atp
    file.write(atp.name)
    file.write("\n")
    file.write(atp.compartment)
    file.write("\n")
    atp.charge
    file.write(atp.formula)
    file.write("\n")
    len(atp.reactions)
    model.metabolites.get_by_id("g6p_c").reactions
    gpr = pgi.gene_reaction_rule
    gpr
    pgi.genes
    pgi_gene = model.genes.get_by_id("b4025")
    pgi_gene
    pgi_gene.reactions
    pgi.gene_reaction_rule = "(spam or eggs)"
    pgi.genes
    pgi_gene.reactions
    model.genes.get_by_id("spam")
    cobra.manipulation.delete_model_genes(
        model, ["spam"], cumulative_deletions=True)
    file.write("after 1 KO: %4d < flux_PGI < %4d" % (pgi.lower_bound, pgi.upper_bound))
    file.write("\n")
    cobra.manipulation.delete_model_genes(
        model, ["eggs"], cumulative_deletions=True)
    file.write("after 2 KO: %4d < flux_PGI < %4d" % (pgi.lower_bound, pgi.upper_bound))
    file.write("\n")
    cobra.manipulation.undelete_model_genes(model)
    model = cobra.test.create_test_model('textbook')
    for reaction in model.reactions[:5]:
        with model as model:
            reaction.knock_out()
            model.optimize()
            file.write('%s blocked (bounds: %s), new growth rate %f' %
                (reaction.id, str(reaction.bounds), model.objective.value))
            file.write("\n")
    [reaction.bounds for reaction in model.reactions[:5]]
    with model:
        model.objective = 'ATPM'
        with model:
            model.objective = 'ACALD'
    file.close()

def capitulo_2():
    file = open("resultados_capitulo_2.txt","w")
    model = Model('example_model')
    reaction = Reaction('3OAS140')
    reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
    reaction.subsystem = 'Cell Envelope Biosynthesis'
    reaction.lower_bound = 0. # This is the default
    reaction.upper_bound = 1000. # This is the default
    ACP_c = Metabolite(
        'ACP_c',
        formula='C11H21N2O7PRS',
        name='acyl-carrier-protein',
        compartment='c')
    omrsACP_c = Metabolite(
        '3omrsACP_c',
        formula='C25H45N2O9PRS',
        name='3-Oxotetradecanoyl-acyl-carrier-protein',
        compartment='c')
    co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
    malACP_c = Metabolite(
        'malACP_c',
        formula='C14H22N2O10PRS',
        name='Malonyl-acyl-carrier-protein',
        compartment='c')
    h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
    ddcaACP_c = Metabolite(
        'ddcaACP_c',
        formula='C23H43N2O8PRS',
        name='Dodecanoyl-ACP-n-C120ACP',
        compartment='c')
    reaction.add_metabolites({
        malACP_c: -1.0,
        h_c: -1.0,
        ddcaACP_c: -1.0,
        co2_c: 1.0,
        ACP_c: 1.0,
        omrsACP_c: 1.0
    })
    reaction.reaction
    reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
    reaction.genes
    file.write('%i reactions initially' % len(model.reactions))
    file.write("\n")
    file.write('%i metabolites initially' % len(model.metabolites))
    file.write("\n")
    file.write('%i genes initially' % len(model.genes))
    file.write("\n")
    model.add_reactions([reaction])
    file.write('%i reaction' % len(model.reactions))
    file.write("\n")
    file.write('%i metabolites' % len(model.metabolites))
    file.write("\n")
    file.write('%i genes' % len(model.genes))
    file.write("\n")
    file.write("Reactions")
    file.write("\n")
    file.write("---------")
    file.write("\n")
    for x in model.reactions:
        file.write("%s : %s" % (x.id, x.reaction))
        file.write("\n")
    file.write("")
    file.write("\n")
    file.write("Metabolites")
    file.write("\n")
    file.write("-----------")
    file.write("\n")
    for x in model.metabolites:
        file.write('%9s : %s' % (x.id, x.formula))
        file.write("\n")
    file.write("")
    file.write("\n")
    file.write("Genes")
    file.write("\n")
    file.write("-----")
    file.write("\n")
    for x in model.genes:
        associated_ids = (i.id for i in x.reactions)
        file.write("%s is associated with reactions: %s" %
            (x.id, "{" + ", ".join(associated_ids) + "}"))
        file.write("\n")
    model.objective = '3OAS140'
    file.close()

def capitulo_3():
    data_dir = cobra.test.data_dir

    print("mini test files: ")
    print(", ".join(i for i in os.listdir(data_dir) if i.startswith("mini")))

    textbook_model = cobra.test.create_test_model("textbook")
    ecoli_model = cobra.test.create_test_model("ecoli")
    salmonella_model = cobra.test.create_test_model("salmonella")
    cobra.io.read_sbml_model(join(data_dir, "mini_fbc2.xml"))
    cobra.io.write_sbml_model(textbook_model, "test_fbc2.xml")
    cobra.io.read_sbml_model(join(data_dir, "mini_cobra.xml"))
    cobra.io.write_sbml_model(
    textbook_model, "test_cobra.xml", use_fbc_package=False)

def capitulo_4():
    file = open("resultados_capitulo_4.txt","w")
    model = cobra.test.create_test_model("textbook")
    solution = model.optimize()
    print(solution)
    file.write(str(solution.objective_value)); file.write("\n")
    file.write(str(model.optimize().objective_value)); file.write("\n")
    file.write(str(model.slim_optimize())); file.write("\n")
    file.write(str(model.summary())); file.write("\n")
    file.write(str(model.metabolites.nadh_c.summary())); file.write("\n")
    file.write(str(model.metabolites.atp_c.summary())); file.write("\n")
    biomass_rxn = model.reactions.get_by_id("Biomass_Ecoli_core")
    from cobra.util.solver import linear_reaction_coefficients
    file.write(str(linear_reaction_coefficients(model))); file.write("\n")
    model.objective = "ATPM"
    model.reactions.get_by_id("ATPM").upper_bound = 1000.
    file.write(str(linear_reaction_coefficients(model))); file.write("\n")
    file.write(str(model.optimize().objective_value)); file.write("\n")
    from cobra.flux_analysis import flux_variability_analysis
    file.write(str(flux_variability_analysis(model, model.reactions[:10]))); file.write("\n")
    file.write(str(cobra.flux_analysis.flux_variability_analysis(
        model, model.reactions[:10], fraction_of_optimum=0.9)))
    file.write("\n")
    loop_reactions = [model.reactions.FRD7, model.reactions.SUCDi]
    file.write(str(flux_variability_analysis(model, reaction_list=loop_reactions, loopless=False))); file.write("\n")
    file.write(str(flux_variability_analysis(model, reaction_list=loop_reactions, loopless=True))); file.write("\n")
    file.write(str(model.optimize())); file.write("\n")
    file.write(str(model.summary(fva=0.95))); file.write("\n")
    file.write(str(model.metabolites.pyr_c.summary(fva=0.95))); file.write("\n")
    model.objective = 'Biomass_Ecoli_core'
    fba_solution = model.optimize()
    pfba_solution = cobra.flux_analysis.pfba(model)
    file.write(str(abs(fba_solution.fluxes["Biomass_Ecoli_core"] - pfba_solution.fluxes[
        "Biomass_Ecoli_core"])))
    file.write("\n")
    file.close()

def capitulo_5():
    file = open("resultados_capitulo_5.txt","w")
    cobra_model = cobra.test.create_test_model("textbook")
    ecoli_model = cobra.test.create_test_model("ecoli")
    file.write(str(cobra_model.optimize())); file.write("\n")
    with cobra_model:
        cobra_model.reactions.PFK.knock_out()
        file.write(str(cobra_model.optimize())); file.write("\n")
    file.write(str(cobra_model.optimize())); file.write("\n")
    with cobra_model:
        cobra_model.genes.b1723.knock_out()
        file.write(str(cobra_model.optimize())); file.write("\n")
        cobra_model.genes.b3916.knock_out()
        file.write(str(cobra_model.optimize())); file.write("\n")
    deletion_results = single_gene_deletion(cobra_model)
    single_gene_deletion(cobra_model, cobra_model.genes[:20])
    single_reaction_deletion(cobra_model, cobra_model.reactions[:20])
    double_gene_deletion(
        cobra_model, cobra_model.genes[-5:], return_frame=True).round(4)
    start = time() 
    double_gene_deletion(
        ecoli_model, ecoli_model.genes[:300], number_of_processes=2)
    t1 = time() - start
    file.write("Double gene deletions for 200 genes completed in " "%.2f sec with 2 cores" % t1); file.write("\n")
    start = time() 
    double_gene_deletion(
        ecoli_model, ecoli_model.genes[:300], number_of_processes=1)
    t2 = time() - start
    file.write("Double gene deletions for 200 genes completed in " "%.2f sec with 1 core" % t2); file.write("\n")
    file.write("Speedup of %.2fx" % (t2 / t1)); file.write("\n")
    double_reaction_deletion(
        cobra_model, cobra_model.reactions[2:7], return_frame=True).round(4)

    file.close()



if __name__ == '__main__':
    print("---------Calculando resultados capitulo 1---------")
    #capitulo_1()
    print("---------Calculando resultados capitulo 2---------")
    #capitulo_2()
    print("---------Calculando resultados capitulo 3---------")
    #capitulo_3()
    print("---------Calculando resultados capitulo 4---------")
    #capitulo_4()
    print("---------Calculando resultados capitulo 5---------")
    capitulo_5()