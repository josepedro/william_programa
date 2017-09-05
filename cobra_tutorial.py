# Tutorial cobra
from __future__ import print_function
from cobra import Model, Reaction, Metabolite

import cobra
import cobra.test

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

if __name__ == '__main__':
    print("---------Calculando resultados capitulo 1---------")
    #capitulo_1()
    print("---------Calculando resultados capitulo 2---------")
    capitulo_2()
