# Tutorial cobra
from __future__ import print_function

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

if __name__ == '__main__':
    print("Calculando resultados capitulo 1")
    capitulo_1()
