from os import listdir
from basic_mods import setup_model
from min_diet import min_diet_finder
from optlang.gurobi_interface import Constraint

mod_dir = 'Resources/Models/no_diet/'
bmass = 'communityBiomass'

mod_paths = [mod_dir+i for i in listdir(mod_dir)]

def run_setup_checks(model):
    # checks if SV=0 set for all metabolites
    if len([const for const in model.constraints if const.lb == const.ub == 0]) == len(model.metabolites):
        pass
    else:
        raise LookupError
    
    # checks if biomass constraint set
    if len([const for const in model.constraints if const.expression == model.reactions.get_by_id('communityBiomass').flux_expression]) == 1:
        pass
    else:
        raise LookupError

def generate_report(metabs, outfile):
    with open(outfile, 'a+') as f:
        lines = f.readlines()
        for met in metabs:
            if met not in lines:
                f.write(met + '\n')

def run_min_diet(model_path):
    model = setup_model(model_path)

    # adds biomass constraint (0.4 - 1.0)
    model.add_cons_vars(Constraint(model.reactions.get_by_id(bmass).flux_expression, lb=0.4, ub=1))

    run_setup_checks(model)

    vars = min_diet_finder(model) 

    print(len(vars))
    print(vars)

    generate_report(vars, '.tmp/minimal_diet.csv')

paths = mod_paths[2:]

for path in paths:
    print(f'''
------------------
{paths.index(path) + 1} of {len(paths)}
------------------
''')
    # print(path)
    run_min_diet(path)
# for mod_path in [mod_paths[2:]]:
#     print(type(mod_path))l
    # run_min_diet(mod_path)
    # print(mod_path)
