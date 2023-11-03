def model_info(model_obj, biomass_opt:bool = True):
    if biomass_opt:
        with model_obj:
            bmass_val = str(model_obj.slim_optimize())
    else:
        bmass_val = 'N/A'
    print(f'''
-----------
Name: {model_obj.name}
-----------
Species: {len(model_obj.groups)}
Reactions: {len(model_obj.reactions)}
Metabolites: {len(model_obj.metabolites)}
Genes: {len(model_obj.genes)}
-----------
Constraints: {len(model_obj.constraints)}
Variables: {len(model_obj.variables)}
-----------
Biomass Flux: {bmass_val}
-----------
''')

def get_species(abun_path) -> list:
    from csv import DictReader

    with open(abun_path, 'r') as f:
        read = DictReader(f, delimiter=',')

        species = [row['X'].strip() for row in read]

        return species

def get_diet_info(diet_path):
    from csv import reader

    diet_info = dict()

    with open(diet_path, 'r') as f:
        read = reader(f, delimiter='\t')

        for row in read:
            diet_info[row[0].strip()] = float(row[1])

    return diet_info

def check_constraints(model_obj): #TODO log results
    with model_obj:
        def has_sv():
            sv = [const for const in model_obj.constraints if 'sv_' in const.name]

            if len(sv) > 0 and len(sv) == len(model_obj.metabolites):
                return True
            
        def has_bmass():
            bmass = [const for const in model_obj.constraints if 'bmass' in const.name]

            if len(bmass) == 1:
                return True


        sv_check = has_sv()
        bmass_check = has_bmass()

        if sv_check: print('SV=0 set')
        else: print('SV=0 missing')

        if bmass_check: print('Biomass constrained')
        else: print('Missing biomass constraint')

        if sv_check and bmass_check: return True
        else: return False

def model_setup(model_path, abun_path, diet_path):
    ending = model_path.split('.')[-1]
    name = model_path.split('/')[-1].split('_')[0]
    # print(name)

    if 'json' in ending:
        from cobra.io import load_json_model

        model = load_json_model(model_path)

    elif 'mat' in ending:
        from cobra.io import load_matlab_model

        model = load_matlab_model(model_path)

    else: print('File type currently not supported. Use .json or .mat model')

    if model:
        from optlang.gurobi_interface import Constraint

        model.name = name

        print('Model Loaded')
        model_info(model)

        to_add = []
        #SECTION setting biomass constraint
        biomass = model.reactions.get_by_id('communityBiomass')
        print('Biomass constrained to 0.4 - 1')
        to_add.append(
            Constraint(name='bmass_const', expression=biomass.flux_expression, lb=0.4, ub=1)
        )

        #SECTION setting SV=0
        for met in model.metabolites:
            print(f'({model.metabolites.index(met) + 1} / {len(model.metabolites)})\tNow adding Sv=0 for {met.id}')
            rxns = [rxn.flux_expression for rxn in met.reactions]
            coefs = [v.get_coefficient(met) for v in met.reactions]

            mulled = []
            for i,j in zip(rxns, coefs):
                mulled.append(i*j)

            expr = sum(mulled)

            same_flux = Constraint(
                name='sv_' + met.id,
                expression=expr,
                lb=0, ub=0
            )

            to_add.append(same_flux)
        
        #SECTION add species groups
        from cobra.core import Group

        species = get_species(abun_path)
        groups = []

        for spec in species:
            print(f'{species.index(spec) + 1} / {len(species)}\t Grouping {spec}')
            rxns = [rxn for rxn in model.reactions if spec in rxn.id]

            if len(rxns) != 0:
                spec_group = Group(
                    id=spec,
                    members=rxns,
                    kind='collection'
                )

                groups.append(spec_group)

        #SECTION updating model
        model.add_cons_vars(to_add)
        model.add_groups(groups)

        model.solver.update()

        model_info(model)

        #SECTION setting diet #TODO automatic diet selection based on sample location information
        print('Setting diet')

        diet = get_diet_info(diet_path)

        diet_rxns = [rxn for rxn in model.reactions if '[d]' in rxn.id]

        for rxn in diet_rxns:
            rxn.lower_bound = 0

        rxns = [rxn for rxn in diet_rxns if rxn.id in diet.keys()]
        
        for rxn in rxns:
            rxn.lower_bound = diet[rxn.id]


        model_info(model)

        #SECTION return model
        return model