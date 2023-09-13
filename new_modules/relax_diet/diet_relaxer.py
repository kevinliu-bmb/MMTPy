from optlang.gurobi_interface import Model, Constraint, Variable


def min_relax_metab_finder(model: Model, scope: str = "all", relax_value: int = -10):
    """
    Input
        * model: model to be examined
        * scope: either all diet metabolites or only essential diet metabolites to be checked for need for relaxed bounds
        * relax_value: value to try relaxing the bound of the metabolite to

    Output
        list of relaxed metabolites
    """
    essentials = [
        "Diet_EX_12dgr180[d]",
        "Diet_EX_26dap_M[d]",
        "Diet_EX_2dmmq8[d]",
        "Diet_EX_2obut[d]",
        "Diet_EX_3mop[d]",
        "Diet_EX_4abz[d]",
        "Diet_EX_4hbz[d]",
        "Diet_EX_ac[d]",
        "Diet_EX_acgam[d]",
        "Diet_EX_acmana[d]",
        "Diet_EX_acnam[d]",
        "Diet_EX_ade[d]",
        "Diet_EX_adn[d]",
        "Diet_EX_adocbl[d]",
        "Diet_EX_adpcbl[d]",
        "Diet_EX_ala_D[d]",
        "Diet_EX_ala_L[d]",
        "Diet_EX_amet[d]",
        "Diet_EX_amp[d]",
        "Diet_EX_arab_D[d]",
        "Diet_EX_arab_L[d]",
        "Diet_EX_arg_L[d]",
        "Diet_EX_asn_L[d]",
        "Diet_EX_btn[d]",
        "Diet_EX_ca2[d]",
        "Diet_EX_cbl1[d]",
        "Diet_EX_cgly[d]",
        "Diet_EX_chor[d]",
        "Diet_EX_chsterol[d]",
        "Diet_EX_cit[d]",
        "Diet_EX_cl[d]",
        "Diet_EX_cobalt2[d]",
        "Diet_EX_csn[d]",
        "Diet_EX_cu2[d]",
        "Diet_EX_cys_L[d]",
        "Diet_EX_cytd[d]",
        "Diet_EX_dad_2[d]",
        "Diet_EX_dcyt[d]",
        "Diet_EX_ddca[d]",
        "Diet_EX_dgsn[d]",
        "Diet_EX_fald[d]",
        "Diet_EX_fe2[d]",
        "Diet_EX_fe3[d]",
        "Diet_EX_fol[d]",
        "Diet_EX_for[d]",
        "Diet_EX_gal[d]",
        "Diet_EX_glc_D[d]",
        "Diet_EX_gln_L[d]",
        "Diet_EX_glu_L[d]",
        "Diet_EX_gly[d]",
        "Diet_EX_glyc[d]",
        "Diet_EX_glyc3p[d]",
        "Diet_EX_gsn[d]",
        "Diet_EX_gthox[d]",
        "Diet_EX_gthrd[d]",
        "Diet_EX_gua[d]",
        "Diet_EX_h[d]",
        "Diet_EX_h2o[d]",
        "Diet_EX_h2s[d]",
        "Diet_EX_his_L[d]",
        "Diet_EX_hxan[d]",
        "Diet_EX_ile_L[d]",
        "Diet_EX_k[d]",
        "Diet_EX_lanost[d]",
        "Diet_EX_leu_L[d]",
        "Diet_EX_lys_L[d]",
        "Diet_EX_malt[d]",
        "Diet_EX_met_L[d]",
        "Diet_EX_mg2[d]",
        "Diet_EX_mn2[d]",
        "Diet_EX_mqn7[d]",
        "Diet_EX_mqn8[d]",
        "Diet_EX_nac[d]",
        "Diet_EX_ncam[d]",
        "Diet_EX_nmn[d]",
        "Diet_EX_no2[d]",
        "Diet_EX_ocdca[d]",
        "Diet_EX_ocdcea[d]",
        "Diet_EX_orn[d]",
        "Diet_EX_phe_L[d]",
        "Diet_EX_pheme[d]",
        "Diet_EX_pi[d]",
        "Diet_EX_pnto_R[d]",
        "Diet_EX_pro_L[d]",
        "Diet_EX_ptrc[d]",
        "Diet_EX_pydx[d]",
        "Diet_EX_pydxn[d]",
        "Diet_EX_q8[d]",
        "Diet_EX_rib_D[d]",
        "Diet_EX_ribflv[d]",
        "Diet_EX_ser_L[d]",
        "Diet_EX_sheme[d]",
        "Diet_EX_so4[d]",
        "Diet_EX_spmd[d]",
        "Diet_EX_thm[d]",
        "Diet_EX_thr_L[d]",
        "Diet_EX_thymd[d]",
        "Diet_EX_trp_L[d]",
        "Diet_EX_ttdca[d]",
        "Diet_EX_tyr_L[d]",
        "Diet_EX_ura[d]",
        "Diet_EX_val_L[d]",
        "Diet_EX_xan[d]",
        "Diet_EX_xyl_D[d]",
        "Diet_EX_zn2[d]",
    ]

    with model as m:
        if scope == "all":
            rxns = [rxn for rxn in m.reactions if "[d]" in rxn.id]
        elif scope == "essential":
            rxns = [rxn for rxn in m.reactions if rxn.id in essentials]
        else:
            print(
                "The scope can either be all diet metabolites (all) or essential metabolites (essential)"
            )

        to_add = []

        # adds indicators and constraints for either diet or essential metabolites
        for rxn in rxns:
            flux = rxn.flux_expression

            rxn.lower_bound = -1000

            indicator = Variable(name="ind_" + rxn.id, lb=0, ub=1, type="binary")

            constr = Constraint(
                name="relax_const_" + rxn.id,
                expression=flux - (indicator * relax_value + -0.1 * (1 - indicator)),
                lb=0,
            )

            to_add.extend([indicator, constr])

        m.add_cons_vars(to_add)
        # m.solver.update()

        indicators = [v for v in m.variables if "ind_" in v.name]
        expression = sum(indicators)

        # minimizes sum of indicators. indicator = 1 if metabolite is relaxed, = 0 if not
        m.objective = expression
        m.objective_direction = "min"
        m.solver.update()
        m.optimize() # minimizes sum of indicators

        # gets list of only metabolites that need to be relaxed
        relaxed_metabs = [
            v.name.strip("ind_")
            for v in m.variables
            if "ind_" in v.name and v.primal == 1
        ]

        # outputs number of relaxed metabolites to console
        print(f"{len(relaxed_metabs)} metabolites relaxed.")
        # print(relaxed_metabs)

        return relaxed_metabs


def min_mod_adjuster(model: Model, metabolites: list):
    """Gets difference from original bound (alpha) by which it needs to be relaxed. Note: must be a positive value

    Input
        model: model object to be checked
        metabolites: list of metabolites to check for relaxation

    Output
        list of (metabolite, alpha) for input metabolites
    """
    with model as m:
        # getting reaction objects for each metabolite from model
        rxns = [rxn for rxn in m.reactions if rxn.id in metabolites]

        to_add = []

        # adding alphas (change needed for each metabolite bound) and constraints for each input metabolite
        for rxn in rxns:
            flux = rxn.flux_expression

            alpha = Variable(name="alpha_" + rxn.id, lb=0)

            alpha_const = Constraint(
                name="alpha_const_" + rxn.id, expression=flux - (-0.1 + alpha)
            )

            to_add.extend([alpha, alpha_const])

        m.add_cons_vars(to_add)
        m.solver.update()

        alphas = [v for v in m.variables if "alpha_" in v.name]
        expression = sum(alphas)

        # sets the objective to minimize the sum of alphas
        m.objective = expression
        m.objective_direction = "min"
        m.optimize()  # minimizes sum of alphas

        for alpha in alphas:  # outputs results to console
            print(alpha.name, alpha.primal)

        return [(alpha.name.strip("alpha_"), alpha.primal) for alpha in alphas]
