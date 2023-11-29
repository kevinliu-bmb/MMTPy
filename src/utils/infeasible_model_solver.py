import cobra

from utils.essential_metabolites import essential_metabolites


class InfeasibleModelSolver:
    def __init__(self, model):
        self.model = model

    def add_slack_variables_to_constraints(self, constraints):
        """Add slack variables to given constraints."""
        slack_constraints = []
        slack_variables = []
        for constraint in constraints:
            slack_variable_pos = self.model.problem.Variable(
                constraint.name + "_slack_pos", lb=0
            )
            slack_variable_neg = self.model.problem.Variable(
                constraint.name + "_slack_neg", lb=0
            )

            slack_variables.extend([slack_variable_pos, slack_variable_neg])

            new_constraint = self.model.problem.Constraint(
                constraint.expression + slack_variable_pos - slack_variable_neg,
                lb=constraint.lb,
                ub=constraint.ub,
                name=constraint.name + "_slack",
            )

            slack_constraints.append(new_constraint)

        return slack_constraints, slack_variables

    def solve_with_slack_variables(self, constraints):
        """Solve the model with added slack variables to the constraints."""
        slack_constraints, slack_variables = self.add_slack_variables_to_constraints(
            constraints
        )

        self.model.add_cons_vars(slack_constraints)
        self.model.solver.update()

        # Set the objective to be the sum of the absolute values of the slack variables
        self.model.objective = self.model.problem.Objective(sum(slack_variables))

        solution = self.model.optimize(objective_sense="minimize")

        # Remove the slack variables and constraints from the model
        self.model.remove_cons_vars(slack_variables)
        self.model.remove_cons_vars(slack_constraints)
        self.model.solver.update()

        return solution

    # Additional methods can be added here as needed
