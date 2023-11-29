class DisplayUtils:
    def print_logo(self, tool: str, tool_description: str, version: str):
        """Print the logo, tool name, and version for the tool.

        Parameters:
        tool : str
            The name of the tool.
        tool_description : str
            The description of the tool.
        version : str
            The version of the tool.
        """
        logo = r"""
        ____      _          _       _____           _     
        |__  /    | |    __ _| |__   |_   _|__   ___ | |___ 
        / / ____| |   / _` | '_ \    | |/ _ \ / _ \| / __|
        / / |____| |__| (_| | |_) |   | | (_) | (_) | \__ \
        /____|    |_____\__,_|_.__/    |_|\___/ \___/|_|___/
                                                            
        """

        tool_name = f"{tool} ({version})\n{tool_description}"

        output = f"{'#'*80}\n{logo}\n{tool_name}\n\n{'#'*80}\n"

        print(output)
