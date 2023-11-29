import numpy as np
import pandas as pd


class DataframeAnalysis:
    """Class for analyzing dataframes."""
    def get_init_mbx_idx(self, df: pd.DataFrame) -> int:
        """Get the index of the first numerical column in a dataframe.

        Args:
            df (pd.DataFrame): The input dataframe.

        Returns:
            int: The index of the first numerical column.

        Raises:
            ValueError: If no numerical columns are found.
        """
        for index, (col, dtype) in enumerate(df.dtypes.items()):
            if np.issubdtype(dtype, np.number):
                return index
        raise ValueError("No numerical columns found.")
