import re

class StringProcessing:
    def convert_string(self, s: str) -> str:
        """Convert a string to a standard format for matching."""
        s = re.sub(r"\s*\(.*?\)", "", s)
        s = re.sub(r"(\d)\.(\d)", r"\1,\2", s)
        s = re.sub(r"(\d) (\w)", r"\1-\2", s)
        return s