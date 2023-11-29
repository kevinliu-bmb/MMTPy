import sys
import pytest

sys.path.append("src/utils")

from string_processing import StringProcessing

def test_convert_string_converts_dots_to_commas_between_numbers():
    processor = StringProcessing()
    assert processor.convert_string("3.5") == "3,5"

def test_convert_string_adds_hyphen_between_number_and_word():
    processor = StringProcessing()
    assert processor.convert_string("1 mg") == "1-mg"

def test_convert_string_leaves_string_unchanged_when_no_pattern_matches():
    processor = StringProcessing()
    assert processor.convert_string("Metabolite") == "Metabolite"
