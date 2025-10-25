import pytest
from genotyper.SmartString import SmartString

@pytest.fixture
def repeat_units(sample_settings):
    return sample_settings.get("grouping_units", ["CAG"])

def test_smart_string_simple_repeat(repeat_units):
    seq = "CAGCAGCAGAAA"
    smart = SmartString.get_smart_string_from_sequence(seq, 3, repeat_units)
    assert smart == "[CAG]3AAA", "Should group three CAG repeats and trailing AAA"

def test_smart_string_with_interruption(repeat_units):
    seq = "CAGCAGTTTCAGCAG"
    smart = SmartString.get_smart_string_from_sequence(seq, 3, repeat_units)
    assert smart == "[CAG]2TTT[CAG]2", "Should handle interruption TTT"

def test_smart_string_multiple_units():
    seq = "CAGCAGCAA"
    smart = SmartString.get_smart_string_from_sequence(seq, 3, ["CAG", "CAA"])
    assert smart == "[CAG]2CAA", "Should handle multiple unit types"

def test_smart_string_no_matches(repeat_units):
    seq = "TTTTTT"
    smart = SmartString.get_smart_string_from_sequence(seq, 3, repeat_units)
    assert smart == "[TTT]2", "Should group non-matching sequence as is"

def test_smart_string_empty_sequence(repeat_units):
    seq = ""
    smart = SmartString.get_smart_string_from_sequence(seq, 3, repeat_units)
    assert smart == "", "Empty sequence should return empty string"

def test_smart_string_invalid_index(repeat_units):
    seq = "CAGCAGCAGAAA"
    out = SmartString.get_smart_string_from_sequence(seq, 100, repeat_units)
    assert out == seq

def test_smart_string_invalid_units():
    seq = "CAGCAGCAGAAA"
    out = SmartString.get_smart_string_from_sequence(seq, 3, ["XYZ"])  # Non-DNA unit
    assert out == "[CAG]3AAA"

@pytest.mark.parametrize("seq, window_size, units, expected", [
    ("CAGCAGCAG", 3, ["CAG"], "[CAG]3"),
    ("AAACAGCAGCAG", 3, ["CAG"], "AAA[CAG]3"),
    ("CAGCAGTTTCAG", 3, ["CAG"], "[CAG]2TTTCAG"),
    ("CAGCAGCAA", 3, ["CAG", "CAA"], "[CAG]2CAA"),
    ("TTTT", 3, ["CAG"], "TTTT")
])
def test_smart_string_variations(repeat_units, seq, window_size, units, expected):
    units = units or repeat_units
    smart = SmartString.get_smart_string_from_sequence(seq, window_size, units)
    assert smart == expected, f"Expected {expected} for sequence {seq} with units {units}"