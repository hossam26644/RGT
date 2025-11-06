import pytest
import allelesdetector.PeakIdentifier as PI
from allelesdetector.AllelesDetector import AllelesDetector

@pytest.fixture
def two_peaks_data():
    return (
        {8: 5, 10: 50, 12: 30},
        {
            'seqA': [50, 10, 0, 0, 1, 'rawA'],
            'seqB': [30, 12, 0, 0, 1, 'rawB'],
            'seqC': [5, 8, 0, 0, 1, 'rawC']
        }
    )

@pytest.fixture
def one_peak_data():
    return (
        {9: 10, 10: 100, 11: 60},
        {
            'top': [100, 10, 0, 0, 1, 'rawTop'],
            'n_plus_one': [60, 11, 0, 0, 1, 'rawN+1'],
            'n_minus_one': [10, 9, 0, 0, 1, 'rawN-1']
        }
    )

def test_allelesdetector_two_matches(monkeypatch, two_peaks_data):
    counts_table, geno_table = two_peaks_data
    def fake_get_peaks(self):
        return [10, 12]
    monkeypatch.setattr(PI.PeakIdentifier, 'get_peaks', fake_get_peaks)
    detector = AllelesDetector(counts_table, geno_table, minimum_no_of_reads=1)
    assert detector.color_code == 'green', "Should be green for two peaks"
    assert any('Heterozygous' in summary for summary in detector.result_summery), "Should indicate heterozygous"
    assert detector.first_allele.sequence_string == 'seqA', "First allele should be seqA"
    assert detector.second_allele.sequence_string == 'seqB', "Second allele should be seqB"

def test_allelesdetector_one_match_with_neighbour(monkeypatch, one_peak_data):
    counts_table, geno_table = one_peak_data
    def fake_get_peaks(self):
        return [10]
    monkeypatch.setattr(PI.PeakIdentifier, 'get_peaks', fake_get_peaks)
    detector = AllelesDetector(counts_table, geno_table, minimum_no_of_reads=1)
    assert detector.color_code == 'green', "Should be green for one peak with neighbor"
    assert any('Heterozygous' in summary for summary in detector.result_summery), "Should indicate heterozygous"
    assert detector.first_allele.sequence_string == 'top', "First allele should be top"
    assert detector.second_allele.sequence_string == 'n_plus_one', "Second allele should be n+1"

def test_zero_peaks(monkeypatch, two_peaks_data):
    counts_table, geno_table = two_peaks_data
    def fake_get_peaks(self):
        return []
    monkeypatch.setattr(PI.PeakIdentifier, 'get_peaks', fake_get_peaks)
    detector = AllelesDetector(counts_table, geno_table, minimum_no_of_reads=10)
    assert detector.color_code == 'red', "Should flag with no peaks"
    assert any(('no possible alleles' in str(summary).lower()) or ('please check' in str(summary).lower())
               for summary in detector.result_summery), "Should report no possible alleles or ask to please check"

def test_empty_tables():
    import pytest
    with pytest.raises(IndexError):
        AllelesDetector({}, {}, minimum_no_of_reads=1)

def test_no_valid_neighbor(monkeypatch, one_peak_data):
    counts_table, geno_table = one_peak_data
    counts_table = {10: 100, 13: 10}  # 13 too far from 10
    geno_table = {'top': [100, 10, 0, 0, 1, 'rawTop'], 'other': [10, 13, 0, 0, 1, 'rawOther']}
    def fake_get_peaks(self):
        return [10]
    monkeypatch.setattr(PI.PeakIdentifier, 'get_peaks', fake_get_peaks)
    detector = AllelesDetector(counts_table, geno_table, minimum_no_of_reads=1)
    assert detector.color_code == 'green', "Single peak should pass"
    assert any('Homozygous' in summary for summary in detector.result_summery), "Should report homozygous"
    assert detector.first_allele.sequence_string == 'top', "First allele should be top"
    assert detector.second_allele.sequence_string == 'top', "Second allele should be same as first (homozygous)"

@pytest.mark.parametrize("peaks, expected_color, expected_summary", [
    ([10, 12], 'green', 'Heterozygous'),
    ([10], 'green', 'Homozygous'),
    ([], 'red', 'No possible alleles')
])
def test_peak_scenarios(monkeypatch, two_peaks_data, peaks, expected_color, expected_summary):
    counts_table, geno_table = two_peaks_data
    def fake_get_peaks(self):
        return peaks
    monkeypatch.setattr(PI.PeakIdentifier, 'get_peaks', fake_get_peaks)
    detector = AllelesDetector(counts_table, geno_table, minimum_no_of_reads=1)
    assert detector.color_code == expected_color, f"Expected color {expected_color}"
    if expected_summary.lower() == 'no possible alleles':
        assert any(('no possible alleles' in str(summary).lower()) or ('please check' in str(summary).lower())
                   for summary in detector.result_summery), f"Expected {expected_summary} in summary"
    else:
        assert any(str(summary).find(expected_summary) != -1 for summary in detector.result_summery), f"Expected {expected_summary} in summary"