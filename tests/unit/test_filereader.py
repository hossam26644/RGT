import pytest
from filereader.ReadFile import ReadFile

def make_fastq(contents, path):
    with open(path, 'w') as f:
        f.write(contents)

def test_extract_reads_between_flanks_exact_match(tmp_path):
    seq = "AAAACCATGTTTTGACCCCG"
    fastq = f"@r1\n{seq}\n+\nIIII\n"
    p = tmp_path / "test1.fastq"
    make_fastq(fastq, p)
    rf = ReadFile(str(p), start_flank="AAAA", end_flank="CCCC",
                  number_of_allowed_flank_point_mutations=0,
                  discard_reads_with_no_end_flank=True)
    assert rf.reads == ['CCATGTTTTGA'], "Should extract sequence between exact flanks"

def test_extract_reads_with_allowed_mismatch(tmp_path):
    seq = "AAATCCATGTTTTGACCCCG"  # One mismatch in start flank (A->T at pos 3)
    fastq = f"@r2\n{seq}\n+\nIIII\n"
    p = tmp_path / "test2.fastq"
    make_fastq(fastq, p)
    rf = ReadFile(str(p), start_flank="AAAA", end_flank="CCCC",
                  number_of_allowed_flank_point_mutations=1,
                  discard_reads_with_no_end_flank=True)
    assert rf.reads == ['CCATGTTTTG'], "Should extract with one mismatch in start flank"

def test_discard_reads_when_no_start_flank(tmp_path):
    seq = "TTTTGACCCCG"
    fastq = f"@r3\n{seq}\n+\nIIII\n"
    p = tmp_path / "test3.fastq"
    make_fastq(fastq, p)
    rf = ReadFile(str(p), start_flank="AAAA", end_flank="CCCC",
                  number_of_allowed_flank_point_mutations=0,
                  discard_reads_with_no_end_flank=True)
    assert rf.reads == [], "Should discard read with no start flank"
    assert rf.number_of_discarded_reads == 1, "Should count one discarded read"

def test_discard_reads_when_no_end_flank(tmp_path):
    seq = "AAAACCATGTTTTGA"
    fastq = f"@r4\n{seq}\n+\nIIII\n"
    p = tmp_path / "test4.fastq"
    make_fastq(fastq, p)
    rf = ReadFile(str(p), start_flank="AAAA", end_flank="CCCC",
                  number_of_allowed_flank_point_mutations=0,
                  discard_reads_with_no_end_flank=True)
    assert rf.reads == [], "Should discard read with no end flank"
    assert rf.number_of_discarded_reads == 1, "Should count one discarded read"

def test_invalid_fastq(tmp_path):
    seq = "AAAACCATGTTTTGACCCCG"
    fastq = f"@r5\n{seq}\n"  # No quality line
    p = tmp_path / "test5.fastq"
    make_fastq(fastq, p)
    rf = ReadFile(str(p), start_flank="AAAA", end_flank="CCCC",
                  number_of_allowed_flank_point_mutations=0,
                  discard_reads_with_no_end_flank=True)
    assert rf.reads == ['CCATGTTTTGA']

def test_empty_fastq(tmp_path):
    fastq = ""
    p = tmp_path / "test6.fastq"
    make_fastq(fastq, p)
    rf = ReadFile(str(p), start_flank="AAAA", end_flank="CCCC",
                  number_of_allowed_flank_point_mutations=0,
                  discard_reads_with_no_end_flank=True)
    assert rf.reads == [], "Empty file should return no reads"
    assert rf.number_of_discarded_reads == 0, "No reads to discard"

def test_discard_percentage(tmp_path):
    fastq = "@r1\nAAAACCATGTTTTGACCCCG\n+\nIIII\n@r2\nTTTTGACCCCG\n+\nIIII\n"
    p = tmp_path / "test7.fastq"
    make_fastq(fastq, p)
    rf = ReadFile(str(p), start_flank="AAAA", end_flank="CCCC",
                  number_of_allowed_flank_point_mutations=0,
                  discard_reads_with_no_end_flank=True)
    rf.reads  # Trigger extraction
    assert rf.get_discarded_reads_percentage() == 50.0, "50% of reads discarded (1/2)"

@pytest.mark.parametrize("seq, mutations, expected_reads", [
    ("AAAACCATGTTTTGACCCCG", 0, ["CCATGTTTTGA"]),  # Exact match
    ("AAATCCATGTTTTGACCCCG", 1, ["CCATGTTTTG"]),  # One mismatch
    ("AATTCCATGTTTTGACCCCG", 1, []),  # Two mismatches
])
def test_flank_mismatches(tmp_path, seq, mutations, expected_reads):
    fastq = f"@r1\n{seq}\n+\nIIII\n"
    p = tmp_path / "test_flank.fastq"
    make_fastq(fastq, p)
    rf = ReadFile(str(p), start_flank="AAAA", end_flank="CCCC",
                  number_of_allowed_flank_point_mutations=mutations,
                  discard_reads_with_no_end_flank=True)
    assert rf.reads == expected_reads, f"Expected {expected_reads} for sequence {seq}"

def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        ReadFile("nonexistent.fastq", start_flank="AAAA", end_flank="CCCC",
                 number_of_allowed_flank_point_mutations=0,
                 discard_reads_with_no_end_flank=True)