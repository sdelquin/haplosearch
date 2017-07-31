import os
from tempfile import mktemp
import filecmp
from django.conf import settings
from app.haploutils import manage_haplosearch


TESTFILES_PATH = os.path.join(settings.BASE_DIR, "app/static/files")
# test file with haplotypes in nomenclature forensic genetics
TESTFILE_HAP_FOR = os.path.join(TESTFILES_PATH,
                                "mtDNA.haplotypes.forensic.txt")
# test file with haplotypes in nomenclature population genetics
TESTFILE_HAP_POP = os.path.join(TESTFILES_PATH,
                                "mtDNA.haplotypes.population.txt")
# test file with sequences
TESTFILE_SEQ = os.path.join(TESTFILES_PATH, "mtDNA.sequences.txt")


def test_hapfor_to_seq():
    outputfile_path = mktemp()
    manage_haplosearch(TESTFILE_HAP_FOR, outputfile_path, "FOR", "H2S")
    cmp = filecmp.cmp(TESTFILE_SEQ, outputfile_path)
    os.remove(outputfile_path)
    assert cmp


def test_happop_to_seq():
    outputfile_path = mktemp()
    manage_haplosearch(TESTFILE_HAP_POP, outputfile_path, "POP", "H2S")
    cmp = filecmp.cmp(TESTFILE_SEQ, outputfile_path)
    os.remove(outputfile_path)
    assert cmp


def test_seq_to_hapfor():
    outputfile_path = mktemp()
    manage_haplosearch(TESTFILE_SEQ, outputfile_path, "FOR", "S2H")
    cmp = filecmp.cmp(TESTFILE_HAP_FOR, outputfile_path)
    os.remove(outputfile_path)
    assert cmp


def test_seq_to_happop():
    outputfile_path = mktemp()
    manage_haplosearch(TESTFILE_SEQ, outputfile_path, "POP", "S2H")
    cmp = filecmp.cmp(TESTFILE_HAP_POP, outputfile_path)
    os.remove(outputfile_path)
    assert cmp
