from gim import genbank
from gim.extract import find_in_fields, find_in_text


def refs(text):
    return [(r.site, r.obs_id, r.kind) for r in find_in_text(text)]


def test_inat_forms():
    assert refs("isolate S.D. Russell ONT iNaturalist # 98565406 small") == [("inat", 98565406, "inat")]
    assert refs("isolate OMDL iNat # 195601397") == [("inat", 195601397, "inat")]
    assert refs("voucher iNat69337073 small") == [("inat", 69337073, "inat")]
    assert refs("https://www.inaturalist.org/observations/144170346") == [("inat", 144170346, "inat")]


def test_mo_forms():
    assert refs("voucher MushroomObserver.org/68404 internal") == [("mo", 68404, "mo_explicit")]
    assert refs("isolate Mushroom Observer 338435") == [("mo", 338435, "mo_explicit")]
    assert refs("MUOB:123456") == [("mo", 123456, "mo_explicit")]
    assert refs("https://mushroomobserver.org/obs/690633") == [("mo", 690633, "mo_explicit")]
    assert refs("Mushroom Observer MO 123456") == [("mo", 123456, "mo_explicit")]
    for s in ("voucher MO411618 small", "mo411618", "MO_411618", "MO # 411618", "MO-411618", "MO obs. #411618"):
        assert refs(s) == [("mo", 411618, "mo_prefix")], s


def test_rejects():
    for s in (
        "Nectria haematococca strain Mo1792-1",
        "isolate DAMO 171046 18S",
        "Waitea circinata 5.8S",
        "specimen MO:2345678",  # 7 digits: Missouri herbarium style, not MO
        "MO 1234",
        "iNat 2021 project",
        "coinat 123456",
        "Box 8501, St. Louis, MO 63108, USA Tel: 314 286 1800",
        "4566 Scott Ave., St. Louis,MO 63110",
        "Cape Girardeau MO 63701-4799",
        "Morchella sp. strain MP-Mo-090514",
    ):
        assert refs(s) == [], s


GB = """\
LOCUS       OL342397                 650 bp    DNA     linear   PLN 01-JAN-2022
DEFINITION  Tricholoma sp. voucher DUKE:0351597 internal transcribed spacer
            region, partial sequence.
ACCESSION   OL342397
VERSION     OL342397.1
SOURCE      Tricholoma sp.
  ORGANISM  Tricholoma sp.
            Eukaryota; Fungi; Dikarya.
FEATURES             Location/Qualifiers
     source          1..650
                     /organism="Tricholoma sp."
                     /specimen_voucher="DUKE:0351597"
                     /note="Survey-FunDiS; Oconee State Park Foray;
                     MUOB:255212; see https://www.inaturalist.org/observations/
                     98565406"
     CDS             1..10
                     /translation="MOMOMO"
ORIGIN
        1 acgtacgtac
//
"""


def test_parser_and_fields():
    (rec,) = genbank.parse_text(GB)
    assert rec.accession == "OL342397" and rec.version == "OL342397.1"
    assert rec.organism == "Tricholoma sp."
    assert rec.length == 650
    assert rec.qualifier("specimen_voucher") == "DUKE:0351597"
    got = {(r.site, r.obs_id, r.field) for r in find_in_fields(rec.fields)}
    assert got == {("mo", 255212, "source/note"), ("inat", 98565406, "source/note")}
