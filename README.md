# genbank-inat-mo

Finds every GenBank nucleotide record that cites an iNaturalist or Mushroom
Observer observation number, then checks each observation on iNaturalist and
MO: does it exist, does it list the GenBank accession, is it the same genus.

It only reads from NCBI, iNaturalist and MO and never writes to them. It uses
only the Python standard library (3.10+).

```bash
cd ~/code/genbank-inat-mo
export NCBI_API_KEY=...        # optional; 10 requests/s instead of 3
python3 -m gim run             # search + fetch + verify + report
xdg-open report/index.html
```

Or step by step (each step resumes where it stopped and skips finished work):

```bash
python3 -m gim search          # Entrez searches -> candidate records
python3 -m gim fetch           # download candidates as GenBank text, extract numbers
python3 -m gim verify          # look the observations up (iNat ~1 req/s, MO 1 req/5 s)
python3 -m gim report          # report/{summary.txt,links.csv,problems.csv,index.html}
```

State lives in `data/gim.sqlite` (`--db` to change it). `gim verify --max-age 30`
re-checks observations looked up more than 30 days ago; `--refresh` re-checks all.

## What counts as a reference

Matched in every text part of the record except the sequence: definition line,
all feature qualifiers (`/specimen_voucher`, `/isolate`, `/bio_material`,
`/culture_collection`, `/note`, ...), COMMENT and reference titles.

| kind | examples |
|---|---|
| `inat` | `iNaturalist # 98565406`, `iNat69337073`, `iNat:123456`, `inaturalist.org/observations/123456`, `naturalista` |
| `mo_explicit` | `Mushroom Observer 338435`, `MushroomObserver.org/68404`, `MUOB:123456`, `mushroomobserver.org/obs/690633` |
| `mo_prefix` | `MO123456`, `mo555555`, `MO_559147`, `MO # 346490`, `MO-123456`, `MO obs. 123456` |

MO numbers must be 5 or 6 digits; iNat numbers 5 to 10. `MO` is also the
Index Herbariorum code of the Missouri Botanical Garden herbarium, so a bare
`MO` match whose MO observation does not exist, or is a different genus and
lacks the accession, is marked `doubtful` in the report.

## How the search finds records

GenBank has about 23 million fungal records, too many to download. The search
relies on how Entrez indexes words instead:

- `MO # 346490`, `MO_559147`, `MO-559147` and `MUOB:123456` are indexed under
  the word `mo` or `muob` (plus a phrase like "mo 559147"), so keyword
  searches for `mo`, `muob`, `inat`, `inaturalist`, `mushroomobserver`,
  `mushroom AND observer` and `naturalista` find them. `mo` is limited to fungi
  and slime moulds. The other words search every organism.
- A number written straight after the letters (`MO411618`, `iNat69337073`) is
  a word of its own, which only a wildcard finds. Entrez expands a wildcard to
  at most 600 words and says so, so each root is split into longer prefixes
  (`mo4*`, `mo41*`, `mo411*`, ...) until no expansion is cut short.
- WGS contigs are left out (each contig repeats the voucher); WGS master
  records stay in.

The candidates are then downloaded as full GenBank text and matched with the
patterns above, so a candidate that only matched a keyword by accident
produces no reference.

### Checking recall

- `gim fasta FILE.fasta.gz` reads a FASTA file's titles (for example the
  weekly ITS download in `/opt/its-blast/db/current/genbank_its.fasta.gz` on
  the mycomap server) and records every title that cites iNat or MO. The next
  `gim fetch` downloads the searches' records first, then any title record
  no search reached, so nothing is lost, and the report lists those misses
  in `missed_by_search.csv`. (`gim run --fasta FILE` does the same.)
- `gim scan 'QUERY'` downloads and checks *every* record matching an Entrez
  query, with no keyword shortcut. For example,
  `gim scan 'fungi[organism] AND 2024:2026[PDAT] AND 200:3000[SLEN] NOT wgs[filter]'`.
  It is slow (millions of records) but resumable.

## Report statuses

| status | meaning |
|---|---|
| `linked` | the observation lists the accession (iNat observation field, MO sequence record) |
| `mentioned` | the accession is only in notes, description or comments |
| `missing_has_seq` | the observation has a sequence or DNA barcode field but not this accession |
| `missing_no_seq` | the observation has no sequence and no accession |
| `not_found` | the observation does not exist: wrong number in GenBank, or deleted |
| `not_returned` | MO knows the number but returns another observation (merged occurrence) |

`taxon_check` compares genera (`different_genus` is often just a name change,
but also catches a mistyped number). `problems.csv` holds every row that is
not `linked`/`mentioned`, plus genus mismatches.
