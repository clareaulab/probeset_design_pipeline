import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from utils import MismatchScanner  # noqa: E402

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def count_mismatches(query: str, target: str) -> int:
    """The reference rule: an 'N' on either side matches, ACGT must agree, anything else is a mismatch."""
    mismatches = 0
    for q, t in zip(query.upper(), target.upper()):
        if q == "N" or t == "N":
            continue
        if q != t or q not in "ACGT":
            mismatches += 1
    return mismatches


def brute_force(background: dict[str, str], query: str, max_mismatches: int) -> dict[str, int]:
    best = {}
    for name, sequence in background.items():
        gene = name.split(" ")[0]
        for strand_query in (query, reverse_complement(query)):
            for start in range(len(sequence) - len(query) + 1):
                mm = count_mismatches(strand_query, sequence[start:start + len(query)])
                if mm <= max_mismatches and mm < best.get(gene, max_mismatches + 1):
                    best[gene] = mm
    return best


def random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(rng: random.Random, sequence: str, n: int) -> str:
    sequence = list(sequence)
    for i in rng.sample(range(len(sequence)), n):
        sequence[i] = rng.choice([b for b in "ACGT" if b != sequence[i]])
    return "".join(sequence)


TARGET = "GATTACAGATTACACCGGTTAAGGCCTTAAGG"  # 32 nt
QUERY = TARGET[3:28]  # 25 nt, an exact match


def scan(background: dict[str, str], query: str, max_mismatches: int = 4) -> dict[str, int]:
    return MismatchScanner(background).search([query], max_mismatches)[0]


def with_n(sequence: str, positions: list[int], n: str = "N") -> str:
    sequence = list(sequence)
    for p in positions:
        sequence[p] = n
    return "".join(sequence)


def test_exact_match():
    assert scan({"GENE t1": TARGET}, QUERY) == {"GENE": 0}


def test_n_in_query_is_a_match():
    assert scan({"GENE t1": TARGET}, with_n(QUERY, [0, 12, 24])) == {"GENE": 0}


def test_n_in_target_is_a_match():
    assert scan({"GENE t1": with_n(TARGET, [3, 15, 27])}, QUERY) == {"GENE": 0}


def test_n_on_both_sides_is_a_match():
    # Same position on both sides, and different positions on each side
    assert scan({"GENE t1": with_n(TARGET, [10, 20])}, with_n(QUERY, [7, 17, 2])) == {"GENE": 0}


def test_lowercase_n_is_a_match():
    assert scan({"GENE t1": with_n(TARGET.lower(), [5], "n")}, with_n(QUERY, [9], "n")) == {"GENE": 0}


def test_n_in_every_seed_position():
    # An 'N' at each position of the query, and of the target, in turn: always 0 mismatches
    for p in range(len(QUERY)):
        assert scan({"GENE t1": TARGET}, with_n(QUERY, [p])) == {"GENE": 0}
        assert scan({"GENE t1": with_n(TARGET, [p + 3])}, QUERY) == {"GENE": 0}


def test_n_with_real_mismatches_counts_only_the_real_ones():
    query = mutate(random.Random(1), QUERY, 3)
    assert scan({"GENE t1": TARGET}, with_n(query, [1, 2])) == {"GENE": count_mismatches(with_n(query, [1, 2]), QUERY)}
    assert scan({"GENE t1": with_n(TARGET, [4, 5])}, query) == {"GENE": 3}


def test_n_does_not_hide_mismatches_above_the_limit():
    query = QUERY
    target = TARGET[:3] + mutate(random.Random(2), QUERY, 5) + TARGET[28:]
    assert scan({"GENE t1": target}, query) == {}
    # One 'N' on a mismatched base brings it down to 4, which is a hit
    first_mismatch = next(i for i in range(25) if target[3 + i] != QUERY[i])
    assert scan({"GENE t1": with_n(target, [3 + first_mismatch])}, query) == {"GENE": 4}
    assert scan({"GENE t1": target}, with_n(query, [first_mismatch])) == {"GENE": 4}


def test_n_respects_a_lower_limit():
    query = mutate(random.Random(3), QUERY, 2)
    matched = next(i for i in range(25) if query[i] == QUERY[i])  # An 'N' here removes no mismatch
    target = with_n(TARGET, [3 + matched])
    assert scan({"GENE t1": target}, query, max_mismatches=1) == {}
    assert scan({"GENE t1": target}, query, max_mismatches=2) == {"GENE": 2}


def test_n_on_the_reverse_strand():
    assert scan({"GENE t1": reverse_complement(with_n(TARGET, [8]))}, with_n(QUERY, [20])) == {"GENE": 0}


def test_run_of_n_in_target_matches_any_query():
    # A run of 'N' as long as the query matches it with 0 mismatches
    assert scan({"GENE t1": "ACGT" + "N" * 30 + "ACGT"}, QUERY) == {"GENE": 0}


def test_other_iupac_codes_are_mismatches():
    assert scan({"GENE t1": with_n(TARGET, [5, 15], "R")}, QUERY) == {"GENE": 2}
    assert scan({"GENE t1": TARGET}, with_n(QUERY, [5], "Y")) == {"GENE": 1}
    # ...except where they face an 'N'
    assert scan({"GENE t1": with_n(TARGET, [8], "N")}, with_n(QUERY, [5], "Y")) == {"GENE": 0}


def test_n_does_not_join_two_sequences():
    # The query spans the end of one sequence and the start of the next: no hit
    background = {"A t1": QUERY[:12], "B t1": QUERY[12:]}
    assert scan(background, with_n(QUERY, [11, 12])) == {}


@pytest.mark.parametrize("seed", range(40))
def test_matches_brute_force(seed):
    rng = random.Random(seed)
    background = {}
    for g in range(4):
        sequence = random_sequence(rng, rng.randint(40, 120))
        if g == 0:
            # Plant near-copies of a query, so there are hits to find
            sequence = sequence[:10] + "QUERY" + sequence[10:]
        background[f"G{g} t{g}"] = sequence
    query = random_sequence(rng, 25)
    planted = mutate(rng, query, rng.randint(0, 6))
    background["G0 t0"] = background["G0 t0"].replace("QUERY", planted)
    background["G1 t1"] = background["G1 t1"] + reverse_complement(mutate(rng, query, rng.randint(0, 6)))
    # Sprinkle 'N' (and a few other non-ACGT bases) over both sides
    for name in background:
        background[name] = with_n(background[name], rng.sample(range(len(background[name])), rng.randint(0, 8)))
    if seed % 5 == 0:
        background["G2 t2"] = with_n(background["G2 t2"], list(range(20, 32)))  # A long run of 'N'
    if seed % 7 == 0:
        background["G3 t3"] = with_n(background["G3 t3"], rng.sample(range(len(background["G3 t3"])), 3), "R")
    query = with_n(query, rng.sample(range(25), rng.randint(0, 4)))
    for max_mismatches in (2, 4):
        assert scan(background, query, max_mismatches) == brute_force(background, query, max_mismatches)


def test_many_queries_in_one_search():
    rng = random.Random(99)
    background = {f"G{g} t{g}": with_n(random_sequence(rng, 300), rng.sample(range(300), 10)) for g in range(5)}
    queries = []
    for _ in range(30):
        name = rng.choice(list(background))
        start = rng.randint(0, 275)
        queries.append(with_n(mutate(rng, background[name][start:start + 25].replace("N", "A"), rng.randint(0, 5)),
                              rng.sample(range(25), rng.randint(0, 2))))
    results = MismatchScanner(background).search(queries, 4)
    assert results == [brute_force(background, q, 4) for q in queries]
