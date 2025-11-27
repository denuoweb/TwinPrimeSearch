# Hotspot-Guided Search for a 1000-Digit Twin Prime Pair

This repository contains a parallel, hotspot-guided search for large twin primes (and prime quadruplets) based on residue-class sieving and Miller–Rabin testing. The current implementation has produced a 1000-digit twin prime pair inside a “4-prime decade” of the form

\[(n+1,\ n+3,\ n+7,\ n+9),\]

where all four entries are prime and thus give two twin prime pairs and a prime constellation of type (2,4,2).

The main search script is:

```text
parallel_hotspot_twinprime_search.py
```

Below is a precise description of the algorithm, the resulting 1000-digit twin prime pair, and instructions for reproducing and verifying the computation.

---

## 1. The 1000-Digit Twin Prime Pair

The script found a 4-prime decade in which

- q1 = n + 1 and q2 = n + 3 form a twin prime pair, and  
- `digits(q1) = digits(q2) = 1000`.

In particular, at `decade_offset = 392` along one of the searched progressions, the log contains:

```text
TWIN_PAIR (q1,q2) at decade_offset=392: (p, q), digits=1000
```

with

- p = q1  
- q = q2 = p + 2.

### 1.1. Compact description

Let

\[
p = 10^{999} + a, \quad q = p + 2,
\]

where

```text
a = 219191175840524707078943676751
```

Thus

```text
p = 10^999 + 219191175840524707078943676751
q = 10^999 + 219191175840524707078943676753
```

Both p and q have exactly 1000 decimal digits and differ by 2.

### 1.2. Full decimal expansions

For completeness, here are the full 1000-digit values. Each number is given as a single line of decimal digits.

<details>
<summary>Click to expand the full 1000-digit twin prime pair</summary>

```text
p =
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000219191175840524707078943676751

q =
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000219191175840524707078943676753
```

</details>

### 1.3. SHA-256 checksums

To make transcription errors easy to detect, the SHA-256 hashes of the **raw decimal strings** (no whitespace) are:

```text
SHA-256(p) = 8290668d7bff63e3a18bad2ca89bd01c02c73a05a0847e999400c1baea906d88
SHA-256(q) = 8d40de659a649bc0820b8498af2e9bb69a0976468a7aaa699c8c033f3a90abff
```

Anyone can independently recompute these hashes from the decimal expansions above to confirm that they are working with exactly the same integers.

---

## 2. Algorithm and Implementation

The core idea is to exploit the structure of “4-prime decades”

\[[n, n+9]\]

with candidate positions

\[
q_1 = n+1,\quad q_2 = n+3,\quad q_3 = n+7,\quad q_4 = n+9.
\]

If all four are prime, we obtain:

- two twin prime pairs: (q1, q2) and (q3, q4),
- a prime constellation of shape (2, 4, 2).

Instead of testing arbitrary large integers, the script pushes as much work as possible into:

1. modular constraints modulo small primes, and  
2. a precomputed progression via the CRT,

and then performs probabilistic primality tests only on candidates that survive this filtering.

### 2.1. Small-prime sieving and residue patterns

Given a bound B (default `--small-prime-max 71`), the script computes all small primes

```python
P = small_primes_up_to(B)
```

For each small prime p in P, it selects a residue r_p in {0,1,...,p-1}. The goal is to find residue patterns

```python
rpat: Dict[int, int]  # p -> r_p
```

such that none of the four target positions is divisible by any small prime:

\[
(r_p + 1) \not\equiv 0,\quad
(r_p + 3) \not\equiv 0,\quad
(r_p + 7) \not\equiv 0,\quad
(r_p + 9) \not\equiv 0 \pmod{p}
\]

for all p in P. This is enforced by

```python
def passes_prime_slot_constraints(rpat: Dict[int, int]) -> bool:
    for p, r in rpat.items():
        for off in PRIME_OFFSETS:  # [1, 3, 7, 9]
            if (r + off) % p == 0:
                return False
    return True
```

Only residue patterns that satisfy this constraint are considered further.

### 2.2. Hotspot scoring

Around each candidate decade `[n, n+9]`, the script considers a window of offsets

\[[n-H, \dots, n+9+H]\]

with `H = --window-H` (default `10`). For each offset m in this window that is **not** one of the prime-positions {1,3,7,9}, the script counts how often that offset is divisible by the small primes under the given residue pattern.

Intuitively, residue patterns that force many neighbors to be divisible by small primes create “hotspots” where the composite numbers are concentrated in the neighborhood, while the four target positions stay as “holes” in that sieve.

The hotspot score is:

```python
def hotspot_score(rpat: Dict[int, int],
                  composite_offsets: List[int],
                  weight_func=None) -> int:
    if weight_func is None:
        weight_func = lambda p: 1
    score = 0
    for m in composite_offsets:
        for p, r in rpat.items():
            if (r + m) % p == 0:
                score += weight_func(p)
    return score
```

Out of many randomly sampled residue patterns, the script retains only the top `--top-k-patterns` according to this score.

### 2.3. CRT and arithmetic progression

Once a residue pattern is fixed, the script uses the Chinese remainder theorem to construct a congruence

\[
n \equiv a \pmod{M}, \quad M = \prod_{p \le B} p,
\]

such that the chosen residue pattern holds for all small primes simultaneously. This is done by

```python
base_residue, modulus = crt_from_residues(rpat)
```

Given a target digit size `--digits = D` (here `1000`), the script chooses the smallest n0 in this arithmetic progression with n0 having at least D-1 digits (so that n0 + 9 will have D digits), then scans decades

\[
n = n_0 + tM,\quad t = 0,1,2,\dots
\]

up to `--max-decades-per-pattern`.

### 2.4. Primality testing and parallelism

For each decade starting at n, the script forms:

```python
q1 = n + 1
q2 = n + 3
q3 = n + 7
q4 = n + 9
```

Each q_i is tested with a standard Miller–Rabin probable-prime test, including trial division by a fixed small set of primes:

```python
def is_probable_prime(n: int, rounds: int = 16) -> bool:
    if n < 2:
        return False

    # Trial division by small primes
    for p in _SMALL_PRIMES_TRIAL:
        if n == p:
            return True
        if n % p == 0:
            return False

    # n-1 = 2^s * d
    ...
    # standard Miller–Rabin loop with random bases
```

Decades that yield twin pairs (or full quadruplets) are written to per-pattern log files. Patterns are processed in parallel via `ProcessPoolExecutor`, with each worker responsible for scanning one residue pattern:

```python
with ProcessPoolExecutor(max_workers=workers) as ex:
    ...
    ex.submit(search_along_pattern, score, rpat, digits, max_decades,
              rank, log_prefix, H)
```

The 1000-digit twin prime pair above arose as `(q1, q2)` in one such pattern at `decade_offset = 392`.

---

## 3. Running the Search

The script has a CLI interface based on `argparse`. The key options are:

```text
--digits                  Target number of digits for q1..q4
--small-prime-max         Upper bound B for the small primes used in residue patterns
--num-pattern-samples     Number of random residue patterns to sample
--top-k-patterns          Number of best hotspot patterns to retain and scan
--window-H                Half-width H of the composite-neighborhood window
--max-decades-per-pattern Max decades to scan along each CRT progression
--log-prefix              Prefix for per-pattern log filenames
--workers                 Number of worker processes (default: os.cpu_count())
```

A typical run for 1000-digit candidates might look like:

```bash
python parallel_hotspot_twinprime_search.py \
    --digits 1000 \
    --small-prime-max 71 \
    --num-pattern-samples 20000 \
    --top-k-patterns 10 \
    --max-decades-per-pattern 2000 \
    --window-H 10 \
    --log-prefix search_1000d \
    --workers 22
```

Each retained pattern will produce a log of the form:

```text
search_1000d_pattern1.log
search_1000d_pattern2.log
...
```

These logs contain lines like:

```text
PRIME q1 at decade_offset=...
TWIN_PAIR (q1,q2) at decade_offset=...
QUADRUPLET (2,4,2) at decade_offset=...
```

The 1000-digit twin pair described in Section 1 comes from one such `TWIN_PAIR` line with `decade_offset=392`.

---

## 4. Verifying the 1000-Digit Twin Prime Pair

For independent verification, the following steps are recommended:

1. **Digit count and spacing**  
   Confirm that p and q each have 1000 digits and that q - p = 2.

2. **Probable-prime checks**  
   Run a Miller–Rabin test with a sufficiently strong base set, or use a higher-level library such as `sympy`.

3. **Rigorous primality proof** (if desired)  
   For a full proof, run an ECPP implementation (e.g. Primo, Pari/GP, Sage) to produce certificates for both p and q.

The following Python snippet (using `sympy`) verifies the basic arithmetical properties and performs probable-prime tests:

```python
from sympy import isprime

p_str = "1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000219191175840524707078943676751"
q_str = "1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000219191175840524707078943676753"

p = int(p_str)
q = int(q_str)

assert len(p_str) == 1000
assert len(q_str) == 1000
assert q - p == 2

print("p and q are 1000-digit integers with q - p = 2.")
print("isprime(p) =", isprime(p))
print("isprime(q) =", isprime(q))
```

To confirm the hashes:

```python
import hashlib

def sha256_hex(s: str) -> str:
    return hashlib.sha256(s.encode("ascii")).hexdigest()

print("SHA-256(p) =", sha256_hex(p_str))
print("SHA-256(q) =", sha256_hex(q_str))
```

The outputs should match the values listed in Section 1.3.

---

## 5. Log Format and Search Statistics

Each worker pattern writes to its own log file:

```text
<log-prefix>_pattern<rank>.log
```

Within each log, lines of interest include:

- Prime hits at individual positions:

  ```text
  PRIME q1 at decade_offset=J: n=..., q1=..., digits=...
  PRIME q2 at decade_offset=J: ...
  ...
  ```

- Twin prime pairs:

  ```text
  TWIN_PAIR (q1,q2) at decade_offset=J: (p, q), digits=...
  TWIN_PAIR (q3,q4) at decade_offset=J: (p', q'), digits=...
  ```

- Full prime quadruplets:

  ```text
  QUADRUPLET (2,4,2) at decade_offset=J: (q1, q2, q3, q4), digits=...
  ```

At the end of each pattern’s search, a summary line aggregates statistics:

```text
# finished pattern: total_decades=..., total_candidates=...,   primes=..., twin_pairs=..., quadruplets=...
```

The main driver prints a global summary aggregating over all patterns:

```text
[summary] all patterns done.
  total candidates tested = ...
  total primes found      = ...
  total twin pairs        = ...
  total quadruplets       = ...
```

This makes it straightforward to track hit rates, compare different parameter sets, and identify which patterns produced interesting large-prime configurations.

---

## 6. Possible Extensions

Some natural extensions and experiments supported by this framework:

- Increasing `--digits` to push beyond 1000 digits while keeping the same structural approach.
- Modifying `PRIME_OFFSETS` to search for other constellations (e.g. prime triplets or longer k-tuples).
- Adjusting the hotspot scoring function to weight small primes differently or to prioritize different neighborhood structures.
- Integrating a deterministic primality prover (ECPP) into the post-processing pipeline to automatically certify interesting hits.

The current script already separates the fast probabilistic search (Miller–Rabin) from the final rigorous certification step, so integrating external provers can be done cleanly by post-processing the `TWIN_PAIR` and `QUADRUPLET` lines from the logs.
