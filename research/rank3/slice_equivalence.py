"""Exact cyclic mutation-loop signatures for the sixteen representatives.

Canonical labels follow the terminal permutation orbits, rooted at the three
mutation events. Exhaust cyclic rotation and adjacent commuting mutations.
Also negate the exchange matrix for exchange of plus and minus.
"""
from __future__ import annotations

from collections import deque
import hashlib
import json

from classify_rank3 import identity, mutate
from solve_lifts import HERE


def connected_component(b, start):
    seen = {start}
    while True:
        more = seen | {j for i in seen for j, x in enumerate(b[i]) if x}
        if more == seen:
            return seen
        seen = more


def slice_loop(cert):
    b, p, mutation = cert["B"], cert["relabel_old_to_new"], cert["mutation_vertices"]
    component = connected_component(b, 0)
    current = set(component)
    t = 0
    while True:
        t += 1
        current = {p[i] for i in current}
        if current == component:
            break
        assert not current & component
    assert t * len(component) == len(b)
    order = sorted(component)
    index = {v: i for i, v in enumerate(order)}
    forward = list(range(len(b)))
    word = []
    for k in range(t):
        word.extend(index[v] for v in order if forward[v] in mutation)
        forward = [p[v] for v in forward]
    permutation = tuple(index[forward[v]] for v in order)
    base = tuple(tuple(b[i][j] for j in order) for i in order)
    assert len(word) == 3
    return t, base, tuple(word), permutation


def canonical_state(b, word, permutation):
    order = []
    lengths = []
    for root in word:
        cycle = [root]
        v = permutation[root]
        while v != root:
            cycle.append(v)
            v = permutation[v]
        lengths.append(len(cycle))
        order.extend(cycle)
    assert sorted(order) == list(range(len(b)))
    index = {v: i for i, v in enumerate(order)}
    cb = tuple(tuple(b[i][j] for j in order) for i in order)
    cw = tuple(index[v] for v in word)
    cp = tuple(index[permutation[v]] for v in order)
    signature = tuple(lengths)+tuple(x for row in cb for x in row)
    return signature, (cb, cw, cp)


def mutation_matrix(b, k):
    return tuple(tuple(row) for row in mutate(b, identity(len(b)), k)[0])


def equivalence_states(b, word, permutation):
    pending = deque([(b, word, permutation)])
    seen = set()
    while pending:
        b, word, permutation = pending.popleft()
        signature, state = canonical_state(b, word, permutation)
        if signature in seen:
            continue
        seen.add(signature)
        b, word, permutation = state
        inv = tuple(permutation.index(i) for i in range(len(b)))
        pending.append((mutation_matrix(b, word[0]), word[1:] + (inv[word[0]],), permutation))
        prefix = b
        for k in range(2):
            if prefix[word[k]][word[k+1]] == 0:
                other = list(word)
                other[k], other[k+1] = other[k+1], other[k]
                pending.append((b, tuple(other), permutation))
            prefix = mutation_matrix(prefix, word[k])
        if len(seen) > 20000:
            raise RuntimeError("Unexpectedly large equivalence orbit")
    return seen


def main():
    witnesses = [json.loads(line) for line in (HERE / "lift_feasibility.jsonl").read_text().splitlines()]
    results, minima = [], {}
    for item in witnesses:
        if item["status"] != "sat":
            continue
        t, b, word, p = slice_loop(item["certificate"])
        states = equivalence_states(b, word, p)
        opposite_b = tuple(tuple(-x for x in row) for row in b)
        states |= equivalence_states(opposite_b, word, p)
        minimum = min(states)
        assert minimum not in minima, (item["id"], minima.get(minimum))
        minima[minimum] = item["id"]
        entry = {"id": item["id"], "components": t, "slice_vertices": len(b),
                 "slice_B": b, "slice_mutation_word": word, "slice_relabel_old_to_new": p,
                 "equivalent_canonical_states_including_sign_exchange": len(states),
                 "canonical_signature": minimum,
                 "signature_sha256": hashlib.sha256(repr(minimum).encode()).hexdigest()}
        results.append(entry)
        print(item["id"], "components", t, "vertices", len(b), "states", len(states), flush=True)
    (HERE / "slice_signatures.json").write_text(json.dumps(results, indent=2)+"\n")
    print("All sixteen cyclic mutation-loop equivalence classes are distinct.")


if __name__ == "__main__":
    main()
