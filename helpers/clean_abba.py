#!/usr/bin/env python3
import sys

filename=sys.argv[1]
print(filename)
with open(filename, 'rt') as fh:
    datasets = set()
    for ln in fh:
        toks = ln.rstrip().split(',')
        for dataset in toks[4:]:
            datasets.add(dataset)


names = list(sorted(datasets))


with open(filename, 'rt') as fh:
    print(','.join(['perm', 'step', 'delta', 'deltadelta', 'added'] + names))
    last_delta = 0
    last_perm = -1
    step = 0
    for ln in fh:
        toks = ln.rstrip().split(',')
        if toks[0] == 'ngen':
            continue
        k, delta, perm = int(toks[1]), float(toks[2]), int(toks[3])
        if perm != last_perm:
            last_delta = 0
            step = 0
        datasets = toks[4:-1]
        added_dataset = toks[-1]
        dataset_vector = []
        for name in names:
            if name == added_dataset:
                dataset_vector.append(0)
            elif name not in datasets:
                dataset_vector.append('-1')
            else:
                assert name in datasets
                dataset_vector.append(1)
        assert delta > last_delta
        line = [perm, step, delta, delta - last_delta, added_dataset] + dataset_vector
        print(','.join(map(str, line)))
        last_delta = delta
        last_perm = perm
        step += 1
