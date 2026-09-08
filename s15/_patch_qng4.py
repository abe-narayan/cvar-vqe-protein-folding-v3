import io

p = 's15/qgeom_qng.py'
s = io.open(p, encoding='utf-8').read()

old = '''LRS = {"sgd": (0.03, 0.1, 0.3, 1.0, 3.0),
       "adam": (0.02, 0.05, 0.1, 0.2, 0.4),
       "qng": (0.03, 0.1, 0.3, 1.0),
       "qng_diag": (0.03, 0.1, 0.3, 1.0),
       "qng_adam": (0.02, 0.05, 0.1, 0.2),
       "ng_shots": (0.03, 0.1, 0.3, 1.0)}
LAMS = {"qng": (1e-4, 1e-3, 1e-2, 1e-1), "qng_diag": (1e-3,), "qng_adam": (1e-3, 1e-2),
        "ng_shots": (1e-2, 1e-1)}'''
new = '''# The BASELINES keep the widest sweep -- a negative result about QNG must never rest on an
# under-tuned SGD or Adam.  The QNG family's grid was trimmed for throughput on a saturated
# box; the retained values bracket the optimum found in the full sweep of the first run
# (qng chose lr in {0.1, 0.3} and lam in {1e-3, 1e-2} in every completed cell).
LRS = {"sgd": (0.03, 0.1, 0.3, 1.0, 3.0),
       "adam": (0.02, 0.05, 0.1, 0.2, 0.4),
       "qng": (0.03, 0.1, 0.3),
       "qng_diag": (0.03, 0.1, 0.3),
       "qng_adam": (0.05, 0.2),
       "ng_shots": (0.03, 0.1, 0.3)}
LAMS = {"qng": (1e-3, 1e-2), "qng_diag": (1e-3,), "qng_adam": (1e-2,),
        "ng_shots": (1e-1,)}'''
assert old in s
s = s.replace(old, new)

old2 = '''LEAN = [("ring", 1), ("all_to_all", 2), ("ring", 2), ("block", 2), ("chain", 2),
        ("chain", 3)]


def main(lean=True):
    G.wait_mem(0.8, "qgeom_qng")
    G.ck_load(TAG)
    t0 = time.time()
    lad = LEAN if lean else None
    sd = (0, 1) if lean else (0, 1, 2)'''
new2 = '''# Four rungs spanning the measured conditioning range: 1.0 / 9.6 / 86 / 2004.
LEAN = [("ring", 1), ("ring", 2), ("chain", 2), ("chain", 3)]


def main(lean=True, iters=90):
    G.wait_mem(0.8, "qgeom_qng")
    G.ck_load(TAG)
    t0 = time.time()
    lad = LEAN if lean else None
    sd = (0, 1) if lean else (0, 1, 2)'''
assert old2 in s
s = s.replace(old2, new2)

s = s.replace('''    vqe_experiment(seeds=sd, ladder=lad)
    vqe_equal_cost(seeds=sd, ladder=lad)''',
'''    vqe_experiment(seeds=sd, ladder=lad, iters=iters)
    vqe_equal_cost(seeds=sd, ladder=lad)''')

io.open(p, 'w', encoding='utf-8', newline='\n').write(s)
print('patched')
