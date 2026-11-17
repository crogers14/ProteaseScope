import numpy as np
i64 = np.array([1, 2, 0, 3], dtype=np.int64)
x = np.arange(10.0)
tests = {
    "repeat(int64 counts)": lambda: np.repeat(np.arange(4), i64),
    "bincount(int64)": lambda: np.bincount(i64, minlength=5),
    "bincount(int64, minlength=np.int64)": lambda: np.bincount(i64.astype(np.intp), minlength=np.int64(5)),
    "reduceat(int64 idx)": lambda: np.add.reduceat(x, np.array([0, 2, 5], dtype=np.int64)),
    "add.at(int64 idx)": lambda: np.add.at(np.zeros(5), i64, 1),
    "minimum.at": lambda: np.minimum.at(np.zeros(5), i64, 1),
    "bitwise_or.at": lambda: np.bitwise_or.at(np.zeros(5, np.int64), i64, 1),
    "fancy index int64": lambda: x[i64],
    "fancy assign int64": lambda: x.__setitem__(i64, 1.0),
    "searchsorted": lambda: np.searchsorted(x, i64),
    "lexsort": lambda: np.lexsort((i64, i64)),
    "argsort int64": lambda: np.argsort(i64),
    "isin": lambda: np.isin(i64, i64),
    "unique": lambda: np.unique(i64),
    "cumsum/where": lambda: np.where(i64 > 0, i64, 0).cumsum(),
    "take": lambda: np.take(x, i64),
    "pandas factorize": lambda: __import__("pandas").factorize(["a", "b"]),
}
print("numpy", np.__version__, "intp =", np.dtype(np.intp))
for k, f in tests.items():
    try:
        f(); print("ok  ", k)
    except Exception as e:
        print("FAIL", k, type(e).__name__, e)
