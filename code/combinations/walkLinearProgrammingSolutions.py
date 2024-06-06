
from itertools import product

def branchCondition(p, N):
    a, b, c, d, e, f = p
    sumCondition = sum(p) == N
    minimumCondition = a > 0
    descendantCondition = f == 0 or c > 0
    middleCondition = (a == b+1) and (c == d)
    return sumCondition and minimumCondition and descendantCondition and middleCondition

def rootCondition(p, N):
    a, b, c, d, e, f = p
    sumCondition = sum(p) == N
    minimumCondition = a > 0
    descendantCondition = f == 0 or c > 0
    topCondition = (a == b) and (c == d)
    return sumCondition and minimumCondition and descendantCondition and topCondition

def leafCondition(p, N):
    a, b, c, d, e, f = p
    sumCondition = sum(p) == N
    minimumCondition = a > 0
    descendantCondition = f == 0 or c > 0
    bottomCondition = (a == b+1) and (c == d+1)
    return sumCondition and minimumCondition and descendantCondition and bottomCondition

conditions = {
    "root": rootCondition,
    "branch": branchCondition,
    "leaf": leafCondition
}

def sufficient(N, c="leaf"):
    values = list(range(N))
    solutions = list(product(values, repeat=6))
    return [s for s in solutions if conditions[c](s, N)]
