import sympy as sm
from functools import reduce


def test0():
    x1, x2, x3, x4 = sm.var("x1 x2 x3 x4")
    m2 = (x1+x2+x3+x4)**3
    coefs = m2.expand().as_coefficients_dict()
    print("coefs:")
    print(coefs)
    res = binom_factors(3, 4)
    print("binom results:")
    print(res)
    return(coefs)


def binom_factors(n, r):
    def func(state):
        fact = sm.factorial(n)/reduce(
            lambda acc, x: acc*sm.factorial(x), state, 1)
        return(fact)
    results = _binom_factors([[0 for i in range(r)]], [], n, r)
    fresults = [func(state) for state in results]
        
    return(list(zip(results, fresults)))
    

def _binom_factors(states, results, n, r):
    # finish:
    if len(states) == 0:
        return(results)

    state = states.pop(0)
    for i in range(r):
        new_state = state[:i]+[state[i]+1] + state[i+1:]
        if sum(new_state) == n:
            if new_state not in results:
                results.append(new_state)
        elif(sum(new_state) < n):
            states.append(new_state)

    return(_binom_factors(states, results, n, r))


if __name__ == "__main__":
    test0()
    # res = binom_factors(3, 4)
    # print("res:")
    # print(res)
    # print("len(res):", len(res))
