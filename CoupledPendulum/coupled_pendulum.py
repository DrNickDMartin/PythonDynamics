from scipy.integrate import solve_ivp


sol = solve_ivp(fun=lambda t, y: cp_pen(t, y, l, a, k), [0, 10], [2, 4, 8])
print(sol.t)
print(sol.y)

def cp_pen(t, y, l, a, k):
