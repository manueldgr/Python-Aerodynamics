import sympy as sp
import numpy as np

x = sp.symbols('x')

"""Global variables"""
c = 0.8  # chord
y_ux = -4.2179 * (x ** 6) + 14.217 * (x ** 5) - 18.83 * (x ** 4) + 12.525 * (x ** 3) - 4.6582 * (
            x ** 2) + 0.958 * x + 0.0047  # Function of the upper surface
y_lx = 3.6716 * (x ** 6) - 12.228 * (x ** 5) + 15.902 * (x ** 4) - 10.189 * (x ** 3) + 3.3337 * (
                x ** 2) - 0.4815 * (
               x) - 0.0076  # Function of the lower surface
"""Slopes"""
dy_ux = sp.diff(y_ux, x)  # Slope for upper surface
dy_lx = sp.diff(y_lx, x)  # Slope for lower surface

""" Cp_ux,Cp_uy FOR 23.5 KIAS"""
cp_ux_1 = np.array([474.44*x**6-1362.5*x**5+1551.2*x**4-885.22*x**3+261.66*x**2-36.126*x+0.9424,
                    1611.5*x**6-4228.5*x**5+4286.6*x**4-2101*x**3+505.72*x**2-52.573*x+0.7802,
                    -159.37*x**6+386.55*x**5-320.86*x**4+80.713*x**3+21.976*x**2-12.493*x+0.7801])

cp_lx_1 = np.array([216.26*x**6-579.86*x**5+599.45*x**4-295.98*x**3+68.358*x**2-5.349*x-0.2893,
                    115.57*x**6-279.97*x**5+251.42*x**4-104.28*x**3+21.339*x**2-2.9081*x+0.2426,
                    -778.32*x**6+1993*x**5-1951.1*x**4+909.16*x**3-206.45*x**2+22.581*x-1.6012])

""" Cp_ux,Cp_uy FOR 43.45 KIAS"""
cp_ux_2 = np.array([766.33*x**6-2178.3*x**5+2457.6*x**4-1387.8*x**3+405.28*x**2-54.986*x+1.0673,
                    1691.1*x**6-4638.9*x**5+4943.7*x**4-2562.7*x**3+657.94*x**2-72.28*x+0.2706,
                    249.8*x**6-658.11*x**5+694.94*x**4-383.37*x**3+120.05*x**2-20.109*x+0.7235])
cp_lx_2 = np.array([-13.194*x**6+29.089*x**5-22.063*x**4+7.1925*x**3-2.2882*x**2+1.3143*x-0.3542,
                    77.432*x**6-177.19*x**5+155.47*x**4-68.718*x**3+18.012*x**2-3.4564*x+0.5511,
                    -129.06*x**6+333.31*x**5-332.88*x**4+164.17*x**3-43.967*x**2+7.4045*x-0.8805])

""" Angle Of Attack """
alpha = np.array([0, 4, -4])


def c_axial_normal_1(cp_u, cp_l, dy_u, dy_l, a):
    """Calculate el coefficient axial y normal para 23.5 KIAS"""
    #  c = 0.8
    fa = cp_u * dy_u
    fb = cp_l * dy_l

    print("RESULTADOS PARA Cn A 23.5 KIAS")
    for i, value_1 in np.ndenumerate(cp_u):
        for j, value_2 in np.ndenumerate(cp_l):
            for k, l in np.ndenumerate(a):
                if i == j == k:
                    cn = (1 / c) * sp.integrate(value_2 - value_1, (x, 0, c))
                    print(f"Cn (AoA {l}°)\t= {cn:.3f}")

    print('-' * 31)
    print("RESULTADOS PARA Ca A 23.5 KIAS")
    for ii, value_11 in np.ndenumerate(fa):
        for jj, value_22 in np.ndenumerate(fb):
            for k, l in np.ndenumerate(a):
                if ii == jj == k:
                    ca = (1 / c) * sp.integrate(value_11 - value_22, (x, 0, c))
                    print(f"Ca (AoA {l}°)\t= {ca:.3f}")
    print('\n')


def c_axial_normal_2(cp_u, cp_l, dy_u, dy_l, a):
    """Calculate el coefficient axial para 43.45 KIAS"""
    #  c = 0.8
    fa = cp_u * dy_u
    fb = cp_l * dy_l

    print("RESULTADOS PARA Cn A 43.45 KIAS")
    for i, value_1 in np.ndenumerate(cp_u):
        for j, value_2 in np.ndenumerate(cp_l):
            for k, l in np.ndenumerate(a):
                if i == j == k:
                    cn = (1 / c) * sp.integrate(value_2 - value_1, (x, 0, c))
                    print(f"Cn (AoA {l}°)\t= {cn:.3f}")

    print('-' * 31)
    print("RESULTADOS PARA Ca A 43.45 KIAS")
    for ii, value_11 in np.ndenumerate(fa):
        for jj, value_22 in np.ndenumerate(fb):
            for k, l in np.ndenumerate(a):
                if ii == jj == k:
                    ca = (1 / c) * sp.integrate((value_11 - value_22), (x, 0, c))
                    print(f"Ca (AoA {l}°)\t= {ca:.4f}")
    print('\n')


def c_moment_le(cp_u1, cp_l1, cp_u2, cp_l2, dy_u, dy_l, a):
    """ Calculate los coefficient de momento con respect al leading edge """
    #  c = 0.8

    fa = cp_u1 * dy_u
    fb = cp_l1 * dy_l

    print("RESULTADOS PARA Cm-LE A 23.5 KIAS")
    for i, value_1 in np.ndenumerate(cp_u1):
        for j, value_2 in np.ndenumerate(cp_l1):
            for k, value_3 in np.ndenumerate(fa):
                for q, value_4 in np.ndenumerate(fb):
                    for m, value_5 in np.ndenumerate(a):
                        if i == j == k == q == m:
                            aa = sp.integrate((value_1 - value_2) * x, (x, 0, c))
                            bb = sp.integrate(value_3 * y_ux, (x, 0, c))
                            cc = sp.integrate((-value_4) * y_lx, (x, 0, c))
                            cm1 = (1 / c ** 2) * (aa + bb + cc)
                            print(f"Cm (AoA {value_5}°)\t= {cm1:.3f}")

    print('-' * 31)
    fc = cp_u2 * dy_ux
    fd = cp_l2 * dy_lx

    print("RESULTADOS PARA Cm-LE A 43.45 KIAS")
    for n, value_6 in np.ndenumerate(cp_u2):
        for o, value_7 in np.ndenumerate(cp_l2):
            for p, value_8 in np.ndenumerate(fc):
                for q, value_9 in np.ndenumerate(fd):
                    for r, value_10 in np.ndenumerate(a):
                        if n == o == p == q == r:
                            dd = sp.integrate((value_6 - value_7) * x, (x, 0, c))
                            ee = sp.integrate(value_8 * y_ux, (x, 0, c))
                            ff = sp.integrate((-value_9) * y_lx, (x, 0, c))
                            cm2 = (1 / c ** 2) * (dd + ee + ff)
                            print(f"Cm (AoA {value_10}°)\t= {cm2:.3f}")


c_axial_normal_1(cp_ux_1, cp_lx_1, dy_ux, dy_lx, a=alpha)
c_axial_normal_2(cp_ux_2, cp_lx_2, dy_ux, dy_lx, a=alpha)
c_moment_le(cp_ux_1, cp_lx_1, cp_ux_2, cp_lx_2, dy_ux, dy_lx, a=alpha)
