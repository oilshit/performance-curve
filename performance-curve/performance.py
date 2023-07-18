from math import pow, log10, log, pi
import matplotlib.pyplot as plt

import calculations.conversion as cvt

BARREL_CONV = 5.615
DAY_TO_SECONDS = 86400

STANDARD_PRESSURE = 14.7
STANDARD_TEMPERATURE = cvt.fahrenheit_to_rankine(60)

def first_graph(n_L):
    """
    Relation of graph between c_NL and n_L
    """

    X = log10(n_L + 3)

    Y1 = -2.6951 + 0.15841 * X - 0.551 * pow(X, 2)
    Y2 = 0.54785 * pow(X, 3) - 0.12195 * pow(X, 4)
    Y = Y1 + Y2

    c_NL = pow(10, Y)

    return c_NL


def second_graph(n_vl, n_vg, p, c_NL, n_D):
    A = n_vl / pow(n_vg, 0.575)
    B = pow(p / STANDARD_PRESSURE, 0.1)
    C = c_NL / n_D

    X = A * B * C

    log_X = log10(X) + 6

    Y1 = -0.10307 + 0.61777 * log_X
    Y2 = -0.63295 * pow(log_X, 2) + 0.29598 * pow(log_X, 3)
    Y3 = -0.0401 * pow(log_X, 4)

    Y = Y1 + Y2 + Y3

    yl_per_phi = Y

    return yl_per_phi

def third_graph(n_vg, n_L, n_D):
    """
    Relation of third graph
    phi versus A function
    """

    A = n_vg * pow(n_L, 0.38) / pow(n_D, 2.14)

    if (A <= 0.01):
        phi = 1
    
    else:
        phi1 = 0.91163 - 4.82176 * A + 1_232.25 * pow(A, 2)
        phi2 = -22_253.6 * pow(A, 3) + 116_174.3 * pow(A, 4)

        phi = phi1 = phi2
    
    return phi

def gas_density(p, sg_g, z, T):
    rho_g = 28.97 * sg_g * p / (z * 10.73 * T)

    return rho_g

def fanning_friction_factor(n_Re, roughness):
    A = roughness / 3.7605
    B = -5.0452 / n_Re
    
    C = pow(roughness, 1.1098) / 2.8257
    D = pow(7.149 / n_Re, 0.8981)
    E = log10(C + D)

    F = A + B * E

    f_F = pow(1 / (4 * log10(F)), 2)

    return f_F

# initial variable
## flow rate properties
q_oil = 2000                             # in bbl/d
q_gas = 1_000_000                        # in stbd

## tubing properties
id = cvt.in_to_feet(2.259)               # inner diameter, in feet
roughness = 0.0006                       # tubing relative roughness

## pressure and temperature properties
p_wellhead = 800                         # in psia
t = cvt.fahrenheit_to_rankine(175)       # in Rankine degree

## oil properties
rho_o = cvt.gcc_to_lbm_per_ft3(0.8)      # oil density, in lbm/ft^3
mu_o = 2                                 # oil viscosity, in cp
sigma_o = 30                             # interfacial tension, dynes/cm

## gas properties
sg_gas = 0.709                           # gas specific gravity
z = 0.935                                # gas compressibility factor
mu_g = 0.0131                            # gas viscosity, in cp

# AREA CALCULATION
# ====================================
section_area = 0.25 * pi * pow(id, 2)
print("Section Area (A):", section_area, "\n")
# ====================================

# SUPERFICIAL VELOCITY
# ====================================
us_liquid = q_oil / section_area * (
    BARREL_CONV / DAY_TO_SECONDS
)
print("Liquid superficial velocity (ft/s):", us_liquid)

us_gas = q_gas / section_area * (
    z * t / STANDARD_TEMPERATURE * STANDARD_PRESSURE / p_wellhead
        / DAY_TO_SECONDS
)
print("Gas superficial velocity (ft/s):", us_gas)

u_mixture = us_gas + us_liquid
print("Mixture velocity (ft/s):", u_mixture)
print()
# ====================================

# CHECKING OF FLOW TYPE
# ====================================
gas_input_griction = us_gas / u_mixture
print("Gas input friction:", gas_input_griction)

l_b = 1.071 - 0.2218 * (pow(u_mixture, 2) / id)
print("L_B (based on modified Hagedorn-Brown):", l_b)

if (l_b < 0.13):
    print("Since the L_B constant is less than 0.13, then L_B is considered as 0.13")
    l_b = 0.13

if (gas_input_griction < l_b):
    print("Bubble flow!")
else:
    print("Not a bubble flow (since gas friction > L_B constant)!")

print()

# ====================================

# CALCULATION OF MODIFIED HAGEDORN-BROWN CONSTANT
# ====================================
n_vl = 1.938 * us_liquid * pow(rho_o / sigma_o, 0.25)
n_vg = 1.938 * us_gas * pow(rho_o / sigma_o, 0.25)
n_D = 120.872 * id * pow(rho_o / sigma_o, 0.5)
n_L = 0.1572 * mu_o * pow(1 / (rho_o * pow(sigma_o, 3)), 0.25)

print("n_VL, n_VG, n_D, n_L =", n_vl, n_vg, n_D, n_L)
# ====================================

c_NL = first_graph(n_L)
print("According to first graph, c_NL =", c_NL)

yl_per_phi = second_graph(
    n_vl=n_vl,
    n_vg=n_vg,
    p=p_wellhead,
    c_NL=c_NL,
    n_D=n_D
)
print("According to second graph, yL per phi =", yl_per_phi)

phi = third_graph(
    n_vg=n_vg,
    n_L=n_L,
    n_D=n_D
)
print("According to third graph, phi =", phi)

y_L = phi * yl_per_phi
print("So that, y_L =", y_L)
print()

# MIXTURE DENSITY
# ====================================
rho_g = gas_density(p_wellhead, sg_gas, z, t)
print("Gas density (lbm/ft^3) =", rho_g)

mt = section_area * (us_liquid * rho_o + us_gas * rho_g) * DAY_TO_SECONDS
print("Mt (in lbm/day) =", mt)
print()
# ====================================

# REYNOLDS NUMBER
# ====================================
n_Re = 2.2E-2 * mt / (id * pow(mu_o, y_L) * pow(mu_g, (1 - y_L)))
print("Reynolds Number =", n_Re)
# ====================================

# FRICTION FACTOR
# ====================================
f_f = fanning_friction_factor(n_Re, roughness)
print("Fanning Friction Factor =", f_f)
print()
# ====================================

# INSITU AVERAGE DENSITY
# ====================================
rho_avg = y_L * rho_o + (1 - y_L) * rho_g
print("Insitu average density (lbm/ft^3) =", rho_avg)
print()
# ====================================

# INSITU AVERAGE DENSITY
# ====================================
dp_dz = (1 / 144) * (
    rho_avg + (
        f_f * pow(mt, 2) / (7.143E10 * pow(id, 5) * rho_avg)
    )
)

print("Pressure gradient at top of tubing (wellhead, psi/ft) =", dp_dz)
# ====================================

# Using pressure gradient, we can calculate
## the pressure for each segment
segment = 100 # per 100 ft

depth_data = [i * segment for i in range(0, 21)]
pressure_data = []

for i in range(len(depth_data)):
    if (i == 0):
        pressure_data.append(p_wellhead)
    
    else:
        pressure_data.append(
            pressure_data[i -1] + (dp_dz) * segment
        )

plt.plot(pressure_data, depth_data, label="Pressure gradient")

plt.ylabel("Depth (ft)")
plt.text((pressure_data[len(pressure_data) // 3]), -250, 'Pressure (psia)')

plt.subplots_adjust(top=0.85)

plt.ylim(depth_data[-1] + 300, 0)
plt.tick_params(axis='both', which='major',
               labelsize=10, labelbottom=False,
               bottom=False, top=True, labeltop=True)

plt.show()