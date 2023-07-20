from performance_curve.pipeline import Pipeline
import performance_curve.calculations.conversion as cvt

import matplotlib.pyplot as plt

pl = Pipeline({
    "q_oil": 2_000,
    "mu_oil": 2,
    "rho_oil": 49.9,

    "q_gas": 1_000_000,
    "mu_gas": 0.0131,
    "rho_gas": None,
    "sg_gas": 0.709,
    "z": 0.935,

    "theta": 0,

    "sigma": 30,

    "id": cvt.in_to_feet(2.259),

    "pressure": 800,
    "temperature": cvt.fahrenheit_to_rankine(175),
})

rho_gas = pl.rho_gas
print("Gas density (lbm/ft^3):", rho_gas)

dp_dz = pl.dp_dz_friction()

print(pl.friction_factor_no_slip())
print(pl.friction_factor_two_phase())

print("Pressure gradient (friction):", dp_dz, "psi/ft")

# Using pressure gradient, we can calculate
## the pressure for each segment
segment = 100 # per 100 ft

depth_data = [i * segment for i in range(0, 21)]
pressure_data = []

for i in range(len(depth_data)):
    if (i == 0):
        pressure_data.append(pl.pressure)
    
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