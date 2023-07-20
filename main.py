from performance_curve.pipeline import Pipeline
import performance_curve.calculations.conversion as cvt

pl = Pipeline({
    "q_oil": 2_000,
    "mu_oil": 2,
    "rho_oil": 49.9,

    "q_gas": 1_000_000,
    "mu_gas": 0.0131,
    "sg_gas": 0.709,
    "z": 0.935,

    "id": cvt.in_to_feet(2.259),

    "pressure": 800,
    "temperature": cvt.fahrenheit_to_rankine(175),
})

print(pl.sg_gas, pl.rho_gas)
print(pl.superficial_velocity("liquid"))
print(pl.superficial_velocity("gas"))

u_sl = pl.superficial_velocity("liquid")
u_sg = pl.superficial_velocity("gas")
u_m = pl.mixture_velocity(u_sl, u_sg)
print(u_m)

lambda_l, lambda_g = pl.input_fraction()
print(lambda_l, lambda_g)

holdup = pl.holdup_constant()
print(holdup["type"])
print(holdup["constants"])