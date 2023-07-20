from typing import Union, Dict, List, Tuple
from math import log10, log, pi, pow, sin, exp

# import calculations.conversion as cvt

DAYS_TO_SECONDS = 86400
BARREL_TO_CUBIC_FEET = 5.615

# STANDARD_TEMPERATURE = cvt.fahrenheit_to_rankine(60)
STANDARD_TEMPERATURE = (60 + 460)
STANDARD_PRESSURE = 14.7
GRAVITY = 32.17

numeric = Union[Union[int, float], None] 

class Pipeline:
    def __init__(self, props: Dict):
        """
        constructor of pipeline properties
        """

        self.q_oil = props["q_oil"]                     # oil flow rate (bpd)
        self.rho_oil = props["rho_oil"]                 # oil density (lbm/ft^3)
        self.mu_oil = props["mu_oil"]                   # oil viscosity (cp)

        self.sigma = props["sigma"]                   # oil viscosity (cp)

        self.q_gas = props["q_gas"]                     # gas flow rate (SCFD)
        self.rho_gas = props["rho_gas"]                 # gas density (lbm/ft^3)
        self.mu_gas = props["mu_gas"]                   # gas viscosity (lbm/ft^3)
        self.sg_gas = props["sg_gas"]                   # gas specific gravity
        self.z = props["z"]                             # gas compressibiliy factor
        
        self.id = props["id"]                           # pipe inner diameter (ft)
        self.theta = props["theta"] or 0                # elevation angle (default horizontal = 0)

        self.pressure = props["pressure"]               # pressure (psia)
        self.temperature = props["temperature"]         # temperatrue (oR)

        if (self.rho_gas == None):
            self.rho_gas = self.gas_density(
                self.sg_gas,
                self.pressure,
                self.z,
                self.temperature,
            )

    def gas_density(
        self,
        sg_gas: numeric,
        P: numeric,
        z: numeric,
        T: numeric
    ) -> numeric:
        """
        get gas density in lbm/ft^3

        input:
            sg_gas (gas specific gravity)  : numeric
            p (pressure, psia)             : numeric
            z (gas compressibility factor) : numeric
            T (temperature, oR)            : numeric

        output:
            rho_gas (gas density, lbm/ft^3): numeric
        """
        
        try:
            rho_gas = 28.97 * sg_gas * P / (10.73 * z * T)
            return rho_gas
        
        except:
            print("Invalid density properties")
            return None
        
    def section_area(self, id: numeric) -> numeric:
        """
        get pipe section area in ft^2

        input:
            id (pipe inner diameter, ft)    : numeric

        output:
            area (ft^2)                     : numeric
        """
        area = pi / 4 * pow(id, 2)
        return area

    def superficial_velocity(self, fluid: str) -> numeric:
        """
        get superficial velocity in ft/s

        input:
            fluid (fluid type): str
        
        output:
            us_l or us_g      : numeric
        """
        area = self.section_area(self.id)

        if (fluid == "liquid"):
            unit_conversion = BARREL_TO_CUBIC_FEET / DAYS_TO_SECONDS

            us_l = self.q_oil * unit_conversion / area
            return us_l
        
        elif (fluid == "gas"):
            t_ratio = self.temperature / STANDARD_TEMPERATURE
            p_ratio = STANDARD_PRESSURE / self.pressure

            unit_conversion = 1 / DAYS_TO_SECONDS

            us_g = self.q_gas * t_ratio * p_ratio * self.z * unit_conversion / area
            return us_g
        
    def mixture_velocity(self, u_sl: numeric, u_sg: numeric) -> numeric:
        """
        get mixture superficial velocity in ft/s

        input:
            u_sl (liquid sup. velocity, ft/s): numeric
            u_sg (gas sup.velocity, ft/s)    : numeric
        
        output
            u_m (mixture sup. velocity, ft/s): numeric
        """
        u_m = u_sg + u_sl

        return u_m
    
    def mixture_density(self) -> numeric:
        """
        get mixture density of liquid and gas in lbm/ft^3
        
        input: None

        output:
            rho_m: numeric
        """

        lambda_l, lambda_g = self.input_fraction()
        rho_m = self.rho_oil * lambda_l + self.rho_gas * lambda_g

        return rho_m
    
    def mixture_viscosity(self) -> numeric:
        """
        get mixture viscosity of liquid and gas in cp
        
        input: None

        output:
            mu_m: numeric
        """

        lambda_l, lambda_g = self.input_fraction()
        mu_m = self.mu_oil * lambda_l + self.mu_gas * lambda_g

        return mu_m
    
    def input_fraction(self) -> Tuple[numeric, numeric]:
        """
        get input fraction of gas and liquid

        input: None

        output:
            (lambda_l, lambda_g): (numeric, numeric)
        """

        u_sg, u_sl = self.superficial_velocity("gas"), self.superficial_velocity("liquid")
        u_m = self.mixture_velocity(u_sl, u_sg)

        lambda_l = u_sl / u_m
        lambda_g = 1 - lambda_l

        # print("liquid:", lambda_l)
        # print("gas:", lambda_g)

        return (lambda_l, lambda_g)

    
    def froude_number(self):
        """
        get Froude Number constant
        
        input: None

        output:
            n_fr (Freude Number)      : numeric
        """


        u_sg, u_sl = self.superficial_velocity("liquid"), self.superficial_velocity("gas")
        u_m = self.mixture_velocity(u_sl, u_sg)

        sq_um = pow(u_m, 2)
        
        n_fr = sq_um / (GRAVITY * self.id)

        return n_fr

    def holdup_constant(self) -> Dict[str, List[numeric]]:
        """
        get holdup constant and flow distribution

        input: None

        output: dict
            type: str
            constants: numeric[] (len = 7)
        """

        flow = {
            "type": "",
            "constants": [0] * 7
        }

        # obtain required input fraction
        ## and Froude Number constant
        lambda_l, _ = self.input_fraction()
        n_fr = self.froude_number()

        L1 = 316 * pow(lambda_l, 0.302)
        L2 = 0.000925 * pow(lambda_l, -2.4684)
        L3 = 0.1 * pow(lambda_l, -1.4516)
        L4 = 0.5 * pow(lambda_l, -6.738)

        SEGREGATED = (lambda_l < 0.01 and n_fr < L1) or (lambda_l >= 0.01 and n_fr < L2)
        TRANSITION = lambda_l >= 0.01 and (L2 < n_fr and n_fr <= L3)
        INTERMITTENT = (
            (
                (0.01 <= lambda_l and lambda_l < 0.4) and (L3 < n_fr and n_fr <= L1)
            ) or (
                (lambda_l >= 0.4) and (L3 < n_fr and n_fr <= L4)
            )
        )
        DISTRIBUTED = (lambda_l < 0.4 and n_fr >= L1) or (lambda_l >= 0.4 and n_fr > L4)

        # print(lambda_l, n_fr, SEGREGATED, TRANSITION, INTERMITTENT, DISTRIBUTED)

        if (SEGREGATED):
            flow = {
                "type": "segregated",
                "constants": [0.98, 0.4846, 0.0868, 0.011, -3.768, 3.539, -1.614]
            }
        
        elif (INTERMITTENT):
            flow = {
                "type": "intermittent",
                "constants": [0.845, 0.5351, 0.0173, 2.96, 0.305, -0.4473, 0.0978]
            }

        elif (DISTRIBUTED):
            flow = {
                "type": "distributed",
                "constants": [1.065, 0.5824, 0.0609, 1, 0, 0, 0]
            }

        a, b, c, d, e, f, g = flow["constants"]
        n_vl = self.superficial_velocity("liquid") * pow(self.rho_oil / (32.17 * self.sigma), 0.25)
        flow["C"] = (1 - lambda_l) * log(d * pow(lambda_l, e) * pow(n_vl, f) * pow(n_fr, g))

        C = flow["C"]
        flow["phi"] = 1 + C * (sin(1.8 * self.theta) - 0.333 * pow(sin(1.8 * self.theta), 3))

        return flow
    
    def liquid_holdup(self) -> numeric:
        """
        calculate y_l (liquid holdup)

        input: None

        output:
            y_l: numeric
        """
        
        flow = self.holdup_constant()
        lambda_l, _ = self.input_fraction()
        n_fr = self.froude_number()

        a, b, c, _, _, _, _ = flow["constants"]
        phi = flow["phi"]

        y_lo = a * pow(lambda_l, b) / pow(n_fr, c)
        y_l = y_lo * phi

        return y_l
    
    def reynolds_number_mixture(self) -> numeric:
        """
        calculation of Reynolds Number

        input: None

        output:
            n_re: numeric
        """
        u_sl, u_sg = self.superficial_velocity("liquid"), self.superficial_velocity("gas")
        conversion = 1488
        n_re = self.mixture_density() * self.mixture_velocity(u_sl, u_sg) * self.id * conversion / self.mixture_viscosity()

        return n_re
    
    def friction_factor_no_slip(self) -> numeric:
        """
        calculate no-slip friction factor

        input: None

        output:
            
        """
        n_re = self.reynolds_number_mixture()

        f_ns = 0.0056 + (0.5 / pow(n_re, 0.32))

        return f_ns
    
    def friction_factor_two_phase(self) -> numeric:
        lambda_l, _ = self.input_fraction()
        y_l = self.liquid_holdup()

        f_ns = self.friction_factor_no_slip()

        x = lambda_l / pow(y_l, 2)
        ln_x = log(x)

        A = -0.0523 + 3.182 * ln_x - 0.875 * pow(ln_x, 2) + 0.01853 * pow(ln_x, 4)
        S = ln_x / A

        f_tp = f_ns * exp(S)

        return f_tp
    
    def dp_dz_friction(self):
        """
        calculate pressure gradient based on friction

        input: None

        output:
            dp_dz (psi/ft): numeric
        """

        f_tp = self.friction_factor_two_phase()
        rho_m = self.mixture_density()
        
        u_sl, u_sg = self.superficial_velocity("liquid"), self.superficial_velocity("gas")
        u_m = self.mixture_velocity(u_sl, u_sg)

        conversion = 0.006711409

        dp_dz = 2 * f_tp * rho_m * pow(u_m, 2) * conversion / (GRAVITY * self.id)

        return dp_dz