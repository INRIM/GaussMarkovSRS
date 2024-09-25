import decimal

import numpy as np


def unc_to_decimal_exp(u, digits=3):
    exp = decimal.Decimal(10) ** (int(np.floor(np.log10(abs(float(u))))) - digits + 1)
    return exp


def decimal_to_string(d, exp):
    return f"{d.quantize(exp)}"


if __name__ == "__main__":
    x = unc_to_decimal_exp(decimal.Decimal(1.42123), 3)
    print(decimal_to_string(decimal.Decimal("1000.123456789"), x))
