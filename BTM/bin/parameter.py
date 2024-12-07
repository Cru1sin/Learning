import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator

def parameter_of_cooling(air_velocity_input, coolant_flow_input):
    # 图中提取的实验数据（大致估计值，实际可用更精确数据）
    air_velocity = [0, 2, 4, 6, 8, 10]  # Air velocity (m/s)
    coolant_flow = [0, 50, 100, 150, 250]  # Coolant flow (L/h)
    heat_exchange = [  # Heat exchange (W x 10^4)
        [200, 300, 350, 400, 450, 500],  # For coolant flow = 0 L/h
        [210, 320, 380, 450, 480, 550],  # For coolant flow = 50 L/h
        [215, 320, 390, 460, 500, 600],  # For coolant flow = 100 L/h
        [215, 320, 390, 480, 550, 620],  # For coolant flow = 150 L/h
        [200, 320, 390, 500, 600, 620],  # For coolant flow = 250 L/h
    ]

    air_velocity_array = np.array(air_velocity)
    coolant_flow_array = np.array(coolant_flow)
    heat_exchange_array = np.array(heat_exchange)

    interpolator = RegularGridInterpolator(
        (air_velocity_array, coolant_flow_array),
        heat_exchange_array.T,  # 注意矩阵需要转置匹配坐标
        method='cubic'
    )

    heat_exchange_pred = interpolator([air_velocity_input, coolant_flow_input])

    return heat_exchange_pred

from CoolProp.CoolProp import PropsSI
import math
def coefficient_heat_exchange(T, P, massflow):
    k = PropsSI('L', 'T', T, 'P', P, 'R134a')  # Thermal conductivity (W/m·K)
    mu = PropsSI('V', 'T', T, 'P', P, 'R134a')  # Dynamic viscosity (Pa·s)
    cp = PropsSI('C', 'T', T, 'P', P, 'R134a')  # Specific heat (J/kg·K)
    rho = PropsSI('D', 'T', T, 'P', P, 'R134a')  # Density (kg/m³)
    print("rho = ", rho)
    # Step 2: Define flow parameters
    D = 0.01  # Pipe diameter in meters
    v = massflow / (math.pi/4 * rho * D**2)

    # Step 3: Calculate Reynolds and Prandtl numbers
    Re = rho * v * D / mu
    Pr = cp * mu / k

    # Step 4: Use ht library to compute Nusselt number and heat transfer coefficient
    # Step 4: Determine the flow regime and calculate Nu
    if Re < 2300:  # Laminar flow
        Nu = 3.66  # Fully developed laminar flow
    elif Re > 4000:  # Turbulent flow
        Nu = 0.023 * (Re**0.8) * (Pr**0.3)  # Dittus-Boelter equation
    h = Nu * k / D  # Convective heat transfer coefficient
    return h