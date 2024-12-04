import numpy
import math
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

class ElectricVehicleModel:
    def __init__(self):
        '''
        32650 LiFePO4
        The number of batteries connected in series is 100, 
        and the number of batteries connected in parallel is 20

        '''
        self.M_bat = 40  # Battery thermal mass (kg)
        self.capacity_bat_thermal = 1350  # Battery specific heat capacity (J/kg·K)
        self.A_bat = 12   # 电池包的表面积(m^2)
        self.entropy_coefficient = 0.06  # 电池的熵系数 (V/K)

        self.U_oc = 320 # 假定开路电压 (V)
        self.R_bat = 6 * 1e-3  # 假定电池包内阻 (Ω)

        '''
        Coolant
        '''
        self.rho_clnt = 1069.5  # Density of coolant (kg/m^3)
        self.capacity_clnt = 3330 # Specific heat capacity of coolant (J/kg·K)
        self.V_pump = 33*1e-6 # Displacement volume of pump (m^3/rev)
        self.h_bat = 300 # Heat transfer coefficient between the battery and the coolant (W/m^2/K)
        self.A_bat = 1 # Heat transfer sectional area between the battery and the coolant (m^2)
        
        '''
        Refrigerant circle
        R-134a
        Compressor, Condenser, Expansion valve, Evaporator
        '''
        self.rho_rfg = 27.8   # Density of refrigerant (kg/m^3)
        self.capacity_rfg = 1117 # Specific heat capacity of refrigerant (J/kg·K)
        self.V_comp = 33*1e-6 # Displacement volume of compressor (m^3/rev)
        self.h_eva = 1000 # Heat transfer coefficient between evaporator and refrigerant/coolant (W/m^2/K)
        self.A_eva = 0.3 # Heat transfer sectional area between evaporator and refrigerant/coolant (m^2)
        self.PR = 5 # Compression ratio of the compressor

        self.h_comp_out = 284.3 # Enthalpy at the outlet of compressor (kJ/Kg)
        self.P_comp_out = 1500  # Pressure at the outlet of compressor (kPa)

        self.h_eva_out = 250.9  # Enthalpy at the outlet of evaporator (kJ/Kg)
        self.P_comp_in = 300    # Pressure at the inlet of compressor (kPa)

        self.h_cond_out = 131.7 # Enthalpy at the outlet of condenser (kJ/Kg)
        self.P_cond_out = 1500  # Pressure at the outlet of condenser (kPa)

        self.P_eva_in = 300

        
        '''
        Control variables
        '''
        self.massflow_rfg = 0   # Mass flow rate of refrigerant (Kg/s)
        self.massflow_clnt = 0  # Mass flow rate of coolant (Kg/s)
        self.T_clnt_eva_in = 22 # 假定的初始的进入蒸发器的冷却液温度 (℃)

        '''
        Vehicle_Dynamics_Model
        '''
        self.v_previous = 0 # 储存上一时刻的车速，初始化从0开始
        self.rho_air = 1.225 # 空气密度 (Kg/m^3)

        '''
        中间变量储存
        '''
        self.P_cooling = 0 # 储存电池发热模型中计算的中间变量冷却系统功率P_cooling

    def compute_massflow_rfg(self, omega_comp):
        eta_comp = 1e-5 * omega_comp + 0.9
        self.massflow_rfg = self.V_comp * eta_comp * omega_comp * self.rho_rfg * 2 * math.pi / 60

    def power_compressor(self, omega_comp):
        eta_isen = 0.75

        s_eva_out = PropsSI('S', 'P', self.P_comp_in * 1e3, 'H', self.h_eva_out * 1e3, 'R134a')  # 入口熵
        h_isen = PropsSI('H', 'P', self.P_comp_out * 1e3, 'S', s_eva_out, 'R134a') * 1e-3 # 等熵出口焓值
        self.h_comp_out = self.h_eva_out + (h_isen - self.h_eva_out)/eta_isen
        if omega_comp == 0:
            P_comp = 0
        else: 
            torque_comp = self.massflow_rfg * (self.h_comp_out - self.h_eva_out) / eta_isen / omega_comp
            eta_c = max(1, 0.02 * torque_comp + 0.6)
            P_comp = torque_comp * omega_comp / eta_c
        return P_comp
    
    def condenser(self):
        Q_c = self.massflow_rfg * (self.h_comp_out - self.h_cond_out)
        '''
        冷凝板待完善，需要添加风扇散热模型
        P_fan = alpha_fan * ( m_{a,fan} - m^{lb}_{a,fan})^2
        用到radiator函数
        目前是假设散热能力足够，可以达到饱和液体状态，同时不考虑风扇能耗
        '''
        return Q_c
    
    def compute_massflow_clnt(self, omega_pump):
        eta_pump = 1e-5 * omega_pump + 0.9
        self.massflow_clnt = self.V_pump * eta_pump * omega_pump * self.rho_clnt * 2 * math.pi / 60

    def evaporator(self):
        self.P_eva_in = self.P_cond_out - (self.P_comp_out - self.P_comp_in)
        T_rfg_eva_in = PropsSI('T', 'P', self.P_eva_in * 1e3, 'H', self.h_cond_out * 1e3, 'R134a') -273.15

        if self.massflow_clnt == 0:
            A = 0.5
            B = 0.5
        else:
            k1 = self.h_eva * self.A_eva / self.massflow_rfg / self.capacity_rfg
            k2 = self.h_eva * self.A_eva / self.massflow_clnt / self.capacity_clnt
            A = (k2 + k1 * math.exp(-(k1 + k2))) / (k1 + k2)
            B = (k1 + k2 * math.exp(-(k1 + k2))) / (k1 + k2)

        T_clnt_eva_out = A * self.T_clnt_eva_in + (1 - A) * T_rfg_eva_in
        T_rfg_eva_out = B * T_rfg_eva_in + (1 - B) * self.T_clnt_eva_in

        self.h_eva_out = PropsSI('H', 'P', self.P_comp_in * 1e3, 'T', T_rfg_eva_out + 273.15, 'R134a') * 1e-3 # 更新制冷循环中蒸发器出口、压缩机进口的焓值
        
        return T_clnt_eva_out
    
    def battery_cooling(self, T_bat):
        T_clnt_in = self.evaporator()
        if self.massflow_clnt == 0:
            T_clnt_out = T_bat
        else:
            T_clnt_out = (T_clnt_in - T_bat) * math.exp(-(self.h_bat * self.A_bat) / (self.massflow_clnt * self.capacity_clnt)) + T_bat
        self.T_clnt_eva_in = T_clnt_out  # 更新冷却循环中冷却液到电池吸热后出口的温度
        Q_cool = self.massflow_clnt * self.capacity_clnt * (T_clnt_in - T_clnt_out)
        return Q_cool
    
    def power_pump(self, omega_pump):
        Delta_pump = self.massflow_clnt ** 2  # Pressure Drop (kPa)
        if omega_pump == 0:
            P_pump = 0
        else:
            torque_pump = Delta_pump * self.massflow_clnt / self.rho_clnt / omega_pump
            eta_p = max(1, 0.02 * torque_pump + 0.6)
            P_pump = torque_pump * omega_pump / eta_p
        return P_pump
        
    
    def battery_thermal_generation(self, T_bat, velocity, omega_comp, omega_pump):
        P_trac = self.Vehicle_Dynamics_Model(velocity)
        P_comp = self.power_compressor(omega_comp)
        P_pump = self.power_pump(omega_pump)
        # print('牵引功率, 压缩机功率, 泵功率 分别为', P_trac, P_comp, P_pump)
        self.P_cooling = P_pump + P_comp
        P_bat = P_trac + self.P_cooling
        
        I_bat = (self.U_oc - math.sqrt(self.U_oc ** 2 - 4 * self.R_bat * P_bat)) / (2 * self.R_bat)
        '''
        不考虑摩擦制动 Assume that all the brakings can be accommodated through electric brakes
        Note that the regenerative braking, unlike the friction braking, 
        directly affects the battery thermal response during the battery charging.
        实践中也尽可能减少摩擦磨损
        再生制动时牵引功率可能是负的, 导致整个电池输出功率为负, 电流为负
        此处电池充放电发热公式一致, 发热公式需取电流的绝对值
        且SOC, SOH计算方法需将电流正负考虑在内
        '''
        #print('电池电流I = ', I_bat)
        #print('电池内阻发热量 = ', I_bat**2 * self.R_bat)
        #print('电池的可逆热 = ', I_bat * (T_bat + 273.15) * self.entropy_coefficient)
        Q_gen = I_bat**2 * self.R_bat + abs(I_bat) * (T_bat + 273.15) * self.entropy_coefficient
        return Q_gen

    def battery_thermal_model(self, T_bat, velocity, omega_comp, omega_pump):
        Q_cool = self.battery_cooling(T_bat)
        print('电池散热量: ', Q_cool)
        Q_gen = self.battery_thermal_generation(T_bat, velocity, omega_comp, omega_pump)
        print('电池发热量: ', Q_gen)
        T_bat_next = T/(self.M_bat * self.capacity_bat_thermal)*(Q_cool + Q_gen) + T_bat
        return T_bat_next

    def Vehicle_Dynamics_Model(self, v):
        delta = 1.13  # Correction coefficient for rotating mass
        m_veh = 2200
        f = 0.019  # Rolling resistance coefficient
        C_d = 0.3  # Aerodynamic drag coefficient
        A_wind = 3.2  # Vehicle frontal area (m^2)
        eta_p = 0.9  # Powertrain efficiency
        g = 9.8
        a = (v - self.v_previous)/T
        F_r = f * m_veh * g
        F_a = 0.5 * C_d * A_wind * self.rho_air * v**2
        if a >= 0:
            P_trac = (v * (delta * m_veh * a + F_r + F_a)) / eta_p
        else:
            P_trac = (v * (delta * m_veh * a + F_r + F_a)) * eta_p
        # print('此时的速度, 加速度, 牵引功率为', v , a, P_trac)
        self.v_previous = v  # 储存这个时刻的速度v, 来计算下一时刻的加速度
        return P_trac 

    def aging_model(self):
        return 0
    

def movement(time):
    a = 1
    if time <= 10:
        v = a * time
    else:
        v = 10 + math.sin(time)
    
    return v

T = 1  # 采样时间 
t = 0  # 系统开始的时间

T_bat = 30
T_thres_upper=32
T_thres_lower=28
EV = ElectricVehicleModel()
omega_comp = 0.
omega_pump = 0.

time = []
battery_temperature = []
# 用于存储制冷循环的数据
cycle_data = []  # 每个元素存储一组 [h_eva_out, P_eva_out, h_comp_out, P_comp_out, h_cond_out, P_cond_out, h_eva_in, P_eva_in]


while t < 1125:
    v = movement(t)

    if T_bat > T_thres_upper:  # 若电池温度大于目标温度，启动制冷冷却循环
        omega_comp = 3000.0
        omega_pump = 2000.0
    elif T_bat < T_thres_lower: # 若电池温度小于最低温度，关闭制冷冷却循环
        omega_comp = 0.
        omega_pump = 0.
    
    EV.compute_massflow_rfg(omega_comp)  # 更新制冷循环质量流量
    EV.compute_massflow_clnt(omega_pump) # 更新冷却循环质量流量
    

    T_bat = EV.battery_thermal_model(T_bat, v, omega_comp, omega_pump)
    '''
    更新温度: 调用了battery_cooling和battery_thermal_generation
        - battery_thermal_generation计算了电池的发热量
        - battery_cooling调用了evaporator()函数, 更新冷却循环中冷却液到电池吸热后出口的温度T_clnt_eva_in
            ·evaporator()函数中更新了制冷循环中的蒸发器出口焓值h_rfg_eva_out
                h_rfg_eva_out影响压缩机的功率P_comp以及压缩机出口焓值h_comp_out
                    h_comp_out影响冷凝板功率, 冷凝板需要将h_comp_out降低到h_cond_out, 即刚好是饱和液体
            ·更新后的T_clnt_eva_in会影响下一时刻的evaporator()中制冷循环eva出口焓值以及冷却循环电池冷却液入口温度
    '''
    print('当时间为: ', t, ' 时 ', ' 电池温度为 ', T_bat)

    # 存储时间和电池温度数据用于绘图
    time.append(t)
    battery_temperature.append(T_bat)

    # 如果制冷循环启动，记录制冷循环中的焓值和压力
    if omega_comp > 0:
        cycle_data.append([
            EV.h_eva_out, EV.P_comp_in,  # 蒸发器出口
            EV.h_comp_out, EV.P_comp_out,  # 压缩机出口
            EV.h_cond_out, EV.P_cond_out,  # 冷凝器出口
            EV.h_cond_out, EV.P_eva_in  # 蒸发器入口
        ])

    t += T

# --- 第一部分：绘制电池温度随时间的变化 ---
plt.figure(figsize=(10, 5))
plt.plot(time, battery_temperature, label='Battery Temperature')
plt.axhline(y=T_thres_upper, color='r', linestyle='--', label='Upper Threshold')
plt.axhline(y=T_thres_lower, color='b', linestyle='--', label='Lower Threshold')
plt.xlabel('Time (s)')
plt.ylabel('Battery Temperature (°C)')
plt.title('Battery Temperature vs. Time')
plt.legend()
plt.grid()
plt.show()

# --- 第二部分：绘制制冷循环的焓值-压力图 ---
# 确保记录的 cycle_data 至少有 4 组数据
if len(cycle_data) >= 4:
    plt.figure(figsize=(10, 8))
    colors = ['r', 'g', 'b', 'm']  # 用不同颜色绘制四个循环

    for i in range(4):  # 取前 4 次制冷循环的数据
        data = cycle_data[i]
        h_values = [data[0], data[2], data[4], data[6], data[0]]  # 按顺序连接蒸发器出口、压缩机出口、冷凝器出口、蒸发器入口
        P_values = [data[1], data[3], data[5], data[7], data[1]]  # 对应的压力值
        
        plt.plot(h_values, P_values, marker='o', color=colors[i], label=f'Cycle {i+1}')

    plt.xlabel('Enthalpy (h) [J/kg]')
    plt.ylabel('Pressure (P) [Pa]')
    plt.title('Refrigeration Cycle (h-P Diagram)')
    plt.legend()
    plt.grid()
    plt.show()
else:
    print("制冷循环未满足启动条件或数据不足，无法绘制焓值-压力图。")