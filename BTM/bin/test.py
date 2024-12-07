from CoolProp.CoolProp import PropsSI

# 测试 PropsSI 函数
pressure = 1.32e6  # 压力，单位：Pa
enthalpy = 272.52e3  # 焓值，单位：J/kg
temperature = PropsSI('T', 'P', pressure, 'H', enthalpy, 'R134a')
Q_rfg_cond_out = PropsSI('Q', 'P', 2100 * 1e3, 'H', 300 * 1e3, 'R134a')
Q_rfg_eva_in = PropsSI('Q', 'P', 600 * 1e3, 'H', 300 * 1e3, 'R134a')
Q_rfg_comp_in = PropsSI('Q', 'P', 600 * 1e3, 'H', 410 * 1e3, 'R134a')
Q_rfg_cond_in = PropsSI('Q', 'P', 2100 * 1e3, 'H', 430 * 1e3, 'R134a')

T_rfg_cond_out = PropsSI('T', 'P', 2100 * 1e3, 'H', 300 * 1e3, 'R134a') - 273.15
T_rfg_eva_in = PropsSI('T', 'P', 600 * 1e3, 'H', 300 * 1e3, 'R134a') - 273.15
T_rfg_comp_in = PropsSI('T', 'P', 600 * 1e3, 'H', 410 * 1e3, 'R134a') - 273.15
T_rfg_cond_in = PropsSI('T', 'P', 2100 * 1e3, 'H', 430 * 1e3, 'R134a') - 273.15
print(f"温度为: {T_rfg_cond_out:.2f} C, 状态为: {Q_rfg_cond_out:.2f}")
print(f"温度为: {T_rfg_eva_in:.2f} C, 状态为: {Q_rfg_eva_in:.2f}")
print(f"温度为: {T_rfg_comp_in:.2f} C, 状态为: {Q_rfg_comp_in:.2f}")
print(f"温度为: {T_rfg_cond_in:.2f} C, 状态为: {Q_rfg_cond_in:.2f} ")



