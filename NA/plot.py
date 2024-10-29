import matplotlib.pyplot as plt

def parse_output(filename):
    data = {}
    with open(filename, 'r') as file:
        n = None
        x_vals = []
        exact_vals = []
        interp_vals = []
        
        for line in file:
            line = line.strip()
            if line.startswith("n = "):
                # 将当前数据存储到字典中
                if n is not None:
                    data[n] = {'x': x_vals, 'exact': exact_vals, 'interp': interp_vals}
                
                # 解析 n 值并重置临时存储
                n = int(line.split('=')[1].strip(':'))
                x_vals = []
                exact_vals = []
                interp_vals = []
            elif line and len(line.split()) >= 6:
                # 解析 x、Exact 和 Interpolated 值
                parts = line.split()
                try:
                    x = float(parts[1])
                    exact = float(parts[3])
                    interp = float(parts[5])
                    x_vals.append(x)
                    exact_vals.append(exact)
                    interp_vals.append(interp)
                except ValueError:
                    print(f"Skipping invalid line: {line}")
        
        # 存储最后一个 n 的数据
        if n is not None:
            data[n] = {'x': x_vals, 'exact': exact_vals, 'interp': interp_vals}
    
    return data


def plot_interpolations(data):
    plt.figure(figsize=(10, 6))
    
    for n, values in data.items():
        plt.plot(values['x'], values['exact'], label=f'Exact function (n={n})', color='black', linestyle='--')
        plt.plot(values['x'], values['interp'], label=f'Newton interpolation (n={n})')
    
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.title('Newton Interpolation and Runge Phenomenon')
    plt.grid()
    plt.show()

filename = 'Newton_output.txt'
data = parse_output(filename)
plot_interpolations(data)

filename = 'Chebyshev_output.txt'
data = parse_output(filename)
plot_interpolations(data)
