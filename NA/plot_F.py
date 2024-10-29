import matplotlib.pyplot as plt

def load_points(filename):
    x, y = [], []
    with open(filename, 'r') as file:
        for line in file:
            px, py = map(float, line.split())
            x.append(px)
            y.append(py)
    return x, y

# Plot each file
plt.figure(figsize=(10, 8))

for m in [10, 40, 160]:
    filename = f"heart_{m}.txt"
    x, y = load_points(filename)
    plt.plot(x, y, label=f'm = {m}')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Heart Shape Approximated with Cubic Bézier Curves')
plt.legend()
plt.grid()
plt.axis('equal')  # To keep the heart shape proportions correct
plt.show()
