import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("timings.csv")
plt.figure(figsize=(8,6))
plt.plot(df['size'], df['MLCG_XOR'] / 1e9, 'o-', label='MLCG-XOR')
plt.plot(df['size'], df['LFG_XOR'] / 1e9, 's-', label='LFG-XOR')
plt.plot(df['size'], df['Xorshift128_MLT'] / 1e9, '^-', label='Xorshift128-MLT')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of generated numbers')
plt.ylabel('Time (seconds)')
plt.title('PRNG Performance Comparison')
plt.legend()
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.savefig('prng_performance.png', dpi=150)
plt.show()