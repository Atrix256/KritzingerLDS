import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

df = pd.read_csv(sys.argv[1]).drop(['Samples'], axis=1)





fig, ax = plt.subplots()

ax.plot(df['Kritzinger'], label="Kritzinger")
ax.plot(df['GoldenRatio'], label="GoldenRatio")
ax.plot(df['Regular'], label="Regular")
ax.plot(df['RegularWalls'], label="RegularWalls")
ax.plot(df['Stratified'], label="Stratified")

ax.legend()

plt.title('Abs Error for ' + sys.argv[1])
plt.ylabel('Abs Error')
plt.xlabel('Samples')

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig(sys.argv[1] + ".png", bbox_inches='tight')





fig, ax = plt.subplots()

ax.plot(df['Kritzinger'], label="Kritzinger")
ax.plot(df['GoldenRatio'], label="GoldenRatio")
#ax.plot(df['Regular'], label="Regular")
#ax.plot(df['RegularWalls'], label="RegularWalls")
#ax.plot(df['Stratified'], label="Stratified")

ax.legend()

plt.title('Abs Error for ' + sys.argv[1])
plt.ylabel('Abs Error')
plt.xlabel('Samples')

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig(sys.argv[1] + ".seq.png", bbox_inches='tight')
