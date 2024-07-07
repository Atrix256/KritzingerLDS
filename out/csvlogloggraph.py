import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

df = pd.read_csv(sys.argv[1]).drop(['Samples'], axis=1)





fig, ax = plt.subplots()

ax.plot(df['Kritzinger'], label="Kritzinger")
ax.plot(df['GoldenRatio'], label="GoldenRatio")
ax.plot(df['White'], label="White")
ax.plot(df['Stratified'], label="Stratified")
ax.plot(df['Regular'], label="Regular")
ax.plot(df['RegularWalls'], label="RegularWalls")

ax.legend()

plt.title('RMSE for ' + sys.argv[1])
plt.ylabel('RMSE')
plt.xlabel('Samples')

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig(sys.argv[1] + ".png", bbox_inches='tight')





fig, ax = plt.subplots()

ax.plot(df['Kritzinger'], label="Kritzinger")
ax.plot(df['GoldenRatio'], label="GoldenRatio")
ax.plot(df['White'], label="White")
#ax.plot(df['Stratified'], label="Stratified")
#ax.plot(df['Regular'], label="Regular")
#ax.plot(df['RegularWalls'], label="RegularWalls")

ax.legend()

plt.title('RMSE for ' + sys.argv[1])
plt.ylabel('RMSE')
plt.xlabel('Samples')

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig(sys.argv[1] + ".seq.png", bbox_inches='tight')




fig, ax = plt.subplots()

#ax.plot(df['Kritzinger'], label="Kritzinger")
#ax.plot(df['GoldenRatio'], label="GoldenRatio")
#ax.plot(df['White'], label="White")
ax.plot(df['Stratified'], label="Stratified")
ax.plot(df['Regular'], label="Regular")
ax.plot(df['RegularWalls'], label="RegularWalls")

ax.legend()

plt.title('RMSE for ' + sys.argv[1])
plt.ylabel('RMSE')
plt.xlabel('Samples')

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig(sys.argv[1] + ".set.png", bbox_inches='tight')



fig, ax = plt.subplots()

#ax.plot(df['Kritzinger'], label="Kritzinger")
#ax.plot(df['GoldenRatio'], label="GoldenRatio")
#ax.plot(df['White'], label="White")
#ax.plot(df['Stratified'], label="Stratified")
ax.plot(df['Regular'], label="Regular")
ax.plot(df['RegularWalls'], label="RegularWalls")

ax.legend()

plt.title('RMSE for ' + sys.argv[1])
plt.ylabel('RMSE')
plt.xlabel('Samples')

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig(sys.argv[1] + ".reg.png", bbox_inches='tight')