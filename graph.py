import matplotlib.pyplot as plt

gx = []
gy = []

lines = 0

with open("alphas.txt","r") as f:
 for line in f:
  x,y = line.split(",")
  gx.append(float(x))
  gy.append(float(y))
  lines = lines + 1

plt.plot(gx,gy,'-o')
plt.savefig("img/"+str(lines)+".png")