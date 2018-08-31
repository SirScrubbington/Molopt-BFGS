import os
import sys
import matplotlib.pyplot as plt

titles = [
			"Atom By Atom Optimisation - Line Start",
			"Atom By Atom Optimisation - Circle Start",
			"Atom By Atom Optimisation - Square Start",
			"Iterative Atom Optimisation - Line Start",
			"Iterative Atom Optimisation - Circle Start",
			"Iterative Atom Optimisation - Square Start",
			"Brute Force Atom Optimisation - Line Start",
			"Brute Force Atom Optimisation - Circle Start",
			"Brute Force Atom Optimisation - Square Start",
			"Atom By Atom Optimisation - Line Start with BFGS",
			"Atom By Atom Optimisation - Circle Start with BFGS",
			"Atom By Atom Optimisation - Square Start with BFGS",
			"Iterative Atom Optimisation - Line Start with BFGS",
			"Iterative Atom Optimisation - Circle Start with BFGS",
			"Iterative Atom Optimisation - Square Start with BFGS",
			"Brute Force Atom Optimisation - Line Start with BFGS",
			"Brute Force Atom Optimisation - Circle Start with BFGS",
			"Brute Force Atom Optimisation - Square Start with BFGS"
		 ]

def graphCSV(filename):
	if os.path.isfile(filename+".png"):
		return

	with open(filename,"r") as fhandle:
		for file in fhandle:
			data = file.split(",")
			
			x = []
			y = []
			
			trial = int(data[0])
			n = int(data[1])
			
			iter = data[-4]
			opt = data[-3]
			cost = data[-2]
			time = data[-1]
			
			for i in range(2,n+2):
				ix,iy = data[i].split(":")
				x.append(float(ix))
				y.append(float(iy))
		
			plt.suptitle(titles[trial])
			plt.title("Atoms: " + str(n) + " Cost: " + str(cost) + " Best: " + opt + " Time: " + str(time) + "Iters: " + iter,y=1,fontsize=8)
		
			plt.xlabel("x coordinates")
			plt.ylabel("y coordinates")
			
			plt.plot(x,y,'-o')
			
			plt.savefig(filename+".png")
			plt.clf()
		
def checkDirectory(path):
	files = os.listdir(path);
	for file in files:
		if os.path.isdir(path + "/" + file):
			checkDirectory(path + "/" + file)
		else:
			if file.endswith(".csv"):
				graphCSV(path + "/" + file)

if __name__ == '__main__':
	if len(sys.argv) > 1:
		checkDirectory(sys.argv[1]);
	else:
		print("No directory provided!");
