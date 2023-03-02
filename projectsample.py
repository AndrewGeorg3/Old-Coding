import numpy as np
import matplotlib.pyplot as plt        #import all the necessary python extensions
# import data
# xmass = np.loadtxt(sys.argv[1])

f = open("jpsismall.bin" , "r" )
datalist= np.fromfile(f, dtype=np.float32)
# number of events
nevent = len(datalist)/7
#print(nevent)
xdata = np.split(datalist, nevent)

# make list of invariant mass of events
xmass = []
xmom = []
rapidity = []
chiSq = []
minMom = []
minProbNNmu = []
minImpact = []

for i in range(0 , int(nevent)):
    xmass.append(xdata[i][0])
    xmom.append(xdata[i][1])
    rapidity.append(xdata[i][2])
    chiSq.append(xdata[i][3])
    minMom.append(xdata[i][4])
    minProbNNmu.append(xdata[i][5])
    minImpact.append(xdata[i][6])
        
count, bins, normalise = plt.hist(xmass, 200, density=False)     #plots the histogram of the values in a defined number of bins and normalises the gaussian fit
index = np.argmax(count)
print("Invariant mass = " + str(bins[index]))
plt.xlabel('Invariant mass (MeV/c^2)')
plt.ylabel('Candidates')
plt.title('Histogram of the Invariant mass of Jψ candidate')
plt.grid(True)
plt.savefig("MassGraph.png")
plt.show()

count, bins, normalise = plt.hist(xmom, 200, density=False)
plt.xlabel('Transverse Momentum (MeV/c)')
plt.ylabel('Candidates')
plt.title('Histogram of the Transverse Momentum of the di-muon candidate')
plt.grid(True)
plt.savefig("TransMomentumGraph.png") 
plt.show()
     
count, bins, normalise = plt.hist(rapidity, 200, density=False)
plt.xlabel('Rapidity')
plt.ylabel('Candidates')
plt.title('Histogram of the Rapidity (η) of the di-muon candidate')
plt.grid(True)
plt.savefig("RapidityGraph.png") 
plt.show()

count, bins, normalise = plt.hist(chiSq, 100, density=False)
plt.xlabel('Chi-squared')
plt.ylabel('Candidates')
plt.title('Histogram of Chi-squared of the geometric vertex of the dimuon candidate')
plt.grid(True)
plt.savefig("ChiSqGraph.png")
plt.show() 

count, bins, normalise = plt.hist(minMom, 200, density=False)
plt.xlabel('Minimum transverse momentum (MeV/c)')
plt.ylabel('Candidates')
plt.title('Histogram of the Minimum Transverse Momentum of the two muons')
plt.grid(True)
plt.savefig("minMomentumGraph.png") 
plt.show()

count, bins, normalise = plt.hist(minProbNNmu, 200, density=False)
plt.xlabel('Probability that a muon is a muon')
plt.ylabel('Candidates')
plt.title('Histogram of how well the candidate matches the hypothesis')
plt.grid(True)
plt.savefig("minProbNNmugraph.png") 
plt.show()

count, bins, normalise = plt.hist(minImpact, 50, density=False)
plt.xlabel('Minimum impact parameter, Chi-squared, of the two muons')
plt.ylabel('Candidates')
plt.title('Does the particle originate from a proton-proton interaction?')
plt.grid(True)
plt.savefig("minImpactgraph.png") 
plt.show()