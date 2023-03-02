import numpy as np
import scipy
import matplotlib.pyplot as plt        #import all the necessary python extensions
from scipy.optimize import curve_fit
from scipy.optimize import minimize

f = open("jpsismall.bin" , "r" )          #open the specified file and read
datalist= np.fromfile(f, dtype=np.float32)          #change from binary
nevent = len(datalist)/7               #divide the data into 7 sections
xdata = np.split(datalist, nevent)         #split the dtata up into these 7 separate sections

xmass = []                # make list of invariant mass of events

for i in range(0 , int(nevent)):              #only take the mass part of the data and append to list
    xmass.append(xdata[i][0])

lamda = 0.002             #initial guesses for my parameters later on
B = 0.5
sigma = 12.74
N = 1

count, bins, normalise = plt.hist(xmass, 200, density=True)     #plots the histogram of the values in a defined number of bins and normalises the fit
index = np.argmax(count)                                  #finds the point at which the max y value occurs
print("Invariant mass = " + str(bins[index]))                  #prints the mass at which the peak y-value for the bin occurs

PeakSig = []                       #create open lists for the sideband subtraction
sidebandL = []
sidebandR = []

for Bars in range(67, 127):               #over the range of the specified bins, append the counts at each bin
    PeakSig.append(count[Bars])
Sigs = sum(PeakSig)                          #total up the number of counts from the range of the peak signal
#print(Sigs)                           (this is an option to check if the answers seem reasonable)
for side1 in range(0, 30):                      #choose the sidebands to be half the peak width and only in a background region
    sidebandL.append(count[Bars])
Lsignal = sum(sidebandL)
#print(Lsignal)
for side2 in range(164, 194):                     #same magnitude and distance from peak for each sideband
    sidebandR.append(count[Bars])
Rsignal = sum(sidebandR)
#print(Rsignal)

SideSubtraction = Rsignal + Lsignal                 #perform the sideband subtraction in steps
TotalSignal = Sigs - SideSubtraction
print("TotalSignal = " + str(TotalSignal))               #print the signal counts with the sidebands removed

Xbar = np.mean(xmass)    #calculates the actual mean from my list of numbers
xbarDiff = abs((bins[index] - Xbar))       #calculates the difference between the input mean and the calculated mean
print("mean= " + str(Xbar))      #prints the mean to check against the input
print("difference in the mean= " + str(xbarDiff))       #print the difference (to quickly check for any errors)

def Curve(bins, sigma, B, N, Xbar, lamda):               #define a curve for plotting later
    gauss = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(bins-Xbar)**2/(2*sigma**2))
    NegExp = (np.exp(-(lamda)*bins))
    return N*(B*gauss+(1-B)*NegExp)             #curve for plotting is a negative exponent and a gaussian

"""
The next method was attempted to be used to calculate both the curve_fit and the 
Negative log likelihood of the Curve function. This was to be used to let the 
program guess the most accurate values for my initial parameters and then try to minimise 
the uncertainty. 
I could not get this to work and im not sure why, as several methods were attempted, but no 
method was without python pulling up errors. The 2 commented out lines below were my initial attempts
but would not plot a curve, so i have left the best case scenario in, where an (inaccurate) curve
has been plotted for the PDF
"""
#scipy.optimize.minimize(Curve(bins, sigma, B, N, Xbar, lamda), x0=(sigma, B, N, Xbar, lamda), method="Nelder-Mead")
#scipy.optimize.curve_fit(Curve, bins, xmass, p0=[sigma, B, N, Xbar, lamda])

def NLL(parameters):                       #Negative log likelihood method for calculating errors 
    sigma, B, N, Xbar, lamda = parameters
    L = -np.sum(np.log(Curve(bins, sigma, B, N, Xbar, lamda)))
    return L

Likelihood = scipy.optimize.minimize(NLL, np.array([sigma, B, N, Xbar, lamda]))         #attempt to minimise the parameters with the NLL

plt.plot(bins, Curve(bins, sigma, B, N, Xbar, lamda), linewidth=2, color='r')        #plot the PDF over the histogram
plt.xlabel('Invariant mass (MeV/c^2)')
plt.ylabel('Probability')
plt.title('PDF for the Invariant mass of Jψ candidate')
plt.grid(True)
plt.savefig("PDFGraph.png")           #save image
plt.show()