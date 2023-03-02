import random
import numpy as np
import matplotlib.pyplot as plt        #import all the necessary python extensions

mu, sigma = 58, 5  #mean and standard deviation I have set as inputs
numbers = np.random.normal(mu, sigma, 1000)      #randomly choose 1000 values from a gaussian of specified mean and sigma

Xbar = np.mean(numbers)    #calculates the actual mean from my list of numbers
xbarDiff = abs(mu - Xbar)       #calculates the difference between the input mean and the calculated mean
print("mean= " + str(Xbar))      #prints the mean to check against the input
print("difference in the mean= " + str(xbarDiff))       #print the difference

deviation = np.std(numbers)         #calculates the actual standard deviation from my list of numbers
devDiff = abs(deviation - sigma)                    #calculates the difference between the input and the calculated value
print("standard deviation= " + str(deviation))         #prints the standard deviation to check against the input
print("difference in the standard deviation= " + str(devDiff))             #print the difference in the standard deviations
  
count, bins, normalise = plt.hist(numbers, 200, density=True)     #plots the histogram of the values in a defined number of bins and normalises the gaussian fit
plt.plot(bins, 1/(deviation * np.sqrt(2 * np.pi)) * np.exp(-(bins-Xbar)**2/(2*deviation**2)), linewidth=2, color='r')  #plots a red line of best fit for the gaussian
plt.xlabel('Randomly generated value')
plt.ylabel('Probability density of the Gaussian')
plt.title('Gaussian of a Set of Random Numbers')
plt.grid(True)
plt.show()
plt.savefig("Q5graph.pdf")         #saves the graph as a PDF 