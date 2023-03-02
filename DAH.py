import random
import numpy as np
import matplotlib.pyplot as plt
import statistics                       #import all the necessary python extensions

randomlist = []                     #create an open array for my randomly chosen numbers
for i in range(0,100):
    n = random.randint(0,100)              #randomly choose a number between 1 and 100, 100 times
    randomlist.append(n)               #append the randomly chosen number to the array and repeat the loop
    
print(randomlist)                       #print values of the array to calaculate the mean and standard deviation
mu = statistics.mean(randomlist)                     #calculates the mean value
print(mu) 
sigma = statistics.stdev(randomlist, xbar = mu)             #calculates the standard deviation
print(sigma)

"""sets values for what the histogram will plot- the list, the number of bins and i

"""
count, bins, ignored = plt.hist(randomlist, 100, density=True)            
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-(bins-mu)**2/(2*sigma**2)), linewidth=2, color='r')
plt.xlabel('Randomly generated value')
plt.ylabel('Probability density of the Gaussian')
plt.title('Gaussian of a Set of Random Numbers')
plt.grid(True)
plt.show()
#plt.savefig("figure.pdf")         #saves the graph as a PDF