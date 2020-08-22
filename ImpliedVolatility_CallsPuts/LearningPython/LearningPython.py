import numpy as np
from scipy.stats import norm
import math


class option:
    '''parameters and prices of options'''
    def __init__(self, S_0, T, r, K, n, sigma):   
        self.strprice = S_0
        self.T = T
        self.r = r
        self.K = K
        self.n = n
        self.sigma = sigma

    def euroPutPrice(self):
        '''Return put intial value'''
        d1 = ((self.r+(self.sigma**2)/2)*self.T+np.log(self.strprice/self.K))/(self.sigma*math.sqrt(self.T))
        d2 = ((self.r-(self.sigma**2)/2)*self.T+np.log(self.strprice/self.K))/(self.sigma*math.sqrt(self.T))
        p_e = math.e**(-self.r*T)*self.K*(1-norm.cdf(d2))-self.strprice*(1-norm.cdf(d1))
        return p_e

    def payoff(self, stock_price):
        '''Return list of exercise payoffs'''
        y = []
        for s in stock_price:
            x = self.K-s
            y.append(x)

    
        for i in range(len(y)):  
            y[i] = max(y[i],0)    
        
        return y

    def usaPutPrice(self):
        '''Return American put intial value'''
        deltaT = self.T/self.n 
        r_n = math.e**(self.r*deltaT)-1

        u_n = math.e**(self.sigma*math.sqrt(deltaT))
        d_n = 1/u_n
        p = ((1+r_n)-d_n)/(u_n-d_n) 

        treeStock = []

        for t in range(self.n+1):
            for i in range(t+1):
                treeStock.append(self.strprice*u_n**i*d_n**(t-i))

        payoffs = self.payoff(treeStock)


        newestStep = payoffs[-(self.n+1):]
        values = []
        values.append(newestStep)
        binary = []
        step = 1
        partition = self.n+1
        while (self.n+1-step) != 0:
            binny = []
            binary_temp = []

            partition = partition + (self.n+1-step)
            for i in range(self.n+1-step):
                down = newestStep[i]
                up = newestStep[i+1]
                calc_E = (1/(1+r_n))*(p*up+(1-p)*down)
                if calc_E > payoffs[-(partition-i)]:
                    binny.append(calc_E)
                    binary_temp.append(0)
                else:
                    binny.append(payoffs[-(partition-i)])
                    binary_temp.append(1)
            newestStep = binny
            values.insert(0,binny)
            binary.insert(0,binary_temp)
            step = step +1
        return values[0][0]
    
    def callPrice(self):
        '''Return European call intial value'''
        # European and American calls are equivalent products
        F = self.strprice - self.K*math.exp(-self.r*self.T)
        return self.euroPutPrice() + F   # using call-put parity

def implVol_europut(contract, mrktprice, precision=5):
    '''Return implied volatily of a European put option'''
    threshold = 10**-(precision+1)
    opt = contract
    sigma_guess = contract.sigma # place the guess of sigma in the contract object
    p_guess = contract.euroPutPrice()
    if p_guess < mrktprice:
        sigma_low = sigma_guess
        p_low = p_guess
        j = 1
        p_high = p_guess # dummy intialization
        while p_high < mrktprice:
            sigma_high = 2**j*sigma_guess
            opt.sigma = sigma_high
            p_high = opt.euroPutPrice()
            j = j + 1 
    else:
        sigma_high = sigma_guess
        p_high = p_guess
        j=1
        p_low = p_guess # dummy intialization
        while p_low > mrktprice:
            sigma_low = 2**(-j*sigma_guess)
            opt.sigma = sigma_low
            p_low = opt.euroPutPrice()
            j = j + 1 

    while sigma_high-sigma_low > threshold: #if the interval length is too big
        sigma_mid = 0.5*(sigma_high+sigma_low)
        opt.sigma = sigma_mid
        p_mid = opt.euroPutPrice()
        if p_mid < mrktprice:
            sigma_low = sigma_mid
        else:
            sigma_high = sigma_mid
    return round(0.5*(sigma_high+sigma_low), precision)

def implVol_usaput(contract, mrktprice, precision=5):
    '''Return implied volatily of a American put option'''
    threshold = 10**-4
    opt = contract
    sigma_guess = contract.sigma # place the guess of sigma in the contract object
    p_guess = contract.usaPutPrice()
    if p_guess < mrktprice:
        sigma_low = sigma_guess
        p_low = p_guess
        j = 1
        p_high = p_guess # dummy intialization
        while p_high < mrktprice:
            sigma_high = 2**j*sigma_guess
            opt.sigma = sigma_high
            p_high = opt.usaPutPrice()
            j = j + 1 
    else:
        sigma_high = sigma_guess
        p_high = p_guess
        j=1
        p_low = p_guess # dummy intialization
        while p_low > mrktprice:
            sigma_low = 2**(-j*sigma_guess)
            opt.sigma = sigma_low
            p_low = opt.usaPutPrice()
            j = j + 1 

    while sigma_high-sigma_low > threshold: #if the interval length is too big
        sigma_mid = 0.5*(sigma_high+sigma_low)
        opt.sigma = sigma_mid
        p_mid = opt.usaPutPrice()
        if p_mid < mrktprice:
            sigma_low = sigma_mid
        else:
            sigma_high = sigma_mid
    return round(0.5*(sigma_high+sigma_low), precision)

def implVol_call(contract, mrktprice, precision=5):
    '''Return implied volatily of a call option'''
    threshold = 10**-(precision+1)
    opt = contract
    sigma_guess = contract.sigma # place the guess of sigma in the contract object
    p_guess = contract.callPrice()
    if p_guess < mrktprice:
        sigma_low = sigma_guess
        p_low = p_guess
        j = 1
        p_high = p_guess # dummy intialization
        while p_high < mrktprice:
            sigma_high = 2**j*sigma_guess
            opt.sigma = sigma_high
            p_high = opt.callPrice()
            j = j + 1 
    else:
        sigma_high = sigma_guess
        p_high = p_guess
        j=1
        p_low = p_guess # dummy intialization
        while p_low > mrktprice:
            sigma_low = 2**(-j*sigma_guess)
            opt.sigma = sigma_low
            p_low = opt.callPrice()
            j = j + 1 

    while sigma_high-sigma_low > threshold: #if the interval length is too big
        sigma_mid = 0.5*(sigma_high+sigma_low)
        opt.sigma = sigma_mid
        p_mid = opt.callPrice()
        if p_mid < mrktprice:
            sigma_low = sigma_mid
        else:
            sigma_high = sigma_mid
    return round(0.5*(sigma_high+sigma_low), precision)
    
S_0 = 50
T = 0.25
r = 0.02
K = 50 #Strike Price
price = 2.75
n = 1000

x = option(S_0, T, r, K, n, .21)
print(implVol_usaput(x,price))
#print(implVol_europut(x,price))
#print(implVol_call(x,price))
#y = x
#y.sigma = implVol_call(x,price)
#print(y.callPrice())


