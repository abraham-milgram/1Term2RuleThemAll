from copy import deepcopy
import matplotlib.pyplot as plt
import math


class polynomial:
    def __init__(self, coefs):
        '''
        coefs is a list of coeficients where the list position equals 
        '''

        self.coefs = coefs
        self.degree = len(coefs) - 1
        self.delete_zeros()

    # arithmetic operators for polynomials

    def __add__(self, other):
        # Converts the other term to a polynomial if not already
        if not isinstance(other, polynomial):
            other = polynomial([other])
        biggerPoly = self.coefs if self.degree > other.degree else other.coefs
        smallerPoly = other.coefs if self.degree >= other.degree else self.coefs
        for i in range(len(biggerPoly) - len(smallerPoly)):
            smallerPoly.insert(0, 0)
        return polynomial([biggerPoly[i]+smallerPoly[i] for i in range(len(biggerPoly))])

    def __sub__(self, other):
        # add leading zeros to other to make it the same length as self
        if len(self.coefs) > len(other.coefs):
            for i in range(len(self.coefs) - len(other.coefs)):
                other.coefs.insert(0, 0)
        elif len(self.coefs) < len(other.coefs):
            for i in range(len(other.coefs) - len(self.coefs)):
                self.coefs.insert(0, 0)
        return polynomial([self.coefs[i] - other.coefs[i] for i in range(len(self.coefs))])

    def __mul__(self, other):
        result = [0 for i in range(len(self.coefs) + len(other.coefs))]
        for i in range(len(self.coefs)):
            for j in range(len(other.coefs)):
                result[i+j] += self.coefs[i] * other.coefs[j]
        return polynomial(result[:-1])

    def __pow__(self, power):
        if power == 0:
            return polynomial([1])
        else:
            return self * (self ** (power - 1))

    def __str__(self):
        return str(self.coefs)

    def __repr__(self):
        return str(self.coefs)

    def __eq__(self, other):
        return self.coefs == other.coefs

    def __ne__(self, other):
        return not self == other
    # function to divide a polynomial by a constant

    def __truediv__(self, divisor):
        result = []
        while 1:
            result.append(
                (self.coefs[0] / divisor.coefs[0], self.degree - divisor.degree))
            subby = divisor * polynomial([result[-1][0]] + [0] * result[-1][1])
            self = self - subby
            if self.degree < divisor.degree:
                break

        return polynomial([i[0] for i in result] + [0] * result[-1][1]), self, divisor

    def __mod__(self, divisor):
        return (self/divisor)[1:]

    def __floordiv__(self, divisor):
        return (self/divisor)[0]

    def plugin(self, x):
        return sum([self.coefs[i] * x**(self.degree - i) for i in range(self.degree + 1)])

    # function to delete the leading zeros in the polynomial
    def delete_zeros(self):
        if self.coefs != []:
            while self.coefs[0] == 0:
                self.coefs.pop(0)
                self.degree -= 1
                if self.coefs == []:
                    break

    # function to find the derivative of a polynomial
    def derivative(self):
        return polynomial([self.coefs[i] * (self.degree - i) for i in range(0, self.degree)])
    # function to find the antiderivative of a polynomial
    def antiderivative(self):
        return polynomial([self.coefs[i] / (self.degree - i + 1) for i in range(0, self.degree)] + [0.0])

    # function to find the riemann sum of a polynomial
    def riemann_sum(self, ll, ul, n, t):
        dx = (ul - ll) / n
        rtypes = {'left': (lambda x: 1, slice(0, -1)), 'right': (lambda x: 1, slice(1, n + 1)), 'trapezoid': (lambda x: 1 if x in (0,n) else 2, slice(0,n+1))}
        if t != 'midpoint':
            return sum([rtypes[t][0](i)*self.plugin(ll + i * dx) for i in range(n + 1)][rtypes[t][1]]) * dx/rtypes[t][0](n+1)
        else:
            return sum([self.plugin(((ll + i * dx) + (ll + (i + 1) * dx))/2) for i in range(n)]) * dx
    
    # function to find the slope of a polynomial at a point
    def slope(self, x):
        return self.derivative().plugin(x)

    # function to find the tangent of a polynomial at a point
    def tangent(self, x):
        return polynomial([self.slope(x), self.plugin(x) - self.slope(x) * x])

    # function to evaluate the sign of a number
    def sign(self, x):
        if x > 0:
            return 1
        elif x < 0:
            return -1
        else:
            return 0

        
    # function to plot a polynomial

    def plot(self, x_range):
        x = [i for i in range(x_range[0], x_range[1])]
        y = [self.plugin(i) for i in x]
        plt.plot(x, y)
        plt.show()

