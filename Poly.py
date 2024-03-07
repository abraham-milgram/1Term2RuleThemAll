import matplotlib.pyplot as plt

class polynomial:
    def __init__(self, coefs):
        """
        A polynomial is defined by a list of coefficients, with the
        highest power of the polynomial in the list at the start. For example,
        the polynomial x^3 + 2x^2 + 3x + 4 is defined by the list [1, 2, 3, 4].
        """
        self.coefs = coefs
        self.degree = len(coefs) - 1
        self.delete_zeros()

    # arithmetic operators for polynomials
    def __add__(self, other):
        # making the polynomials the same size
        polylens = sorted([self, other], key = lambda x: len(x.coefs))
        polylens[0].coefs = [0]*(len(polylens[-1].coefs) - len(polylens[0].coefs)) + polylens[0].coefs
        return polynomial([polylens[0].coefs[i]+polylens[1].coefs[i] for i in range(len(polylens[0].coefs))])

    def __sub__(self, other):
        return self + (-other)

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

    def __truediv__(self, divisor):
        result = []
        while 1:
            # dividing the polynomial by the divisor using the polynomial long division method
            result.append((self.coefs[0] / divisor.coefs[0], self.degree - divisor.degree))
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

    # function to find the definite integral of a polynomial
    def definite_integral(self, ll, ul):
        return self.antiderivative().plugin(ul) - self.antiderivative().plugin(ll)

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

    # function to plot a polynomial
    def plot(self, x_range):
        x = [i for i in range(x_range[0], x_range[1])]
        y = [self.plugin(i) for i in x]
        plt.plot(x, y)
        plt.show()
