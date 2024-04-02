########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Dorival M. Pedroso                                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

from numpy import matrix, linalg, sqrt, zeros

# Vector
# ======
class Vector(matrix):
    # Constructor
    # ===========
    def __init__(self, vals, dtype=float, copy=False):
        super(Vector,self).__init__()
        if self.shape[1]>self.shape[0]:
            self.resize((self.shape[1],self.shape[0]))

    # Access item
    # ===========
    def __getitem__(self, key):
        return matrix.__getitem__(self, (key,0))

    # Euclidian norm
    # ==============
    def norm(self):
        return sqrt((self.T*self)[0])

    # Nice Print
    # ==========
    def write(self, nf='%10g', Tol=1.0e-14):
        m = self.shape[0] # number of rows
        lin = ''          # empty string
        for i in range(m):
            if abs(self[i])<Tol: lin += nf % 0
            else:                lin += nf % self[i]
            lin += '\n'
        print lin,


# Matrix
# ======
class Matrix(matrix):
    # Determinant
    # ===========
    def det(self): return linalg.det(self)

    # Access item
    # ===========
    def __getitem__(self, key):
        if isinstance(key, int):
            if self.shape[1] == 1: # column matrix
                return matrix.__getitem__(self, (key,0))
            else: return matrix.__getitem__(self, key)
        else: return matrix.__getitem__(self, key)

    # Nice Print
    # ==========
    def write(self, nf='%10g', Tol=1.0e-14):
        m = self.shape[0] # number of rows
        n = self.shape[1] # number of columns
        lin = ''          # empty string
        for i in range(m):
            for j in range(n):
                if abs(self[i,j])<Tol: lin += nf % 0
                else:                  lin += nf % self[i,j]
            lin += '\n'
        print lin,

    # Nice Print (complex numbers)
    # ============================
    def writec(self, nf='%10g', Tol=1.0e-14):
        m = self.shape[0] # number of rows
        n = self.shape[1] # number of columns
        lin = ''          # empty string
        nf = '('+nf+' + '+nf+'j) '
        for i in range(m):
            for j in range(n):
                if abs(self[i,j].real)<Tol and abs(self[i,j].imag)<Tol:
                    lin += nf % (0,0)
                else:
                    lin += nf % (self[i,j].real,self[i,j].imag)
            lin += '\n'
        print lin,


# Create vector with zeros
# ========================
def ZeroVec(m): return Vector(zeros((m,1)))


# Create matrix with zeros
# ========================
def ZeroMat(m,n): return Matrix(zeros((m,n)))


# Vector dot product
# ==================
def Dot (U,V): return (U.T*V)[0]


# Vector dyadic product
# =====================
def Dyad(U,V):
    m = U.shape[0]
    n = V.shape[0]
    M = Matrix(zeros(shape=(m,n)))
    for i in range(m):
        for j in range(n):
            M[i,j] = U[i]*V[j]
    return M


# test
# ====
if __name__=="__main__":

    K = Matrix([[1., 2., 3., 4.],
                [5., 6., 7., 8.]])

    v = Vector([1., 2., 3., 4.])

    C = Matrix([[1.],
                [2.],
                [3.]])

    vdyv = Dyad(v,v)

    print 'K ='; K.write()
    print 'v ='; v.write()
    print 'v dyad v ='; vdyv.write()
    print 'C ='; C.write()
    print 'C[0] =', C[0]
    print 'C[1] =', C[1]
    print 'C[2,0] =', C[2,0]
