from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.matrix.all import *
from sage.functions.generalized import kronecker_delta
from sage.misc.flatten import flatten
from sage.arith.functions import lcm
from sage.arith.misc import GCD
from sage.combinat.root_system.dynkin_diagram import DynkinDiagram
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.misc.cachefunc import cached_method
from functools import reduce

__all__=['RSG','SG','TamelyLaced','UntwistedAffine','Unknown','LengthOne']

class TDatumConstructor(object):
    def __init__(self,variable_name='z'):
        self._variable_name =variable_name
        
    def variable_name(self):
        return self._variable_name
    
    def indices(self):
        raise NotImplementedError()
    
    def size(self):
        return len(self.indices())
        
    def _N0_entry(u,v):
        raise NotImplementedError()
    
    def _N_plus_entry(self,u,v):
        raise NotImplementedError()
    
    def _N_minus_entry(self,u,v):
        raise NotImplementedError()
    
    def variable(self):
        return LaurentPolynomialRing(QQ,self.variable_name()).gen()
    
    def t_datum(self):
        N0=matrix(self.size() , lambda u,v: self._N0_entry(self.indices()[u],self.indices()[v]))
        Ap = N0-matrix(self.size() , lambda u,v: self._N_plus_entry(self.indices()[u],self.indices()[v]))
        Am = N0-matrix(self.size() , lambda u,v: self._N_minus_entry(self.indices()[u],self.indices()[v]))
        return (Ap,Am)
    
        
class RSG(TDatumConstructor):
    def __init__(self,n_list,variable_name='z'):
        super(RSG,self).__init__(variable_name)
        self._n_list = n_list
        
    def n(self,a):
        return self._n_list[a-1]
    def n_tilde(self,a):
        return self.n(a) - 2*kronecker_delta(a,1)
    def F(self):
        return len(self._n_list)
    def p(self,a):
        if a==1:
            return 1
        elif a==2:
            return self.n(1)
        else:
            return self.q(a-1)
    def q(self,a):
        if a==0:
            return 1
        elif a==1:
            return self.n(1)
        else:
            return self.n(a)*self.q(a-1)+ self.p(a-1)
        
    def epsilon(self,a):
        return (-1)**(a-1)
    
    @cached_method
    def indices(self):
        am = [ range(1,self.n(1)-1) ] + [ range(1,1+self.n(a)) for a in range(2,1+self.F()) ]
        return reduce(lambda a,b: a+b,[[ (a,m) for m in am[a-1] ] for a in range(1,1+self.F())])
    
    def _N0_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(b,k):
            result += z**(2*self.p(a)) + 1
        return result
    
    def _N_plus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if self.n(1)==2 and (a,m)==(2,1) and (b,k)==(2,1):
            result+=z**(self.p(a))
        if m==1 and self.epsilon(a)==-1:
            if (b,k) == (a-2 , self.n_tilde(a-2)):
                result +=  z**(self.p(a))
            elif b==a-1:
                result += z**(self.p(a))*(z**(self.p(a) -(self.n_tilde(a-1) +1 -k)*self.p(a-1) ) + z**(-self.p(a) +(self.n_tilde(a-1) +1 -k)*self.p(a-1)))
        if m==self.n_tilde(a) and self.epsilon(a)==1:
            if (b,k)==(a+1,1):
                result+= z**self.p(a)
        if abs(m - k)==1 and a == b and self.epsilon(a)==-1:
            result+= z**self.p(a)     
        return result
    
    def _N_minus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(2,1) and (b,k) == (1,1):
            result+= z**self.p(a)
        if m==1 and self.epsilon(a)==1:
            if (b,k) == (a-2 , self.n_tilde(a-2)):
                result+= z**self.p(a)
            elif b==a-1:
                result += z**self.p(a)*(z**(self.p(a) -(self.n_tilde(a-1) +1 -k)*self.p(a-1) ) + z**(-self.p(a) +(self.n_tilde(a-1) +1 -k)*self.p(a-1)))
        if m==self.n_tilde(a) and self.epsilon(a)==-1:
            if (b,k)==(a+1,1):
                result+= z**self.p(a)
        if abs(m - k)==1 and a == b and self.epsilon(a)==1:
            result+= z**self.p(a)               
        return result
            
    
class SG(TDatumConstructor):
    def __init__(self,n_list,variable_name='z'):
        super(SG,self).__init__(variable_name)
        self._n_list = n_list
        
    def n(self,a):
        return self._n_list[a-1]
    def F(self):
        return len(self._n_list)
    def p(self,a):
        if a==1:
            return 1
        elif a==2:
            return self.n(1)
        else:
            return self.q(a-1)
    def q(self,a):
        if a==0:
            return 1
        elif a==1:
            return self.n(1)
        else:
            return self.n(a)*self.q(a-1)+ self.p(a-1)
        
    def epsilon(self,a):
        return (-1)**(a-1)
    
    @cached_method
    def indices(self):
        am = [[-2,-1,0]+list(range(1,self.n(1)-1))] + [ list(range(1,1+self.n(a))) for a in range(2,1+self.F())]
        return reduce(lambda a,b: a+b,[[ (a,m) for m in am[a-1] ] for a in range(1,1+self.F())])
    
    
    def _N0_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(b,k):
            result += z**(2*self.p(a)) + 1
        return result
    
    def _N_plus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(2,1):
            if (b,k) == (2,2) or (b,k) == (1,-1) or (b,k) == (1,-2):
                result+= z**(self.p(a))                
            elif b == 1:
                result+=  z**self.p(a)*(z**(k+1) + z**(-k-1))
        elif a>2 and m==1 and self.epsilon(a)==-1:
            if (b,k) ==(a,2):
                result +=  z**(self.p(a))
            elif (b,k) == (a-2 , self.n(a-2) - 2*kronecker_delta(a,3)):
                result +=  z**(self.p(a))
            elif b==a-1:
                result +=  z**(self.p(a))*(z**(self.p(a) -(self.n(a-1) +1 -k)*self.p(a-1) ) + z**(-self.p(a) +(self.n(a-1) +1 -k)*self.p(a-1) ))
        elif a>2 and m==1 and self.epsilon(a)==1:
            if self.n(a)==1 and (b,k)==(a+1,1) :
                result+= z**self.p(a)
        elif abs(self.indices().index((a,m)) - self.indices().index((b,k)))==1:
            if m>=0 and k>=0:
                if b%2==0:
                    result+= z**self.p(a)   
        return result
    
    def _N_minus_entry(self,u,v):
        z=self.variable()
        a,m=u
        b,k=v
        result=0
        if (a,m)==(2,1):
            if self.n(2)==1 and (b,k)==(3,1):
                result+= z**self.p(a)
        elif a>2 and m==1 and self.epsilon(a)==1:
            if (b,k) ==(a,2):
                result+= z**self.p(a)
            elif (b,k) == (a-2 , self.n(a-2) - 2*kronecker_delta(a,3)):
                result+= z**self.p(a)
            elif b==a-1:
                result+= z**self.p(a)*(z**(self.p(a) -(self.n(a-1) +1 -k)*self.p(a-1) ) + z**(-self.p(a) +(self.n(a-1) +1 -k)*self.p(a-1) ))
        elif a>2 and m==1 and self.epsilon(a)==-1:
            if self.n(a)==1 and (b,k)==(a+1,1) :
                result+= z**self.p(a)
        elif (m,k)==(-1,0) or (m,k)==(0,-1) or (m,k)==(-2,0) or (m,k)==(0,-2):
            result+=z**self.p(a)
        elif abs(self.indices().index((a,m)) - self.indices().index((b,k)))==1:
            if m>=0 and k>=0:
                if b%2==1:
                    result+= z**self.p(a)       
        return result
    
    