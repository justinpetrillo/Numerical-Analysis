import math
from itertools import product

def f1(x):
	y = x**4 - 4*x**3 + 8*x**2 - 2*x - 4
	return y
f_display='f(x) = x^4 - 4x^3 + 8x^2 - 2x - 4'

#Tol=[p_diff,p_diff/p,f(p)]
def binary_search(a,b,Tol,N):
	print('finding the zero of F(x)=',f_display)
	print('between',a,b)
	print('Tol=[p_max_diff,p_max_diff/p_max,f(p_max)]', Tol)
	print('N=',N)
	Xa=[0 for x in range(N)]
	Xb=[0 for x in range(N)]
	Xp=[0 for x in range(N)]

	Xa[0]=a
	Xb[0]=b

	found = 0

	for i in range(N):
		Xp[i] = float((Xa[i]+Xb[i])/2)
		print(Xa[i],Xp[i],Xb[i])

		f_error=abs(f(Xp[i]))
		p_diff_error=abs(Xp[i]-Xp[i-1])
		p_prop_error=p_diff_error/abs(Xp[i])

		if f_error<Tol[2]:
			found = 1
		if p_diff_error<Tol[0]:
			found = 1
		if p_prop_error<Tol[1]:
			found = 1

		if found == 1:
			p = Xp[i]
			print('selection:p=',p)
			print('errors(abs):')
			print('f-0:',f_error)
			print('last p differences:',p_diff_error)
			print('p seq difference / max_p :',p_prop_error)
			return p
		else:
			if f(Xp[i])*f(Xa[i])<0:
				Xa[i+1]=Xa[i]
				Xb[i+1]=Xp[i]	
			if f(Xp[i])*f(Xb[i])<0:
				Xa[i+1]=Xp[i]
				Xb[i+1]=Xb[i]
	if found == 0:
		print("Tol", Tol, "not reached in", N, "max iterations")
		return

#binary_search(1.8,2.3,[.1,.01,.1],100)

def g(x):
	y = f(x)
	return x - y

def fixed_point(x0,Tol,N):
	print('Finding the fixed point of G(x) = x - F(x)=','x - ',f_display,'from x0=',x0,'Tol', Tol, 'N=',N)
	X0=[0 for x in range(N+1)]
	X0[0]=x0

	found = 0

	for i in range(N):
		X0[i+1] = float(g(X0[i]))
		print(X0[i+1])
		seq_diff_error = abs(g(X0[i+1])-X0[i])
		if seq_diff_error<Tol:
			p = X0[i+1]
			found = 1
			print('selection:p=',p)
			print('sequence difference error=',seq_diff_error)
			return p
	if found == 0:
		print("Tol", Tol, "not reached in", N, "max iterations")
		return

#fixed_point(1.7,.1,100)

def df_dx(x):
	df_dy = (1/4)*x**3 - (1/3)*4*x**2 + (1/2)*8*x - 2
	return df_dy
df_display='df/dx = x^3/4 - 4x^2/3+ 4x - 2'

def Newton_Method(x0,Tol,N):
	print('Finding the zero of F(x)=',f_display,df_display,'from x0=',x0,'Tol', Tol, 'N=',N)
	X0=[0 for x in range(N+1)]
	X0[0]=x0

	found = 0

	for i in range(N):
		X0[i+1] = float(X0[i]-f(X0[i])/df(X0[i]))
		print(X0[i+1])
		seq_diff_error = abs(g(X0[i+1])-X0[i])
		if seq_diff_error<Tol:
			p = X0[i+1]
			found = 1
			print('selection:p=',p)
			print('sequence difference error=',seq_diff_error)
			return p
	if found == 0:
		print("Tol", Tol, "not reached in", N, "max iterations")
		return

#Newton_Method(1.7,.1,100)

def Secant_Method():
	return 0
 
def interp():
	n=int(input("degree: n=?"))
	x = []
	for i in range(n+1):
		x = x + [0]
	y = []
	for i in range(n+1):
		y = y + [0]
	yprime=[]
	for i in range(n+1):
		yprime=yprime+[0]
	print ('enter all values of x,y,yprime')
	for i in range(n+1):
		print('i=',i)
		x[i]=float(input('x=?'))
		y[i]=float(input('y=?'))
		yprime[i]=float(input('yprime=?'))
	print()
	print('the data points are:')
	for i in range(n+1):
		print("(",x[i],",",y[i],")")
	print()
	def f(w):
		return y[w]
	def ddf(t):
		if x[t]==x[t+1]:
			b=yprime
		else:
			b=(f(t+1)-f(t))/(x[t+1]-x[t])
		return b
	def sdf(u):
		c=(ddf(u+1)-ddf(u))/(x[u+2]-x[u])
		return x
	def tdf(v):
		e=(sdf(v+1)-sdf(v))/(x[v+3]-x[v])
		return e
	print('The Divided Difference Values:')
	for i in range(n):
		print (ddf(i), ' ')
	print()
	print('The 2nd Divided Difference Values:')
	for i in range(n-1):
		print (sdf(i), ' ')	
	print()
	print('The 3rd Divided Difference Values:')
	for i in range(n-2):
		print (tdf(i), ' ')	
	def L(p):
		poly=f(0)+ddf(0)*(p-x[0])+sdf(0)*(p-x[0])*(p-x[1])
		poly=poly+tdf(0)*(p-x[0])*(p-x[1])*(p-x[2])
		return poly
	return L(p)


def Taylor_Approximation(x0,f):
	return 0

def product( iterable ):
    a=1
    for n in iterable:
        a *= n
    return a

#x,f are n+1 range lists
#zeta is f^(n+1)(Zeta); use NA for Not-available
def Lagrange_Interpolation(x,f,X,zeta):
	n = len(x) - 1
	L=[0 for k in range(n+1)]
	for k in range(n+1):
		L[k]=product([(X-x[i])/(x[k]-x[i]) for i in range(n+1) if i!=k])
	
	P=sum([f[k]*L[k] for k in range(n+1)])
	print(P)
	if zeta != 'NA':
		error = (zeta/math.factorial(n+1))*product([(X-x[i]) for i in range(n+1)])
		print('error=',error)


	print('lagrange polynomial:')
	for i in range(n+1):
		for j in range(n+1):
			if i!=j:
				print('(X-',x[j],')/(',x[i],'-',x[j],')',end='')
				if j<n:
					print('*',end='')
		if i<n:
			print('+',end='')
	print('')
	print('at X=',X)
	for i in range(n+1):
		for j in range(n+1):
			if i!=j:
				print('(',X,'-',x[j],')/(',x[i],'-',x[j],')',end='')
			if j<n:
				print('*',end='')
		if i<n:
			print('+',end='')
	print('')
	print(P)

	return P
#Lagrange_Interpolation([1,2,3],[4,5,6], 3.5,'NA')


	#m is a list of input registers
def Neviles_Method(x,f,X,m):
	#m is a list of input registers
	xm=[x[i] for i in m]
	P=Lagrange_Interpolation(xm,f,X,NA)




def f2(x,y):
	z=1/(x^2 + y^2)
	return z

'''
#x,f are n+1 range lists
def Newton_DivDiff_Method(x,f,X):
	n = len(x)
	F=[[i for i in range(n)] for j in range(n)]
	for i in range(n):
		for j in range(i):
			if i>0 and j>0:
				F[i][j]=(F[i][j-l] - F[i-1][j-1])/(x[i] - x[i-1])
			else:
				F[0][0]=f[0]
				F[1][0]=(f[1]-f[0])/(x[1]-x[0])
	def P(X):
		P=sum([F[i][i]*product([(X - x[j]) for j<i]) for  i in n])
		return P
	return P(X)
'''


