#################################################################
# Author: Debasis Das 
# This code solves sLLG for a PMA nanomagnet
# The effective field consists of : 
# Uniaxial field, Demag field, and Magnetoelectric field.
#
# The parameters are taken from https://doi.org/10.1038/s41598-017-11732-w 
#
# To solve the sLLG Huen method is used. 
# Details of the solution method can be found in 
# https://arxiv.org/pdf/1607.04596v4.pdf
#################################################################

import numpy as np
from random import gauss
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 

plt.rcParams.update({'font.size':24})
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams.update({'font.weight':'bold'})
plt.rcParams["font.family"] = "Times New Roman"
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
#################### Real unit  LLG ###########################
def deterministic_LLG(m,h):
	precision=-np.cross(m,h)
	damping=-(alpha)*np.cross(m,np.cross(m,h))
	
	dmdt=precision+damping
	return dmdt
#################################################################
def stochastic_LLG(m,h_T):
	s_precision=-np.cross(m,h_T)
	s_damping=-(alpha)*np.cross(m,np.cross(m,h_T))
	
	s_dmdt=s_precision+s_damping
	return s_dmdt
#################################################################
def dimless_determinstic_field_cal(m_vec):
	global hk, hd, h_ME_val
	h_uni=hk*np.array([0,0,m_vec[2]])
	h_demag=np.array([hd[0]*m_vec[0], hd[1]*m_vec[1], hd[2]*m_vec[2]])
	h_ME=h_ME_val*np.array([0,0,1])
	
	h_eff_det=h_uni+h_demag+h_ME
	return h_eff_det
#################################################################
def get_random_normal_vector():
	vec=np.random.normal(0,1,3)
	norm_vec=normalize(vec)
	return norm_vec
##################### Vector magnitude #############################
def normalize(M):
	magnitude = np.sqrt((M[0])**2+(M[1])**2+(M[2])**2)
	norm_m=M/magnitude
	return norm_m
#################################################################

#################### constant parameters ############################
gamma=1.76e11;           # Gyromagnetic ratio [(rad)/(s.T)]
mu0=4*np.pi*1e-7 ;      # in T.m/A

q=1.6e-19;               # in Coulomb
hbar=1.054e-34;          # Reduced Planck's constant (J-s)
K_B=1.38064852e-23    #in J/K
#################### parameters related to nanomagnet ################
alpha=0.15              # Gilbert damping parameter
Ms=250*1e3            # in A/m
Length=16e-9
Width=8e-9
t_FL=0.9e-9
A_MTJ=Length*Width 
Area=A_MTJ
V=A_MTJ*t_FL;             # Volume [m^3]
Ki=0.068*1e-3
Ku2=Ki/t_FL              # in J/m^3
print('Ku2 = ' + str(Ku2))
################# Anisotropy field ###################################
Hk=2*Ku2/(mu0*Ms)             # Uniaxial field [in A/m]
hk=Hk/Ms                                # Normalized uniaxial field
print('Uniaxial = ' + str(hk))
################## Demagnetization field related #######################
Nxx=0.0
Nyy=0
Nzz=1.0

hdx=-Nxx*Ms	# in A/m
hdy=-Nyy*Ms	# in A/m
hdz=-Nzz*Ms		# in A/m
Hd=np.array([hdx, hdy, hdz])
hd=Hd/Ms              #normalized demag field
################## Magneto-electric field related #######################
alpha_ME=0.03*1e-7
tox_ME=5e-9

Vpos=np.linspace(0.5,4,15) 
Vneg=np.linspace(-4,-0.5,15)  
Voltage=[]
Voltage=np.append(Voltage,Vneg)   
Voltage=np.append(Voltage,Vpos)

print('Voltage = ' + str(Voltage))	
#exit()



###################### time related portion ###########################
dim_less_time_fac = gamma*mu0*Ms

stop_time=2.0e-9             # in s
stop_tau=stop_time*dim_less_time_fac
# delta_t should be modified according to the desired accuracy
delta_t=2e-12
delta_tau = delta_t*dim_less_time_fac

n=int(stop_tau/delta_tau)
t=np.linspace(0,stop_time,n) 
tau=np.linspace(0,stop_tau,n) # in ns
delta_tau=tau[1]-tau[0]
##################### Themal field info #############################
T_K=300     # in Kelvin
const=(2*alpha*K_B*T_K)/(gamma*Ms*V*delta_t)  
H_T=(1/mu0)*np.sqrt(const)
H_T_norm=(H_T/Ms)

avg_N_mz=np.zeros(len(Voltage))

for iV in range(len(Voltage)):
	V_IN=Voltage[iV]
	
	H_ME_coeff=alpha_ME*(V_IN/tox_ME)  # in Tesla
	H_ME=H_ME_coeff/mu0         # in A/m
	h_ME_val=H_ME/Ms            # Dimensionless Magneto-electric field constant

	loop=n-1
	N_runs=30
	avg_mz=np.zeros(N_runs)
	
	for iN in range(N_runs):
		################## unit Magnetization vector ##########################
		m=np.zeros((n,3))            
		#  Initial magnetization
		if V_IN>=0:
			mz0=-1.0
		else:
			mz0=1.0	
		my0=np.sqrt(1-mz0**2)
		mx0=0
		m[0,:]=[mx0,my0,mz0]
		####################
		eta_n=np.sqrt(delta_tau)

		g_Xn_tn=np.zeros(3)
		g_X_var_nplus1_tnplus1=np.zeros(3)
		for i in range(loop):
			print('---------------------------------------------------------')
			print('V= %f; iN = %d;  n = %d ; i = %d' %(V_IN,iN,n,i))
			X_n=m[i,:]
			h_eff_det_i=[]
			h_eff_det_i =dimless_determinstic_field_cal(X_n)
			f_Xn_tn=(deterministic_LLG(X_n,h_eff_det_i))
			G01=get_random_normal_vector()
			#print('H_T_Norm = ' + str(H_T_norm))
			h_T_n=H_T_norm*G01
			if T_K>0:
				g_Xn_tn=(stochastic_LLG(X_n,h_T_n))

			X_var_nplus1=normalize(X_n+f_Xn_tn*delta_tau+g_Xn_tn*eta_n*G01)
			h_eff_det_m_var=dimless_determinstic_field_cal(X_var_nplus1)
			f_X_var_nplus1_tnplus=(deterministic_LLG(X_var_nplus1,h_eff_det_m_var))
			Deterministic_update = (1/2.0)*(f_X_var_nplus1_tnplus+f_Xn_tn)*delta_tau
			G01_nplus1=get_random_normal_vector()
			h_T_nplus1=H_T_norm*G01_nplus1
			if T_K>0:
				g_X_var_nplus1_tnplus1=(stochastic_LLG(X_var_nplus1,h_T_nplus1))

			Stochastic_update=(1/2.0)*(g_X_var_nplus1_tnplus1+g_Xn_tn)*eta_n*G01
			X_n_plus1=normalize(X_n+Deterministic_update+Stochastic_update)
			m[i+1,:]=X_n_plus1

		avg_mz[iN]=np.average(m[:,2])
		
	avg_N_mz[iV]=np.average(avg_mz)
	
#print('Voltage = ' + str(Voltage))	
#print('Average mz = ' + str(avg_N_mz))

fig = plt.figure(figsize=(14,7))
lw=2.2
plt.plot(Voltage, avg_N_mz, '*-', markersize=15, linewidth=lw)
plt.ylim([-1.2,1.2])
plt.yticks([-1.0,0,1.0])
plt.xlabel('Voltage (V)')
plt.ylabel(r"$<m_z>$")
plt.grid()
plt.savefig('Average_mz_voltage_variation_plot.pdf', bbox_inches='tight', pad_inches=0.2)
plt.show()

'''
t=t*1e9 # time converted to ns scale

fig = plt.figure(figsize=(14,7))
lw=2.2

#plt.plot(t[0:loop],m[0:loop,0], 'b', linewidth=lw, label='mx')
#plt.plot(t[0:loop],m[0:loop,1], 'k', linewidth=lw, label='my')
plt.plot(t[0:loop],m[0:loop,2], 'r', linewidth=lw, label='mz')
plt.title('V_IN = ' +str(V_IN))
plt.grid()
plt.legend()
plt.xlabel('Time(ns)')
plt.ylabel(r"$m$")
plt.ylim([-1.2,1.2])
plt.yticks([-1.0,0,1.0])
#plt.savefig('Stochastic_Huen_pbit_result_V_IN_pos_4.pdf', bbox_inches='tight', pad_inches=0.2)
plt.show()
'''
