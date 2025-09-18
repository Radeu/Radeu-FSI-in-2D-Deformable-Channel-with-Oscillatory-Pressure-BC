
from dolfin import *
from ufl import indices, Jacobian, Min, Max, shape
from mshr import *
import numpy as np
import os


####### Domain and mesh setup #######

# Parameters defining domain geometry:
SOLID_TOP = 9e-5
SOLID_RIGHT = 1.0256e-3
SOLID_BOTTOM = 4.5e-5
SOLID_LEFT = 0.0
OMEGA_H = 9e-5
OMEGA_W = SOLID_RIGHT
REF_MARGIN = 0.1

# Define the mshr geometrical primitives for this domain:
r_Omega = Rectangle(Point(0,0),Point(OMEGA_W,OMEGA_H))
r_Omega_s = Rectangle(Point(SOLID_LEFT,SOLID_BOTTOM),
                      Point(SOLID_RIGHT,SOLID_TOP))

# Enumerate subdomain markers
FLUID_FLAG = 0
SOLID_FLAG = 1
# Zero is the default flag, and does not need to be
# explicitly set for the fluid subdomain.
r_Omega.set_subdomain(SOLID_FLAG,r_Omega_s)

# Parameters defining refinement level:
N = 600

# Generate mesh of Omega, which will have a fitted
# subdomain for Omega_s.
mesh = generate_mesh(r_Omega, N)

# Mesh-derived quantities:
d = mesh.geometry().dim()
n_y = FacetNormal(mesh)
I = Identity(d)

# Define subdomains for use in boundary condition definitions:
class Wall_b(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary 
                and (x[1] < DOLFIN_EPS))

class Wall_t(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary 
                and (x[1] > SOLID_TOP - DOLFIN_EPS))        

class Wall_l(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary 
                and (x[0] < DOLFIN_EPS)
				and x[1] > SOLID_BOTTOM - DOLFIN_EPS)

class Wall_r(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary 
                and (x[0] > SOLID_RIGHT - DOLFIN_EPS)
				and x[1] > SOLID_BOTTOM - DOLFIN_EPS)				

class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary 
                and (x[0] < DOLFIN_EPS)
				and x[1] < SOLID_BOTTOM + DOLFIN_EPS)

#class Inflow(SubDomain):
   # def inside(self, x, on_boundary):
		#return (on_boundary and x[0] < DOLFIN_EPS and x[1] < SOLID_BOTTOM + DOLFIN_EPS)

class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary and (x[0] > (OMEGA_W - DOLFIN_EPS))
				and x[1] < SOLID_BOTTOM + DOLFIN_EPS)
		
class SolidDomainClosure(SubDomain):
    def inside(self, x, on_boundary):
        return (x[1] > SOLID_BOTTOM - DOLFIN_EPS
                and x[0] > SOLID_LEFT - DOLFIN_EPS
                and x[0] < SOLID_RIGHT + DOLFIN_EPS)
class SolidDomainInterior(SubDomain):
    def inside(self, x, on_boundary):
        
        return (x[0] > SOLID_LEFT - DOLFIN_EPS
                and x[0] < SOLID_RIGHT + DOLFIN_EPS
                and x[1] < SOLID_TOP + DOLFIN_EPS
				and x[1] > SOLID_BOTTOM + DOLFIN_EPS)
 
# Set up integration measures, with flags to integrate over
# subsets of the domain.
markers = MeshFunction('size_t', mesh, d, mesh.domains())
bdry_markers = MeshFunction('size_t', mesh, d-1, 0)
OUTFLOW_FLAG = 2
INFLOW_FLAG = 3
Outflow().mark(bdry_markers,OUTFLOW_FLAG)
Inflow().mark(bdry_markers,INFLOW_FLAG)
dx = dx(metadata={'quadrature_degree': 2},
        subdomain_data=markers)
ds = ds(metadata={'quadrature_degree': 2},
        subdomain_data=bdry_markers)

####### Elements and function spaces #######

# Define function spaces (equal order interpolation):
cell = mesh.ufl_cell()
Ve = VectorElement("CG", cell, 1)
Qe = FiniteElement("CG", cell, 1)
VQe = MixedElement((Ve,Qe))
# Mixed function space for velocity and pressure:
W = FunctionSpace(mesh,VQe)
# Function space for mesh displacement field, 
# which will be solved for separately in a 
# quasi-direct scheme:
V = FunctionSpace(mesh,Ve)

####### Set up time integration variables #######

TIME_INTERVAL = 0.084*2
N_STEPS = 84000*2
Dt = Constant(TIME_INTERVAL/N_STEPS)

# Mesh motion functions:
uhat = Function(V)
uhat_old = Function(V)
du = TestFunction(V)
vhat = (uhat-uhat_old)/Dt

# Fluid--structure functions:
(dv, dp) = TestFunctions(W)
w = Function(W)
v,p = split(w)
p_mean = Function(W)
w_old = Function(W)
v_old, p_old = split(w_old)
dv_dr = (v - v_old)/Dt
dv_ds = dv_dr # (Only valid in solid)

# This is the displacement field used in the 
# solid formulation; notice that it is NOT 
# in the space V; it is an algebraic object 
# involving the unknown fluid--structure velocity 
# field v.
u = uhat_old + Dt*v

# This will need to be updated to match u, for 
# purposes of setting the boundary condition 
# on the mesh motion subproblem.
u_func = Function(V)

####### Changes of variable #######

# Follow notation from Bazilevs et al., where y is 
# the coordinate in the reference domain, x is the 
# coordinate in the spatial domain, and X is the 
# coordinate in the material domain.  Note that 
# the built-in differential operators (e.g., grad, 
# div, etc.) and integration measures (e.g., dx, ds, 
# etc.) are w.r.t. the reference configuration, y, 
# which is the mesh that FEniCS sees.  
dX = dx(SOLID_FLAG)
dy = dx
ds_y = ds
grad_y = grad
grad_X = grad # (Only valid in solid)
y = SpatialCoordinate(mesh)
x = y + uhat
det_dxdy = det(grad_y(x))
def grad_x(f):
    return dot(grad_y(f),inv(grad_y(x)))
def div_x(f): # (For vector-valued f)
    return tr(grad_x(f))
def div_x_tens(f): # (For (rank-2) tensor-valued f)
    i,j = indices(2)
    return as_tensor(grad_x(f)[i,j,j],(i,))

# Note:  Trying to define dx = det_dxdy*dy would 
# result in an object of type Form, which could no 
# longer be used as an integration measure.
# Thus, integration over the spatial configuration 
# is done with det_dxdy*dy directly in the fluid 
# formulation.  

####### Boundary conditions #######

# BCs for the fluid--structure subproblem:
bc_b_fs = DirichletBC(W.sub(0), Constant((0.0,0.0)), Wall_b())
bc_l_fs = DirichletBC(W.sub(0), Constant((0.0,0.0)), Wall_l())
bc_r_fs = DirichletBC(W.sub(0), Constant((0.0,0.0)), Wall_r())
bc_t_fs = DirichletBC(W.sub(0), Constant((0.0,0.0)), Wall_t())
# Note that "x" in this Expression is really y 
# in the kinematics described above, but, because the 
# mesh motion problem has a zero Dirichlet BC on the inflow
# boundary, there happens to be no difference.
FluidDensity = 3301*4
rho_f = Constant(FluidDensity)
viscosity = 1.2e-1
mu_f = Constant(viscosity)
nu_f = viscosity/FluidDensity
#alpha = 5.0
omega = 17951.94 #((alpha*alpha)*nu_f)/(SOLID_BOTTOM*SOLID_BOTTOM)
alpha = (SOLID_BOTTOM)/sqrt(nu_f/omega)
p_0 = Constant(7584.3)

p_in = Expression("p_0*cos(omega*t)", p_0=p_0, omega=omega,t=0.0,degree=1)
               
bc1_fs = DirichletBC(W.sub(1), p_in, Inflow())
bc2_fs = DirichletBC(W.sub(1), Constant(0), 
                     SolidDomainInterior())                     

#bc3_fs = DirichletBC(W.sub(0).sub(1), Constant(0.0), Inflow())
bcs_fs = [bc_b_fs, bc_l_fs,bc_r_fs,bc_t_fs, bc1_fs, bc2_fs]

# BCs for the mesh motion subproblem:
bc_m_wall_b = DirichletBC(V.sub(1), Constant(0), Wall_b())
bc_m_wall_l = DirichletBC(V, Constant(d*(0,)), Wall_l())
bc_m_wall_r = DirichletBC(V, Constant(d*(0,)), Wall_r())
bc_m_wall_t = DirichletBC(V, Constant(d*(0,)), Wall_t())
bc_m_inflow = DirichletBC(V.sub(0), Constant(0), Inflow())
bc_m_outflow = DirichletBC(V.sub(0), Constant(0), Outflow())
bc_m_struct = DirichletBC(V, u_func, SolidDomainClosure())
bcs_m = [bc_m_struct,bc_m_wall_b,bc_m_wall_l,bc_m_wall_r,bc_m_wall_t,bc_m_inflow,bc_m_outflow]


####### Formulation of mesh motion subproblem #######

# Residual for mesh, which satisfies a fictitious elastic problem:
F_m = grad_y(uhat) + I
E_m = 0.5*(F_m.T*F_m - I)
m_jac_stiff_pow = Constant(3)
# Artificially stiffen the mesh where it is getting crushed:
K_m = Constant(1)/pow(det(F_m),m_jac_stiff_pow)
mu_m = Constant(1)/pow(det(F_m),m_jac_stiff_pow)
S_m = K_m*tr(E_m)*I + 2.0*mu_m*(E_m - tr(E_m)*I/3.0)
res_m = (inner(F_m*S_m,grad_y(du)))*dy
Dres_m = derivative(res_m, uhat)


####### Formulation of the solid subproblem #######

# Elastic properties
YoungsModulus = 5.9e+5
PoissonR = 0.45
MaterialBulk = YoungsModulus/(3*(1-2*PoissonR))
ShearModuli = YoungsModulus/(2*(1+PoissonR))
mu_s = Constant(ShearModuli)
K = Constant(MaterialBulk)
rho_s0 = Constant(1e+3)

# Kinematics:
F = grad_X(u) + I  
E = 0.5*(F.T*F - I)
S = K*tr(E)*I + 2.0*mu_s*(E - tr(E)*I/3.0)
f_s = Constant(d*(0,))
res_s = rho_s0*inner(dv_ds - f_s,dv)*dX \
        + inner(F*S,grad_X(dv))*dX

####### Formulation of the fluid subproblem #######

# Galerkin terms:
rho_f = Constant(3301*4)
viscosity = 1.2e-1
mu_f = Constant(viscosity)
sigma_f = 2.0*mu_f*sym(grad_x(v)) - p*I
v_adv = v - vhat
DvDt = dv_dr + dot(grad_x(v),v_adv)
resGal_f = (rho_f*dot(DvDt,dv) + inner(sigma_f,grad_x(dv))
            + dp*div_x(v))*det_dxdy*dy(FLUID_FLAG)

# Stabilization:

# Deformed mesh size tensor in the spatial configuration:
dxi_dy = inv(Jacobian(mesh))
dxi_dx = dxi_dy*inv(grad_y(x))
G = (dxi_dx.T)*dxi_dx

# SUPG/PSPG:
resStrong_f = rho_f*DvDt - div_x_tens(sigma_f)
Cinv = Constant(1.0)
tau_M = 1.0/sqrt(((2*rho_f/Dt)**2) 
                 + inner(rho_f*v_adv,G*(rho_f*v_adv))
                 + Cinv*(mu_f**2)*inner(G,G))
resSUPG_f = inner(tau_M*resStrong_f,
                  rho_f*dot(grad_x(dv),v_adv)
                  + grad_x(dp))*det_dxdy*dy(FLUID_FLAG)
# LSIC/grad-div:
tau_C = 1.0/(tr(G)*tau_M)
resLSIC_f = tau_C*div_x(v)*div_x(dv)*det_dxdy*dy(FLUID_FLAG)

# Stable Neumann BC term, assuming advective 
# form of material time derivative term and transforming
# normal and area element by Nanson's formula:
dsx_dsy_n_x = det_dxdy*inv(grad_y(x).T)*n_y
v_adv_minus = Min(dot(v_adv,dsx_dsy_n_x),Constant(0))
resOutflow_f = -dot(rho_f*v_adv_minus*dv,v)*ds_y(OUTFLOW_FLAG)

traction = -p_in*n_y
gamma = Constant(1.0)
Consistency_term = -(inner(det_dxdy*inv(grad_y(x).T)*traction,dv))*ds_y(INFLOW_FLAG) #+ gamma*rho_f*Min(inner(v-vhat,n_y),Constant(0.0))*inner(det_dxdy*inv(grad_y(x).T)*(v),dv))*ds_y(INFLOW_FLAG)
sgn=1
sigma_adjoint = (2.0*mu_f*sym(grad_x(dv)) + dp*I)
AdjointConsistency = -sgn*dot(det_dxdy*inv(grad_y(x).T)*(sigma_adjoint*n_y),v)*ds_y(INFLOW_FLAG)
ThirdTerm = -dot(rho_f*dot(v_adv,dsx_dsy_n_x)*dv,v)*ds_y(INFLOW_FLAG)
C_pen = Constant(6.0)
penalty = C_pen*(mu_f/sqrt(dot(dsx_dsy_n_x,G*dsx_dsy_n_x)))*dot(dsx_dsy_n_x,(v))*dot(dsx_dsy_n_x,dv)*ds_y(INFLOW_FLAG)


# # Full fluid residual
# res_weakBC = weakDirichletBC(v=v,p=p,dv=dv,dp=dp,g=None,rho=rho_f,mu=mu_f,mesh=mesh,uhat=uhat,vhat=vhat,ds=ds_y(INFLOW_FLAG),
                    # sym=True,C_pen=Constant(1e3),
                    # overPenalize=False,G=G)
res_f = resGal_f + resSUPG_f + resLSIC_f +resOutflow_f+ Consistency_term #+ ThirdTerm + penalty

# Residual of fluid--structure coupled problem:
res_fs = res_f + res_s 
Dres_fs = derivative(res_fs, w)


####### Nonlinear solver setup #######

# Nonlinear solver parameters; set relative tolerances 
# for subproblems tighter than tolerance for coupled problem, 
# to prevent stagnation.
REL_TOL_FSM = 1e-5
REL_TOL_FS = REL_TOL_M = REL_TOL_FSM*1e-1
MAX_ITERS_FSM = 100
MAX_ITERS_M = 100
MAX_ITERS_FS = 100

# Set up nonlinear problem for mesh motion:
problem_m = NonlinearVariationalProblem(res_m, uhat, 
                                        bcs_m, Dres_m)
solver_m = NonlinearVariationalSolver(problem_m)
solver_m.parameters['newton_solver']\
                   ['maximum_iterations'] = MAX_ITERS_M
solver_m.parameters['newton_solver']\
                   ['relative_tolerance'] = REL_TOL_M

# Create variational problem and solver for 
# the fluid--structure problem:
problem_fs = NonlinearVariationalProblem(res_fs, w, 
                                         bcs_fs, Dres_fs)
solver_fs = NonlinearVariationalSolver(problem_fs)
solver_fs.parameters['newton_solver']\
                    ['maximum_iterations'] = MAX_ITERS_FS
solver_fs.parameters['newton_solver']\
                    ['relative_tolerance'] = REL_TOL_FS

####### Time stepping loop #######

# Create files for storing solution:
DataFolder = "FSIGm0pt5Wo2"
vfile = File(DataFolder+"/velocity.pvd")
pfile = File(DataFolder+"/pressure.pvd")
mfile = File(DataFolder+"/mesh.pvd")


# Initialize time and step counter.
t = float(Dt)
count = 0
# Prevent divide-by-zero in relative residual on first
# iteration of first step.
time_array = []
for bc in bcs_fs:
    bc.apply(w.vector())
while t <= TIME_INTERVAL:

    print(80*"=")
    print("  Time step "+str(count+1)+" , t = "+str(t))
    print(80*"=")
    
    # Use the current time in the inflow BC definition.
    p_in.t = t

    # "Quasi-direct" coupling: the fluid and structure 
    # are solved in one system, but the mesh is solved 
    # in a separate block.
    for i in range(0,MAX_ITERS_FSM):

        # Check fluid--structure residual on the moved
        # mesh, and terminate iteration if this residual 
        # is small:
        res_fs_vec = assemble(res_fs)
        for bc in bcs_fs:
            bc.apply(res_fs_vec,w.vector())
        res_norm = norm(res_fs_vec)
        if(i==0):
            res_norm0 = res_norm
        res_rel = res_norm/res_norm0
        print(80*"*")
        print("  Coupling iteration: "+str(i+1)
              +" , Relative residual = "+str(res_rel))
        print(80*"*")
        if(res_rel < REL_TOL_FSM):
            break

        # Solve for fluid/structure velocity and 
        # pressure at current time:
        solver_fs.solve()
        
        # Update Function in V to be used in mesh 
        # motion BC.  (There are workarounds to avoid 
        # this projection (which requires a linear
        # solve), but projection is most concise for 
        # illustration.)
        u_func.assign(project(u,V))

        # Mesh motion problem; updates uhat at current 
        # time level:
        solver_m.solve()
        
   
    # Extract solutions:
    (v, p) = w.split()

    Qin = []
    Qin1 = assemble(dot(v,n_y)*ds_y(INFLOW_FLAG))
    Qin2 = assemble(dot(v,dsx_dsy_n_x)*ds_y(INFLOW_FLAG))
    Qout1 = assemble(dot(v,n_y)*ds_y(OUTFLOW_FLAG))
    Qout2 = assemble(dot(v,dsx_dsy_n_x)*ds_y(OUTFLOW_FLAG))    
    Qin.append([Qin1,Qin2,Qout1,Qout2])
    
    ypoints_post = np.linspace(0.0,SOLID_BOTTOM,500)
    xpoints_post = np.linspace(0.0,SOLID_RIGHT,500)

    ypoints_Vmid = np.linspace(0.0,SOLID_BOTTOM,351)
    xpoints_Vmid = [0.0,SOLID_RIGHT/2,SOLID_RIGHT]
    
    xpoints_disp = np.linspace(0.0,SOLID_RIGHT,500)
    ypoints_disp = [SOLID_BOTTOM]    

    p_mean = []
    if (count%1==0 and (count+1)>=2100):
        for j in range(len(xpoints_post)):
           
            
            POI_p = ([[0.0] * 2] * 500);
            p_output = []
            
            for k in range(len(ypoints_post)):
                POI_p[k] = (xpoints_post[j],ypoints_post[k])
                p_output.append(p(POI_p[k]))
                
            p_mean.append(sum(p_output)/len(p_output))
        
        pmeanfname = "pmean" + str("%0.6d" % (count/1)) + ".csv"
        pathname = os.path.join(os.getcwd(),DataFolder,pmeanfname)
        np.savetxt(pathname, p_mean, delimiter=',')
     
        qinfname = "qin" + str("%0.6d" % (count/1)) + ".csv"
        pathname = os.path.join(os.getcwd(),DataFolder,qinfname)
        np.savetxt(pathname, Qin, delimiter=',')  
    	
    if (count%1==0 and (count+1)>=2100):
        for jj in range(len(xpoints_Vmid)):
           
            
            POI_v = ([[0.0] * 2] * 351);
            v_output = []
            
            for kk in range(len(ypoints_Vmid)):
                POI_v[kk] = (xpoints_Vmid[jj],ypoints_Vmid[kk])
                
                v_output.append(v(POI_v[kk]))

            if (jj==0):
                vfname = "vDataIn" + str("%0.6d" % (count/1)) + ".csv"
                pathname = os.path.join(os.getcwd(),DataFolder,vfname)
                np.savetxt(pathname, v_output, delimiter=',')  
            if (jj==1):
                vfname = "vDataMid" + str("%0.6d" % (count/1)) + ".csv"
                pathname = os.path.join(os.getcwd(),DataFolder,vfname)
                np.savetxt(pathname, v_output, delimiter=',')  
            if (jj==2):
                vfname = "vDataOut" + str("%0.6d" % (count/1)) + ".csv"
                pathname = os.path.join(os.getcwd(),DataFolder,vfname)
                np.savetxt(pathname, v_output, delimiter=',')  

    disp_data = []
    if (count%1==0 and (count+1)>=2100):
        for jjjj in range(len(ypoints_disp)):
            POI_disp = ([[0.0] * 2] * 500);         
            for kkkk in range(len(xpoints_disp)):
                POI_disp[kkkk] = (xpoints_disp[kkkk],ypoints_disp[jjjj])
                disp_data.append(uhat(POI_disp[kkkk]))
        
        dispdatafname = "dispdata" + str("%0.6d" % (count/1)) + ".csv"
        pathname = os.path.join(os.getcwd(),DataFolder,dispdatafname)
        np.savetxt(pathname, disp_data, delimiter=',')

    #Save to file
    if ((count+1)==2100 or (count+1)==2275 or (count+1)==2450):
        v.rename("v","v")
        p.rename("p","p")
        uhat.rename("u","u")
        vfile << v
        pfile << p
        mfile << uhat

    
    if ((count+1)==2450):
    	break   
    
     
    # Move to next time step:
    uhat_old.assign(uhat)
    w_old.assign(w)
    count += 1
    t += float(Dt)

