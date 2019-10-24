#include <fstream>
#include <iostream>
#include <dolfin.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <unordered_map>
#include "ns.h"
#include "calcd.h"
#include "calQ.h"

#define meshnumk 5031
#define tpi    3.141592653589793
#define distol 1e-5
#define r_sim  1.49385
#define r_sp   0.55
#define xin    0.542234
#define xout   3.542234
#define yin    -1.49385
#define zin    -1.49385

#define out_axis_x -1.0
#define out_axis_y  0.0
#define out_axis_z  0.0
#define cen_point_x 0.0
#define cen_point_y 0.0
#define cen_point_z 0.0
#define blade_x_min 0.8
#define blade_x_max 1.8
#define rot_speed 28.125
#define water_speed 1.5

using namespace dolfin;

// Do not use abs function use std::abs!!!

// Define surface domain
class SurfaceDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		double len=0.0;
		len=x[1]*x[1]+x[2]*x[2];
		len=std::sqrt(len);

		return bool ( on_boundary and len<r_sp and x[0]>blade_x_min and x[0]<blade_x_max);
	}
};

// Define noslip domain
class NoslipDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		double len=0.0;
		len=x[1]*x[1]+x[2]*x[2];
		len=std::sqrt(len);

		return bool ( on_boundary and len<r_sp and x[0]>blade_x_min and x[0]<blade_x_max);
	}
};

// Define side domain
class SymmetricDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		double len=0.0;
		len=x[1]*x[1]+x[2]*x[2];
		len=std::sqrt(len);

		return bool ( on_boundary and len>r_sp and x[0] > (xin + 1e-3) and x[0] < (xout - 1e-3) );
	}
};

// Define pressure domain
class PreDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		return bool ( on_boundary and std::abs(x[0]-xout) < 1e-3 );
	}
};

// Define inlet doamin
class InletDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		return bool ( on_boundary and std::abs(x[0]-xin) < 1e-3 );
	}
};

// Define outflow domain
class OutletDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		return bool ( on_boundary and std::abs(x[0]-xout) < 1e-3 );
	}
};

// User defined nonlinear problem
class NSequation : public NonlinearProblem
{
	public:

	// Constructor
	NSequation(std::shared_ptr<const Form> F,
			std::shared_ptr<const Form> J,
			std::vector<DirichletBC*> bcs) : _F(F), _J(J), _bcs(bcs) {}

	// User defined residual vector
	void F(GenericVector& b, const GenericVector& x)
	{
		assemble(b, *_F);
		for (std::size_t i = 0; i < _bcs.size(); i++)
			_bcs[i]->apply(b, x);
	}

	// User defined assemble of Jacobian
	void J(GenericMatrix& A, const GenericVector& x)
	{
		assemble(A, *_J);
		for (std::size_t i = 0; i < _bcs.size(); i++)
			_bcs[i]->apply(A);
	}

	void bcschange(std::vector<DirichletBC*> bcnew)
	{
		for (std::size_t i = 0; i < _bcs.size(); i++)
                        _bcs[i]=bcnew[i];
	}

	private:

	// Forms
	std::shared_ptr<const Form> _F;
	std::shared_ptr<const Form> _J;
	std::vector<DirichletBC*> _bcs;
};

class Source : public dolfin::Expression
{
	public:
	
	Source(double fx)
	: Expression(3), _fx(fx) {}

	void eval(Array<double>& values, const Array<double>& x) const
	{
		values[0]= _fx;
		values[1]= 0.0;
		values[2]= 0.0;
	}

	private:
		double _fx;
};

class Incal : public dolfin::Expression
{
	public:
  		Incal(std::vector<double>cor, std::vector<double>res,int nums):Expression(3){
			corin.resize(nums);
			resin.resize(nums);
			for(int i=0;i<nums;i++){
				corin[i]=cor[i];
				resin[i]=res[i];
			}
			numk=nums;
		}
  		void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  		{
			double tol=1e-5;
			for(int i=0;i<numk;i=i+3){
				double err=std::abs(x[0]-corin[i])+std::abs(x[1]-corin[i+1])+std::abs(x[2]-corin[i+2]);
				if(err<tol){
					values[0]=resin[i];
					values[1]=resin[i+1];
					values[2]=resin[i+2];
					break;
				}
			}
		}
	private:
		std::vector<double> corin;
		std::vector<double> resin;
		int numk;
};

void normalize(double* x,int num)
{
	double norm;
	norm=0.0;
	for(int i=0;i<num;i++){
		norm+=x[i]*x[i];
	}
	norm=sqrt(1e-10+norm);
	for(int i=0;i<num;i++){
		x[i]=x[i]/norm;
	}
};

void cross_product(double* u,double* v,double* y)
{
	y[0]=u[1]*v[2]-u[2]*v[1];
	y[1]=u[2]*v[0]-u[0]*v[2];
	y[2]=u[0]*v[1]-u[1]*v[0];
};

class urotation : public dolfin::Expression
{
	public:
	
	urotation(double tnow)
	: Expression(3), _tnow(tnow) {}

	void eval(Array<double>& values, const Array<double>& x) const
	{
		double rad_axis[3];
		double out_axis[3];
		double rot_axis[3];
		double tmp[3];
		double omega;
		double xita;
		double rd;

		rad_axis[0]=0.0;
		rad_axis[1]=x[1]-cen_point_y;
		rad_axis[2]=x[2]-cen_point_z;
		out_axis[0]=out_axis_x;
		out_axis[1]=out_axis_y;
		out_axis[2]=out_axis_z;
		omega=rot_speed;
		xita =omega*_tnow;
		rd   =std::sqrt(rad_axis[1]*rad_axis[1]+rad_axis[2]*rad_axis[2]);

		normalize(rad_axis,3);
		normalize(out_axis,3);

		//Assume rotation axis is z axis, counter-clock-wise
		//cos(xita+omega) sin(xita+omega)
		tmp[0]=0.0;
		tmp[1]=rad_axis[1]*cos(xita)+rad_axis[2]*sin(xita);
		tmp[2]=rad_axis[2]*cos(xita)-rad_axis[1]*sin(xita);

		for(int i=0;i<3;i++)
		{
			rad_axis[i]=tmp[i];
		}

		cross_product(out_axis,rad_axis,rot_axis);

		values[0]= omega*rd*rot_axis[0];
		values[1]= omega*rd*rot_axis[1];
		values[2]= omega*rd*rot_axis[2];
	}

	public:
		double _tnow;
};

// Define initial condition
class InitialConditions : public Expression
{
	public:

	InitialConditions() : Expression(4) 
	{
		dolfin::seed(2 + dolfin::MPI::rank(MPI_COMM_WORLD));
	}

	void eval(Array<double>& values, const Array<double>& x) const
	{		
		double r=std::sqrt(x[0]*x[0]+x[1]*x[1]);
		//values[0]= 4.0*(1.0-r/r_sim)*(1.0-r/r_sim);
		values[0]= water_speed;
		values[1]= 0.0;
		values[2]= 0.0;
		values[3]= 0.0;
	}

};

void get_inletbc(std::shared_ptr<const FunctionSpace> W, std::vector<DirichletBC*>& bcnew, std::vector<DirichletBC>& bclist,
			std::vector<double>& cor, std::vector<double>& res, int nums, int flag)
{
	auto Inprofile = std::make_shared<InletDomain>();
	auto velin     = std::make_shared<Incal>(cor,res,nums);

	std::vector<DirichletBC>().swap(bclist);

	DirichletBC bcin(W->sub(0), velin, Inprofile);
	bclist.push_back(bcin);
	if(flag==1)
		bcnew.push_back(&(bclist[0]));
	else
	{
		bcnew.pop_back();
		bcnew.push_back(&(bclist[0]));
	}
};

void read_vel(std::vector<double>& dat, std::string name, int& nums)
{
	int numk=0;
	int num =meshnumk*3;
	int rankvalues=dolfin::MPI::rank(MPI_COMM_WORLD);

	std::vector<double>().swap(dat);

	std::ifstream infile;
	infile.open(name, std::ios::in);
	
	if(!infile.is_open() )
		std::cout << "Open file fail." <<std::endl;
	//infile >> num;

	dat.resize(num);
	while(numk<num)
	{
		//The format of inlet data is px py pz or ux uy uz
		infile >> dat[numk] >> dat[numk+1] >> dat[numk+2];
		numk=numk+3;
	}
	nums=num;
	infile.close();
};

void rotate_cal(std::vector<double> &corv, double tinc, int num)
{
	for(int i=0;i<num;i=i+3){
		double tmp_r  =  std::sqrt(corv[i+1]*corv[i+1] + corv[i+2]*corv[i+2]);
		double tmp_cos;
		double tmp_sin;
		double tmp[3];
		double xita=rot_speed*tinc;
		if(tmp_r>1e-4)
		{
			tmp_cos = corv[i+2]/tmp_r;
			tmp_sin = corv[i+1]/tmp_r;
		}
		else
		{
			tmp_cos = 0.0;
			tmp_sin = 0.0;
		}
		corv[i+1]=(tmp_sin*cos(xita) + tmp_cos*sin(xita))*tmp_r;
		corv[i+2]=(tmp_cos*cos(xita) - tmp_sin*sin(xita))*tmp_r;
	}
};

int main(int argc, char* argv[])
{
	init(argc, argv);

	// Backend
	//parameters["linear_algebra_backend"] = "PETSc";

	//auto mesh = std::make_shared<Mesh>("sphere.xml");
	auto mesh = std::make_shared<Mesh>(MPI_COMM_WORLD);
	auto f1 = HDF5File(MPI_COMM_WORLD,"mesh3000.h5", "r");
	f1.read(*mesh, "mesh", false);

	auto meshu = std::make_shared<Mesh>(MPI_COMM_WORLD);
	auto f4 = HDF5File(MPI_COMM_WORLD,"rotor.h5", "r");
	f4.read(*meshu, "mesh", false);

	auto markers = std::make_shared<MeshFunction<std::size_t>>(mesh, mesh->topology().dim()-1, 1);
	*markers = 0;
	auto surface_domain  = std::make_shared<SurfaceDomain>();
	surface_domain->mark(*markers, 1);

	//auto periodic_boundary = std::make_shared<PeriodicBoundary>();
	//auto W = std::make_shared<ns::FunctionSpace>(mesh,periodic_boundary);
	auto W = std::make_shared<ns::FunctionSpace>(mesh);
	auto F = std::make_shared<ns::LinearForm>(W);
	auto J = std::make_shared<ns::JacobianForm>(W, W);

	auto Wg= std::make_shared<calQ::FunctionSpace>(mesh);
	auto Fg= std::make_shared<calQ::LinearForm>(Wg);
	auto Jg= std::make_shared<calQ::BilinearForm>(Wg,Wg);

	auto Wu= std::make_shared<ns::FunctionSpace>(meshu);

	J->ds = markers;
	F->ds = markers;

	double t   = 0.0;
	double tinc= 0.0001;
	double tinv= 1.0/tinc;
	int    tnum= 3000;
	int    tref= tnum + 1;
	int    ttol= tnum + 20000;
	double T   = tinc*ttol;
	double vucal=(1e-3)/(1050.0);

	auto zero 	   = std::make_shared<Constant>(0.0);
	auto zero_vector   = std::make_shared<Constant>(0.0, 0.0, 0.0);
	//auto inflow_vector = std::make_shared<Incal>();
	auto ur            = std::make_shared<urotation>(0.0);

	// Inlet read every time step, noslip use weak boundary condition here:
	//auto inletdomain   = std::make_shared<InletDomain>();
	//auto noslip_domain = std::make_shared<NoslipDomain>();
	auto outletdomain  = std::make_shared<OutletDomain>();
	auto symdomain     = std::make_shared<SymmetricDomain>();

	//DirichletBC noslip (W->sub(0)         , ur           , noslip_domain);
	//DirichletBC inlet  (W->sub(0)         , inflow_vector, inletdomain  );
	DirichletBC outlet (W->sub(1)         , zero         , outletdomain );
	DirichletBC sym_y  (W->sub(0)->sub(1) , zero         , symdomain    );
	DirichletBC sym_z  (W->sub(0)->sub(2) , zero         , symdomain    );
	//std::vector<DirichletBC*> bcs = {{&inlet, &outlet, &sym_y, &sym_z}};
	//std::vector<DirichletBC*> bcs = {{&noslip, &inlet, &outlet, &sym_y, &sym_z}};
	// Use for incoming turbulence here:
	//std::vector<DirichletBC*> bco = {{&outlet, &sym_y, &sym_z}};
	std::vector<DirichletBC*> bco = {{&outlet}};

	DirichletBC outlet_un (Wu->sub(1)         , zero         , outletdomain );
	DirichletBC sym_y_un  (Wu->sub(0)->sub(1) , zero         , symdomain    );
	DirichletBC sym_z_un  (Wu->sub(0)->sub(2) , zero         , symdomain    );
	std::vector<DirichletBC*> bcun = {{&outlet_un, &sym_y_un, &sym_z_un}};

	auto Wn   = std::make_shared<Function>(W);
	auto Wo   = std::make_shared<Function>(W);
	auto Wo1  = std::make_shared<Function>(W);
	
	auto Wrmean= std::make_shared<Function>(W);
	auto Wrtmp = std::make_shared<Function>(W);
	auto Wrrms = std::make_shared<Function>(W);
	auto Wqa   = std::make_shared<Function>(W);

	auto Qcal = std::make_shared<Function>(Wg);

	auto Wu1   = std::make_shared<Function>(Wu);
	auto Wu2   = std::make_shared<Function>(Wu);
	auto Wumean= std::make_shared<Function>(Wu);
	auto Wutmp = std::make_shared<Function>(Wu);
	auto Wurms = std::make_shared<Function>(Wu);
	auto Wpmean= std::make_shared<Function>(Wu);

	// Read initial condition from user define function:
	//InitialConditions win;
	//*Wo1 = win;

	// Read initial condition from Stokes problem:
	/*
	std::string wina="W1";
        auto f2 = HDF5File(MPI_COMM_WORLD,"W1.h5", "r");
        f2.read(*Wo1, wina);
	*/

	// Read initial condition from previous time step:
	
	int datanumk=3000;
	std::string wina="sdata/Wa"+std::to_string(datanumk)+".h5";
        std::string winb="sdata/Wb"+std::to_string(datanumk)+".h5";
	std::string wpa ="Wa"      +std::to_string(datanumk)+".h5";
	std::string wpb ="Wb"      +std::to_string(datanumk)+".h5";
        auto f2 = HDF5File(MPI_COMM_WORLD,wpa, "r");
        f2.read(*Wo,  "Wo");
        auto f3 = HDF5File(MPI_COMM_WORLD,wpb, "r");
        f3.read(*Wo1, "Wo1");

	auto du0= std::make_shared<Function>((*Wo)[0]);
	auto p0 = std::make_shared<Function>((*Wo)[1]);
	auto du1= std::make_shared<Function>((*Wn)[0]);
	auto p1 = std::make_shared<Function>((*Wn)[1]);
	auto u0 = std::make_shared<Function>((*Wo1)[0]);

	auto k  = std::make_shared<Constant>(tinc);
	auto idt= std::make_shared<Constant>(tinv);
	auto vu = std::make_shared<Constant>(vucal);
	double fxval = 0.0;
	auto fx = std::make_shared<Source>(fxval);

	// Collect coefficient into groups
	//std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients
	//= {{"fx", fx}, {"Wn", Wn}, {"du0", du0}, {"u0", u0}, {"k", k}, {"idt", idt}, {"vu", vu}};
	std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients
	= {{"fx", fx}, {"Wn", Wn}, {"du0", du0}, {"u0", u0}, {"k", k}, {"idt", idt}, {"vu", vu}, {"ur", ur}};

	// Add extra coefficient for residual
	std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients_F = coefficients;

	J->set_coefficients(coefficients);
	F->set_coefficients(coefficients_F);

	auto A = std::make_shared<PETScMatrix>();

	IS is[2];
	auto u_dofs = W->sub(0)->dofmap()->dofs();
	auto p_dofs = W->sub(1)->dofmap()->dofs();
	dolfin::cout << "Number of u and p dofs: " << u_dofs.size() << ", "
		<< p_dofs.size() << dolfin::endl;
	ISCreateGeneral(PETSC_COMM_WORLD, u_dofs.size(), u_dofs.data(),
		PETSC_COPY_VALUES, &is[0]);
	ISCreateGeneral(PETSC_COMM_WORLD, p_dofs.size(), p_dofs.data(),
		PETSC_COPY_VALUES, &is[1]);

	auto solver = std::make_shared<PETScKrylovSolver>("gmres");
	KSP ksp = solver->ksp();
	PC pc;

	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCFIELDSPLIT);

	PCFieldSplitSetIS(pc, "u", is[0]);
	PCFieldSplitSetIS(pc, "p", is[1]);
	dolfin::PETScOptions::set("ksp_view");
	dolfin::PETScOptions::set("ksp_monitor_true_residual");
	dolfin::PETScOptions::set("ksp_pc_side", "right");
	dolfin::PETScOptions::set("pc_fieldsplit_type", "schur");
	dolfin::PETScOptions::set("pc_fieldsplit_schur_fact_type", "upper");
	dolfin::PETScOptions::set("pc_fieldsplit_schur_preconditioning", "selfp");
	dolfin::PETScOptions::set("fieldsplit_u_ksp_type", "preonly");
	dolfin::PETScOptions::set("fieldsplit_u_pc_type", "jacobi");
	dolfin::PETScOptions::set("fieldsplit_p_ksp_type", "preonly");
	dolfin::PETScOptions::set("fieldsplit_p_pc_type", "hypre");
	dolfin::PETScOptions::set("fieldsplit_p_pc_hypre_type", "boomeramg");
	dolfin::PETScOptions::set("fieldsplit_p_pc_hypre_boomeramg_coarsen_type", "pmis");
        dolfin::PETScOptions::set("fieldsplit_p_pc_hypre_boomeramg_interp_type", "FF1");
	dolfin::PETScOptions::set("fieldsplit_p_pc_hypre_boomeramg_strong_threshold", "0.5");
	KSPSetFromOptions(ksp);

	//NSequation nseq(F, J, bcs);
	NewtonSolver newton_solver(MPI_COMM_WORLD, solver, PETScFactory::instance());
	newton_solver.parameters["relative_tolerance"] = 1.0e-5;
	newton_solver.parameters["absolute_tolerance"] = 1.0e-5;
	newton_solver.parameters["maximum_iterations"] = 20;

	std::string velname   ="results/vel";
	std::string prename   ="results/pre";
	std::string rmeanname ="results/rmean";
	std::string rrmsname  ="results/rrms";
	std::string fmeanname ="results/fmean";
	std::string frmsname  ="results/frms";
	std::string fpmeanname="results/fpmean";
	std::string Qfname    ="results/Qcal";
	std::string endname="-end";
	std::string e1name =".pvd";
	std::string e2name =".h5";
	std::string m1name =".h5";
	std::string Qename =".pvd";
	std::string waname ="sdata/Wa";
        std::string wbname ="sdata/Wb";
	std::string maname ="sdata/mesh";
        std::string sname,mname,rmname,rrname,fmname,frname,fpmname,was,wbs,meshname,Qname;

	double rhoc = 0.5;
	double am   = 0.5*(3.0-rhoc)/(1.0+rhoc);
	double af   = 1.0/(1+rhoc);
	double gamma= 0.5 + am -af;
	double c1   = tinc*(1-gamma);
	double c2   = tinc*gamma;

	// Force coefficient vector here:
	std::vector<double> fxcof;
	std::vector<double> fycof;
	std::vector<double> fzcof;
	std::vector<double> txcof;
	std::vector<double> tycof;
	std::vector<double> tzcof;

	// Incoming turbulence vector here:
	
	std::vector<double> corv;
	std::vector<double> resv;
	std::unordered_map<double,int> corv_hash;
	std::vector<DirichletBC*> bcinflow;
	std::vector<DirichletBC>  bclist;
	std::vector<DirichletBC*> bcs;
	for(int i=0;i<bco.size();i++){
                bcs.push_back(bco[i]);
        }

	// Read incoming turbulence for first step:
	int numvin;
	std::string pos_name="veldat/pos.txt";
	std::string read_name="veldat/res"+std::to_string(tnum+1)+".txt";
	read_vel(corv,pos_name ,numvin);
	read_vel(resv,read_name,numvin);
	rotate_cal(corv,tnum*tinc,numvin);
        get_inletbc(W, bcs, bclist, corv, resv, numvin, 1);
	NSequation nseq(F, J, bcs);

	while (tnum < ttol)
	{
		// Update for next time step
		t += tinc;
		tnum = tnum + 1;

		// Move the mesh to alpha state here:
		{
			MeshGeometry& geo=mesh->geometry();
			std::array<double,3> xr;
			for(int i=0;i<geo.num_vertices();i++){
				double tmp_r  =  std::sqrt(geo.x(i,1)*geo.x(i,1) + geo.x(i,2)*geo.x(i,2));
				double tmp_cos;
				double tmp_sin;
				double tmp[3];
				double xita=rot_speed*af*tinc;
				if(tmp_r>1e-4)
				{
					tmp_cos = geo.x(i,2)/tmp_r;
					tmp_sin = geo.x(i,1)/tmp_r;
				}
				else
				{
					tmp_cos = 0.0;
					tmp_sin = 0.0;
				}
				tmp[0]=geo.x(i,0);
				tmp[1]=(tmp_sin*cos(xita) + tmp_cos*sin(xita))*tmp_r;
				tmp[2]=(tmp_cos*cos(xita) - tmp_sin*sin(xita))*tmp_r;

				xr[0] =tmp[0];
				xr[1] =tmp[1];
				xr[2] =tmp[2];

				geo.set(i,xr.data());
		
			}
			
		}
		// Update my inlet position here:
		if(tnum==tref){
			rotate_cal(corv,af*tinc,numvin);
		}
		else{
			rotate_cal(corv,tinc,numvin);
		}

		// Rotate here:
		ur->_tnow = 0.0;

		// Read the incoming turbulence data from txt each step here:
		std::string read_st="veldat/res";
		int rankvalues=dolfin::MPI::rank(MPI_COMM_WORLD);
		read_name=read_st+std::to_string(tnum)+".txt";
		read_vel(resv,read_name,numvin);
		get_inletbc(W, bcs, bclist, corv, resv, numvin, 2);
		nseq.bcschange(bcs);

		// Initial Guess and Apply boundary conditions
		du1= std::make_shared<Function>((*Wo)[0]);
		p1 = std::make_shared<Function>((*Wo)[1]);
		*(du1->vector()) *= (gamma-1.0)/gamma;
		assign(Wn,{du1,p1});
		for (std::size_t i = 0; i < bcs.size(); i++)
			bcs[i]->apply(*Wn->vector());

		// Solve the Fixed point iterations
		newton_solver.solve(nseq, *Wn->vector());

		// Calculate the new velocity and Apply boundary conditions
		Wo1->vector()->axpy(c1, *Wo->vector());
		Wo1->vector()->axpy(c2, *Wn->vector());
		for (std::size_t i = 0; i < bcs.size(); i++)
			bcs[i]->apply(*Wo1->vector());

		// Update the velocity and derivative of velocity
		du1= std::make_shared<Function>((*Wn)[0]);
		assign(du0,du1);
		du1= std::make_shared<Function>((*Wo1)[0]);
		assign(u0,du1);

		// Save velocity:
		*Wo->vector() = *Wn->vector();

		// Move the mesh to t_n+1 state
		{
			MeshGeometry& geo=mesh->geometry();
			std::array<double,3> xr;
			for(int i=0;i<geo.num_vertices();i++){
				double tmp_r  =  std::sqrt(geo.x(i,1)*geo.x(i,1) + geo.x(i,2)*geo.x(i,2));
				double tmp_cos;
				double tmp_sin;
				double tmp[3];
				double xita=rot_speed*(1.0-af)*tinc;
				if(tmp_r>1e-4)
				{
					tmp_cos = geo.x(i,2)/tmp_r;
					tmp_sin = geo.x(i,1)/tmp_r;
				}
				else
				{
					tmp_cos = 0.0;
					tmp_sin = 0.0;
				}
				tmp[0]=geo.x(i,0);
				tmp[1]=(tmp_sin*cos(xita) + tmp_cos*sin(xita))*tmp_r;
				tmp[2]=(tmp_cos*cos(xita) - tmp_sin*sin(xita))*tmp_r;

				xr[0] =tmp[0];
				xr[1] =tmp[1];
				xr[2] =tmp[2];

				geo.set(i,xr.data());
		
			}
		}

		if (tnum%500==0)
		{
			sname = velname+std::to_string(tnum)+e1name;
			File ufile(sname);
			auto u1 = (*Wo1)[0];
			ufile << u1;

			mname = prename+std::to_string(tnum)+e1name;
			File pfile(mname);
			auto p1 = (*Wn)[1];
			pfile << p1;
		}

		// Calculate the force coefficient here:
		if(tnum%1==0)
		{
			auto ucof = std::make_shared<Function>((*Wo1)[0]);
			auto pcof = std::make_shared<Function>((*Wo)[1]);

			auto calfx = std::make_shared<calcd::Form_fx>(mesh, ucof, pcof);
			auto calfy = std::make_shared<calcd::Form_fy>(mesh, ucof, pcof);
			auto calfz = std::make_shared<calcd::Form_fz>(mesh, ucof, pcof);
			auto caltx = std::make_shared<calcd::Form_tx>(mesh, ucof, pcof);
			auto calty = std::make_shared<calcd::Form_ty>(mesh, ucof, pcof);
			auto caltz = std::make_shared<calcd::Form_tz>(mesh, ucof, pcof);

			calfx->ds = markers;
			calfy->ds = markers;
			calfz->ds = markers;
			caltx->ds = markers;
			calty->ds = markers;
			caltz->ds = markers;

			double calf1 = assemble(*calfx);
			double calf2 = assemble(*calfy);
			double calf3 = assemble(*calfz);
			double calt1 = assemble(*caltx);
			double calt2 = assemble(*calty);
			double calt3 = assemble(*caltz);

			fxcof.push_back(calf1);
			fycof.push_back(calf2);
			fzcof.push_back(calf3);
			txcof.push_back(calt1);
			tycof.push_back(calt2);
			tzcof.push_back(calt3);
			if(tnum%10==0 and rankvalues==0)
			{
				std::string cof_name_st="cc/cof";
				std::string cof_name_en=".txt";
				std::string cof_name=cof_name_st+std::to_string(tnum)+cof_name_en;
				std::ofstream calcofFile(cof_name);
				for(int i=0;i<fxcof.size();i++){
					calcofFile << fxcof[i] << " " << fycof[i] << " " << fzcof[i] << " " << txcof[i] << " " << tycof[i] << " " << tzcof[i];
					calcofFile << "\n";
				}
				calcofFile.close();
			}
		}

		if (tnum%500==0)
		{
			was     = waname+std::to_string(tnum)+e2name;
                        wbs     = wbname+std::to_string(tnum)+e2name;
                        auto fout1 = HDF5File(MPI_COMM_WORLD,was, "w");
                        fout1.write(*Wo ,"Wo");
                        auto fout2 = HDF5File(MPI_COMM_WORLD,wbs, "w");
                        fout2.write(*Wo1,"Wo1");
			meshname= maname+std::to_string(tnum)+m1name;
			auto meshout = HDF5File(MPI_COMM_WORLD, meshname,"w");
			meshout.write(*mesh, "mesh");

			auto uin= std::make_shared<Function>((*Wo1)[0]);
			auto Ag = std::make_shared<PETScMatrix>();
			auto Bg = std::make_shared<PETScVector>();
			Fg->u   = uin;
			assemble_system(*Ag, *Bg, *Jg, *Fg, {});
			auto grad_solver = std::make_shared<PETScKrylovSolver>("gmres");
			grad_solver->set_operator(Ag);
			grad_solver->solve(*Qcal->vector(), *Bg);

			Qname = Qfname+std::to_string(tnum)+Qename;
			File qfile(Qname);
			auto q1 = *Qcal;
			qfile << q1;

		}

		if (tnum%10==0)
		{
			// Calculate the mean and rms value of pressure on rotor here:
			*Wrmean->vector() += *Wn->vector();
			*Wrtmp->vector()   = *Wn->vector();
			*Wrtmp->vector()  *= *Wn->vector();
			*Wrrms->vector()  += *Wrtmp->vector();

			/*
			// Do lagrange interpolation here (nonfitted mesh):
			mesh->bounding_box_tree()->build(*mesh);
			LagrangeInterpolator::interpolate(*Wu1,*Wo1);
			LagrangeInterpolator::interpolate(*Wu2,*Wn);

			// Apply boundary conditions here:
			for (std::size_t i = 0; i < bcun.size(); i++){
				bcun[i]->apply(*Wu1->vector());
				bcun[i]->apply(*Wu2->vector());
			}

			// calculate the mean and rms value of velocity here:	
			*Wumean->vector() += *Wu1->vector();
			*Wutmp->vector()   = *Wu1->vector();
			*Wutmp->vector()  *= *Wu1->vector();
			*Wurms->vector()  += *Wutmp->vector();
			*Wpmean->vector() += *Wu2->vector();
			*/

			if (tnum%500==0)
			{
				rmname = rmeanname+std::to_string(tnum)+e2name;
				auto m1file = HDF5File(MPI_COMM_WORLD,rmname, "w");
				m1file.write(*Wrmean ,"Wrmean");

				rrname = rrmsname+std::to_string(tnum)+e2name;
				auto m2file = HDF5File(MPI_COMM_WORLD,rrname, "w");
				m2file.write(*Wrrms ,"Wrrms");

				/*
				fmname = fmeanname+std::to_string(tnum)+e2name;
				auto m3file = HDF5File(MPI_COMM_WORLD,fmname, "w");
				m3file.write(*Wumean ,"Wumean");

				frname = frmsname+std::to_string(tnum)+e2name;
				auto m4file = HDF5File(MPI_COMM_WORLD,frname, "w");
				m4file.write(*Wurms ,"Wurms");

				fpmname = fpmeanname+std::to_string(tnum)+e2name;
				auto m5file = HDF5File(MPI_COMM_WORLD,fpmname, "w");
				m5file.write(*Wpmean ,"Wpmean");
				*/
			}
		}

	}

	return 0;
}


