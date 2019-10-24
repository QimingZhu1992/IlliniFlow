#include <fstream>
#include <iostream>
#include <dolfin.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include "ns.h"
#include "calcd.h"
#include "calQ.h"

#define tpi    3.141592653589793
#define distol 1e-5
#define r_sim  2.0
#define r_sp   0.1
#define zin    2.0
#define zout   -2.0

#define out_axis_x 0.0
#define out_axis_y 0.0
#define out_axis_z 1.0
#define cen_point_x 0.0
#define cen_point_y 0.0
#define rot_speed 5.0

using namespace dolfin;

// Do not use abs function use std::abs!!!

// Define surface domain
class SurfaceDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		double len=0.0;
		len=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
		len=std::sqrt(len);
		len=std::abs(len - r_sp);

		return bool ( on_boundary and len<1e-3);
	}
};

// Define noslip domain
class NoslipDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		double len=0.0;
		len=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
		len=std::sqrt(len);
		len=std::abs(len - r_sp);

		return bool ( on_boundary and len<1e-3);
	}
};

// Define pressure domain
class PreDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		return bool ( on_boundary and std::abs(x[2]-zout) < 1e-3 );
	}
};

// Define inlet doamin
class InletDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		return bool ( on_boundary and std::abs(x[2]-zin) < 1e-3 );
	}
};

// Define outflow domain
class OutletDomain : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		return bool ( on_boundary and std::abs(x[2]-zout) < 1e-3 );
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

// body force
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

// Inlet boundary condition here
class Incal : public dolfin::Expression
{
	public:
	
	Incal()
	: Expression(3) {}

	void eval(Array<double>& values, const Array<double>& x) const
	{
		values[0]=0.0;
		values[1]=0.0;
		values[2]=-2.0;
	}

};

// Normalize vector with dimention num
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

// Cross product y = u x v
void cross_product(double* u,double* v,double* y)
{
	y[0]=u[1]*v[2]-u[2]*v[1];
	y[1]=u[2]*v[0]-u[0]*v[2];
	y[2]=u[0]*v[1]-u[1]*v[0];
};

// Rotation of mesh
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

		rad_axis[0]=x[0]-cen_point_x;
		rad_axis[1]=x[1]-cen_point_y;
		rad_axis[2]=0.0;
		out_axis[0]=out_axis_x;
		out_axis[1]=out_axis_y;
		out_axis[2]=out_axis_z;
		omega=rot_speed;
		xita =omega*_tnow;
		rd   =std::sqrt(rad_axis[0]*rad_axis[0]+rad_axis[1]*rad_axis[1]);

		normalize(rad_axis,3);
		normalize(out_axis,3);

		//Assume rotation axis is z axis, counter-clock-wise
		//cos(xita+omega) sin(xita+omega)
		tmp[0]=rad_axis[0]*cos(xita)-rad_axis[1]*sin(xita);
		tmp[1]=rad_axis[1]*cos(xita)+rad_axis[0]*sin(xita);
		tmp[2]=0.0;

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
		values[0]= 0.0;
		values[1]= 0.0;
		values[2]=-2.0;
		values[3]= 0.0;
	}

};

int main(int argc, char* argv[])
{
	init(argc, argv);

	// Backend
	//parameters["linear_algebra_backend"] = "PETSc";

	//auto mesh = std::make_shared<Mesh>("sphere.xml");
	auto mesh = std::make_shared<Mesh>(MPI_COMM_WORLD);
	auto f1 = HDF5File(MPI_COMM_WORLD,"sp-coarse.h5", "r");
	f1.read(*mesh, "mesh", false);

	auto markers = std::make_shared<MeshFunction<std::size_t>>(mesh, mesh->topology().dim()-1, 1);
	*markers = 0;
	auto surface_domain  = std::make_shared<SurfaceDomain>();
	surface_domain->mark(*markers, 1);

	/*
	auto sub_domains = std::make_shared<MeshFunction<std::size_t>>(mesh, mesh->topology().dim() - 1);
	*sub_domains = 0;
	auto noslip_domain  = std::make_shared<NoslipDomain>();
	noslip_domain->mark(*sub_domains, 1);
	*/

	//auto periodic_boundary = std::make_shared<PeriodicBoundary>();
	//auto W = std::make_shared<ns::FunctionSpace>(mesh,periodic_boundary);
	auto W = std::make_shared<ns::FunctionSpace>(mesh);
	auto F = std::make_shared<ns::LinearForm>(W);
	auto J = std::make_shared<ns::JacobianForm>(W, W);

	auto Wg= std::make_shared<calQ::FunctionSpace>(mesh);
	auto Fg= std::make_shared<calQ::LinearForm>(Wg);
	auto Jg= std::make_shared<calQ::BilinearForm>(Wg,Wg);

	//J->ds = sub_domains;
	//F->ds = sub_domains;

	auto zero 	   = std::make_shared<Constant>(0.0);
	auto zero_vector   = std::make_shared<Constant>(0.0, 0.0, 0.0);
	auto inflow_vector = std::make_shared<Incal>();
	auto ur            = std::make_shared<urotation>(0.0);

	auto noslip_domain = std::make_shared<NoslipDomain>();
	auto inletdomain   = std::make_shared<InletDomain>();
	auto outletdomain  = std::make_shared<OutletDomain>();

	// Strong boundary condition is used here
	DirichletBC noslip (W->sub(0)         , zero_vector  , noslip_domain);
	DirichletBC inlet  (W->sub(0)         , inflow_vector, inletdomain  );
	DirichletBC outlet (W->sub(1)         , zero         , outletdomain );
	std::vector<DirichletBC*> bcs = {{&noslip, &inlet, &outlet}};
	// Use for incoming turbulence here:
	//std::vector<DirichletBC*> bco = {{&noslip, &outlet}};

	auto Wn   = std::make_shared<Function>(W);
	auto Wo   = std::make_shared<Function>(W);
	auto Wo1  = std::make_shared<Function>(W);
	auto Wmean= std::make_shared<Function>(W);
	auto Wtmp = std::make_shared<Function>(W);
	auto Wrms = std::make_shared<Function>(W);

	auto Qcal = std::make_shared<Function>(Wg);

	// Read initial condition from user define function:
	InitialConditions win;
	*Wo1 = win;

	// Read initial condition from Stokes problem:
	/*
	std::string wina="W1";
        auto f2 = HDF5File(MPI_COMM_WORLD,"W1.h5", "r");
        f2.read(*Wo1, wina);
	*/

	// Read initial condition from previous time step:
	/*
	int datanumk=11000;
	std::string wina="sdata/Wa"+std::to_string(datanumk)+".h5";
        std::string winb="sdata/Wb"+std::to_string(datanumk)+".h5";
	std::string wpa ="Wa"      +std::to_string(datanumk)+".h5";
	std::string wpb ="Wb"      +std::to_string(datanumk)+".h5";
        auto f2 = HDF5File(MPI_COMM_WORLD,wpa, "r");
        f2.read(*Wo, wina);
        auto f3 = HDF5File(MPI_COMM_WORLD,wpb, "r");
        f3.read(*Wo1, winb);
	*/

	auto du0= std::make_shared<Function>((*Wo)[0]);
	auto p0 = std::make_shared<Function>((*Wo)[1]);
	auto du1= std::make_shared<Function>((*Wn)[0]);
	auto p1 = std::make_shared<Function>((*Wn)[1]);
	auto u0 = std::make_shared<Function>((*Wo1)[0]);

	double t   = 0.0;
	double tinc= 0.002;
	double tinv= 1.0/tinc;
	int    tnum= 0;
	int    tref= tnum + 1;
	int    ttol= tnum + 40000;
	double T   = tinc*ttol;

	auto k  = std::make_shared<Constant>(tinc);
	auto idt= std::make_shared<Constant>(tinv);
	auto vu = std::make_shared<Constant>(1e-3);
	double fxval = 0.0;
	auto fx = std::make_shared<Source>(fxval);

	// Collect coefficient into groups
	// Assign the parameters in Residual and Jacobian
	//std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients
	//= {{"fx", fx}, {"Wn", Wn}, {"du0", du0}, {"u0", u0}, {"k", k}, {"idt", idt}, {"vu", vu}};
	std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients
	= {{"fx", fx}, {"Wn", Wn}, {"du0", du0}, {"u0", u0}, {"k", k}, {"idt", idt}, {"vu", vu}, {"ur", ur}};

	// Add extra coefficient for residual
	std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients_F = coefficients;

	J->set_coefficients(coefficients);
	F->set_coefficients(coefficients_F);

	auto A = std::make_shared<PETScMatrix>();

	// Define dof mapping for velocity and pressure
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

	// Use fieldsplit for preconditioning
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCFIELDSPLIT);

	PCFieldSplitSetIS(pc, "u", is[0]);
	PCFieldSplitSetIS(pc, "p", is[1]);
	dolfin::PETScOptions::set("ksp_view");
	dolfin::PETScOptions::set("ksp_monitor_true_residual");
	dolfin::PETScOptions::set("ksp_pc_side", "right");
	// Schur decomposition define
	dolfin::PETScOptions::set("pc_fieldsplit_type", "schur");
	dolfin::PETScOptions::set("pc_fieldsplit_schur_fact_type", "upper");
	dolfin::PETScOptions::set("pc_fieldsplit_schur_preconditioning", "selfp");
	// Jacobi is chosen for velocity preconditioning
	dolfin::PETScOptions::set("fieldsplit_u_ksp_type", "preonly");
	dolfin::PETScOptions::set("fieldsplit_u_pc_type", "jacobi");
	dolfin::PETScOptions::set("fieldsplit_p_ksp_type", "preonly");
	// AMG is chosen for pressure preconditioning
	dolfin::PETScOptions::set("fieldsplit_p_pc_type", "hypre");
	dolfin::PETScOptions::set("fieldsplit_p_pc_hypre_type", "boomeramg");
	// Long range interpolation is crucial for efficiency of AMG
	dolfin::PETScOptions::set("fieldsplit_p_pc_hypre_boomeramg_coarsen_type", "pmis");
        dolfin::PETScOptions::set("fieldsplit_p_pc_hypre_boomeramg_interp_type", "FF1");
	dolfin::PETScOptions::set("fieldsplit_p_pc_hypre_boomeramg_strong_threshold", "0.5");
	KSPSetFromOptions(ksp);

	// Nonlinear problem define
	NSequation nseq(F, J, bcs);
	// Newton solver define
	NewtonSolver newton_solver(MPI_COMM_WORLD, solver, PETScFactory::instance());
	newton_solver.parameters["relative_tolerance"] = 1.0e-5;
	newton_solver.parameters["absolute_tolerance"] = 1.0e-5;
	newton_solver.parameters["maximum_iterations"] = 10;

	std::string velname   ="results/vel";
	std::string prename   ="results/pre";
	std::string meanname  ="results/mean";
	std::string rmsname   ="results/rms";
	std::string Qfname    ="results/Qcal";
	std::string endname="-end";
	std::string e1name =".pvd";
	std::string e2name =".h5";
	std::string Qename =".pvd";
	std::string waname ="sdata/Wa";
        std::string wbname ="sdata/Wb";
        std::string sname,mname,rname,was,wbs,Qname;

	// Parameters for general alpha method
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

	// Incoming turbulence vector here:
	/*
	std::vector<double> corv;
	std::vector<double> resv;
	std::vector<DirichletBC*> bcinflow;
	std::vector<DirichletBC>  bclist;
	std::vector<DirichletBC*> bcs;
	for(int i=0;i<bco.size();i++){
                bcs.push_back(bco[i]);
        }
	*/

	// Read incoming turbulence for first step:
	/*
	int numvin;
	std::string pos_name="veldat/pos.txt";
	std::string read_name="veldat/res"+std::to_string(tnum+1)+".txt";
	read_vel(corv,pos_name,numvin);
	read_vel(resv,read_name,numvin);
        get_inletbc(W, bcs, bclist, corv, resv, numvin, 1);
	NSequation nseq(F, J, bcs);
	*/

	while (tnum < ttol)
	{
		// Update for next time step
		t += tinc;
		tnum = tnum + 1;

		// Mesh rotation
		{
			MeshGeometry& geo=mesh->geometry();
			std::array<double,3> xr;
			for(int i=0;i<geo.num_vertices();i++){
				double tmp_r  =  std::sqrt(geo.x(i,0)*geo.x(i,0) + geo.x(i,1)*geo.x(i,1));
				double tmp_cos;
				double tmp_sin;
				double tmp[3];
				double xita=rot_speed*af*tinc;
				if(tmp_r>1e-4)
				{
					tmp_cos = geo.x(i,0)/tmp_r;
					tmp_sin = geo.x(i,1)/tmp_r;
				}
				else
				{
					tmp_cos = 0.0;
					tmp_sin = 0.0;
				}
				tmp[0]=(tmp_cos*cos(xita) - tmp_sin*sin(xita))*tmp_r;
				tmp[1]=(tmp_sin*cos(xita) + tmp_cos*sin(xita))*tmp_r;
				tmp[2]=geo.x(i,2);

				xr[0] =tmp[0];
				xr[1] =tmp[1];
				xr[2] =tmp[2];

				geo.set(i,xr.data());					
			}
		}

		// Calcute mesh velocity here:
		ur->_tnow=0.0;

		// Get the rank of processor:
		int rankvalues=dolfin::MPI::rank(MPI_COMM_WORLD);

		// Read the incoming turbulence data from txt each step here:
		/*
		std::string read_st="veldat/res";
		int rankvalues=dolfin::MPI::rank(MPI_COMM_WORLD);
		read_name=read_st+std::to_string(tnum)+".txt";
		read_vel(resv,read_name,numvin);
		get_inletbc(W, bcs, bclist, corv, resv, numvin, 2);
		nseq.bcschange(bcs);
		*/

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

		// Mesh rotate:
		{
			MeshGeometry& geo=mesh->geometry();
			std::array<double,3> xr;
			for(int i=0;i<geo.num_vertices();i++){
				double tmp_r  =  std::sqrt(geo.x(i,0)*geo.x(i,0) + geo.x(i,1)*geo.x(i,1));
				double tmp_cos;
				double tmp_sin;
				double tmp[3];
				double xita=rot_speed*(1.0-af)*tinc;
				if(tmp_r>1e-4)
				{
					tmp_cos = geo.x(i,0)/tmp_r;
					tmp_sin = geo.x(i,1)/tmp_r;
				}
				else
				{
					tmp_cos = 0.0;
					tmp_sin = 0.0;
				}
				tmp[0]=(tmp_cos*cos(xita) - tmp_sin*sin(xita))*tmp_r;
				tmp[1]=(tmp_sin*cos(xita) + tmp_cos*sin(xita))*tmp_r;
				tmp[2]=geo.x(i,2);

				xr[0] =tmp[0];
				xr[1] =tmp[1];
				xr[2] =tmp[2];

				geo.set(i,xr.data());					
			}
		}

		// Save velocity
		if (tnum%500==0)
		{
			sname = velname+std::to_string(tnum)+e1name;
			File ufile(sname);
			auto u1 = (*Wo1)[0];
			ufile << u1;

			//mname = prename+std::to_string(tnum)+e1name;
			//File pfile(mname);
			//auto p1 = (*Wn)[1];
			//pfile << p1;
		}

		// Calculate the force coefficient here:
		if(tnum%1==0)
		{
			auto ucof = std::make_shared<Function>((*Wo1)[0]);
			auto pcof = std::make_shared<Function>((*Wo)[1]);

			auto calfx = std::make_shared<calcd::Form_fx>(mesh, ucof, pcof);
			auto calfy = std::make_shared<calcd::Form_fy>(mesh, ucof, pcof);
			auto calfz = std::make_shared<calcd::Form_fz>(mesh, ucof, pcof);

			calfx->ds = markers;
			calfy->ds = markers;
			calfz->ds = markers;

			double calf1 = assemble(*calfx);
			double calf2 = assemble(*calfy);
			double calf3 = assemble(*calfz);

			fxcof.push_back(calf1);
			fycof.push_back(calf2);
			fzcof.push_back(calf3);
			if(tnum%10==0)
			{
				std::string cof_name_st="cc/cof";
				std::string cof_name_en=".txt";
				std::string cof_name=cof_name_st+std::to_string(tnum)+"_"+std::to_string(rankvalues)+cof_name_en;
				std::ofstream calcofFile(cof_name);
				for(int i=0;i<fxcof.size();i++){
					calcofFile << fxcof[i] << " " << fycof[i] << " " << fzcof[i];
					calcofFile << "\n";
				}
				calcofFile.close();
			}
		}

		// Calculate Q by projection
		if (tnum%500==0)
		{
			was   = waname+std::to_string(tnum)+e2name;
                        wbs   = wbname+std::to_string(tnum)+e2name;
                        auto fout1 = HDF5File(MPI_COMM_WORLD,was, "w");
                        fout1.write(*Wo,was);
                        auto fout2 = HDF5File(MPI_COMM_WORLD,wbs, "w");
                        fout2.write(*Wo1,wbs);

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

		

	}

	return 0;
}
