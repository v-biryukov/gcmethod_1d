#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

class mesh_1d
{
	double left, right;
	double h;
	int nx;

    int N;

	void read_from_file(std::string path);
public:
    mesh_1d(std::string path) {read_from_file(path); nx = static_cast<int>((right-left)/h + 0.5);}
	void create_mesh(){};
	int get_number_of_points() {return nx+1;};
	double get_point(int n) {return left + h * n;};
    double get_additional_point(int segment_num, int j){return left + h*segment_num + h/N*(j+1);}
	int get_number_of_segments() {return nx;};
	double get_left() {return left;};
	double get_right() {return right;};

    void set_mesh(double left, double right, double h)
    {
        this->left = left;
        this->right = right;
        this->h = h;
        nx = static_cast<int>((right-left)/h + 0.5);
    }
};


struct value
{
    double v;
    double s;
    value(){};
    value (double u1, double u2) : v(u1), s(u2) {}
};

struct rvalue
{
    double w1;
    double w2;
    rvalue(){};
    rvalue (double w1, double w2) : w1(w1), w2(w2) {}
    rvalue operator+(rvalue w)
    {
        return rvalue(this->w1+w.w1, this->w2+w.w2);
    }
    rvalue operator-(rvalue w)
    {
        return rvalue(this->w1-w.w1, this->w2-w.w2);
    }
    rvalue operator*(double a)
    {
        return rvalue(a*this->w1, a*this->w2);
    }
};

class gcmethod_1d
{
	int N;
	double tau, h;
	int number_of_steps;
	bool is_monotonic;
	mesh_1d & mesh;

    int saving_frequency;


	// values:
    std::vector<value> main0;
    std::vector<value> main1;
    std::vector<std::vector<value> > additional0;
    std::vector<std::vector<value> > additional1;

    // params
    std::vector<double> c;
    std::vector<double> rho;


	
	void init();
	void read_from_file(std::string path);
	void save_to_vtk();
	void step();

    rvalue approximate(double x, int n_of_segment, std::vector<rvalue> & ws);
    rvalue approximate_linear(double x, int n_of_segment, std::vector<rvalue> & ws);
    rvalue approximate_quadratic(double x, int n_of_segment, std::vector<rvalue> & ws);
    rvalue approximate_cubic(double x, int n_of_segment, std::vector<rvalue> & ws);


    value approximate_linear(int l_segment, int r_segment);


    rvalue get_omega(const value & u, int segn);
    rvalue get_omega(const value & u, int segnR, int segnL);
    value get_u(const rvalue & w, int segn);
    value get_u(const rvalue & w, int segnR, int segnL);

	void save_to_vtk(std::string name);
    void set_initial();
    value initial_conditions(double x) {
     return value(5.0*exp(-4.0*x*x), 5.0*exp(-4.0*x*x)* c[mesh.get_number_of_segments()/2-10]*rho[mesh.get_number_of_segments()/2-10]);}
     //if (x < -1.0 && x > -2.0) return value(3.0, 3.0 * c[mesh.get_number_of_segments()/2-10]*rho[mesh.get_number_of_segments()/2-10]); else return value(0.0, 0.0);}
    void set_params();
    void set_one_border_params(double cl, double cr, double rl, double rr, double x0);
    double find_divergence_q(std::string type, double cl, double cr, double rl, double rr);
public:
    gcmethod_1d(mesh_1d & mesht, std::string path) : mesh(mesht) {read_from_file(path); init();}
	void calculate();
    void calculate_diverge_contour();

};

void gcmethod_1d::init()
{
    c.resize(mesh.get_number_of_segments());
    rho.resize(mesh.get_number_of_segments());
    set_params();

    main0.resize(mesh.get_number_of_points());
    main1.resize(mesh.get_number_of_points());
    additional0.resize(mesh.get_number_of_segments(), std::vector<value>(N-1));
    additional1.resize(mesh.get_number_of_segments(), std::vector<value>(N-1));

    set_initial();
}

void gcmethod_1d::set_initial()
{
    for (int i = 0; i < mesh.get_number_of_points(); i++ )
    {
        main0.at(i) = initial_conditions(mesh.get_point(i));
    }
    for (int i = 0; i < mesh.get_number_of_segments(); i++ )
        for(int j = 0; j < N-1; j++)
        {
            additional0.at(i).at(j) = initial_conditions(mesh.get_additional_point(i, j));
        }
}

void gcmethod_1d::set_one_border_params(double cl, double cr, double rl, double rr, double x0)
{
    for (int i = 0; i < mesh.get_number_of_segments(); i++ )
    {
        double x = mesh.get_left() + i*h + h/2.0;
        if (x < x0)
        {
            c[i] = cl;
            rho[i] = rl;
        }
        else
        {
            c[i] = cr;
            rho[i] = rr;
            //6.75;(c=1)
            // c=0.8: 9.0-9.1,
        }
    }
    c[0] = cr;
    rho[0] = rr;
}


void gcmethod_1d::set_params()
{
    set_one_border_params(1,1,2,1, (mesh.get_left() + mesh.get_right())/2);
}


void gcmethod_1d::step()
{
    for ( int i = 0; i < mesh.get_number_of_points(); i++ )
	{

        // Find segment numbers
        int snL = i-1;
        if (i==0) snL = mesh.get_number_of_segments()-1;
        int snR = i;
        if (i==mesh.get_number_of_segments()) snR = 0;
        // Find x's
        double xL = mesh.get_point(i) - c[snL] * tau;
        if ( i==0 )
            xL += (mesh.get_right() - mesh.get_left());
        double xR = mesh.get_point(i) + c[snR] * tau;
        if ( i==mesh.get_number_of_segments())
            xR -= (mesh.get_right() - mesh.get_left());

        // Find riemann invs

        rvalue rdata_here = get_omega(main0[i], snR, snL);

        // R or -

        std::vector<rvalue> wsR;
        wsR.push_back(get_omega(main0[snR], snR));
        for (int k = 0; k < N-1; k++)
            wsR.push_back(get_omega(additional0[snR][k], snR));
        wsR.push_back(get_omega(main0[snR+1], snR));
        double diff_wR = approximate(xR, snR, wsR).w1 - rdata_here.w1;

        // L or +
        std::vector<rvalue> wsL;
        wsL.push_back(get_omega(main0[snL], snL));
        for (int k = 0; k < N-1; k++)
            wsL.push_back(get_omega(additional0[snL][k], snL));
        wsL.push_back(get_omega(main0[snL+1], snL));
        double diff_wL = approximate(xL, snL, wsL).w2 - rdata_here.w2;

        rvalue diff_w = rvalue(diff_wR, diff_wL);

        value temp = get_u(diff_w, snR, snL);

        main1.at(i).v = main0.at(i).v + temp.v;
        main1.at(i).s = main0.at(i).s + temp.s;
	}

	for ( int i = 0; i < mesh.get_number_of_segments(); i++ )
		for (int j = 0; j < N-1; j++)
		{
            int snL = i;
            int snR = i;
            // Find x's
            double xL = mesh.get_additional_point(i, j) - c[snL] * tau;

            double xR = mesh.get_additional_point(i, j) + c[snR] * tau;

            // Find riemann invs

            rvalue rdata_here = get_omega(additional0[i][j], snR, snL);

            // R or -

            std::vector<rvalue> wsR;
            wsR.push_back(get_omega(main0[snR], snR));
            for (int k = 0; k < N-1; k++)
                wsR.push_back(get_omega(additional0[snR][k], snR));
            wsR.push_back(get_omega(main0[snR+1], snR));
            double diff_wR = approximate(xR, snR, wsR).w1 - rdata_here.w1;

            // L or +
            std::vector<rvalue> wsL;
            wsL.push_back(get_omega(main0[snL], snL));
            for (int k = 0; k < N-1; k++)
                wsL.push_back(get_omega(additional0[snL][k], snL));
            wsL.push_back(get_omega(main0[snL+1], snL));
            double diff_wL = approximate(xL, snL, wsL).w2 - rdata_here.w2;

            rvalue diff_w = rvalue(diff_wR, diff_wL);

            value temp = get_u(diff_w, snR, snL);

            additional1[i][j].v = additional0[i][j].v + temp.v;
            additional1[i][j].s = additional0[i][j].s + temp.s;
		}

    main0.swap(main1);
    additional0.swap(additional1);

}

rvalue gcmethod_1d::get_omega(const value & u, int segnR, int segnL)
{
    rvalue temp;
    temp.w1 = (-c[segnR]*rho[segnR]*u.v + u.s) / 2.0;
    temp.w2 = (+c[segnL]*rho[segnL]*u.v + u.s) / 2.0;
    return temp;
}

rvalue gcmethod_1d::get_omega(const value & u, int segn)
{
    return get_omega(u, segn, segn);
}

value gcmethod_1d::get_u(const rvalue & w, int segn)
{
    return get_u(w, segn, segn);
}

value gcmethod_1d::get_u(const rvalue & w, int segnR, int segnL)
{
    value temp;
    temp.v = -w.w1/ (c[segnR]*rho[segnR]) + w.w2 / (c[segnL]*rho[segnL]);
    temp.s = w.w1 + w.w2;
    return temp;
}

rvalue gcmethod_1d::approximate(double x, int n_of_segment, std::vector<rvalue> & ws)
{
	switch ( N )
	{
        case 1 : return approximate_linear(x, n_of_segment, ws); break;
        case 2 : return approximate_quadratic(x, n_of_segment, ws); break;
        case 3 : return approximate_cubic(x, n_of_segment, ws); break;
        default : return rvalue(0.0, 0.0);
	}
}


rvalue gcmethod_1d::approximate_linear(double x, int n_of_segment, std::vector<rvalue> & ws)
{
    double x1 = mesh.get_point(n_of_segment);
    double x2 = mesh.get_point(n_of_segment+1);
    //assert(x>=x1 && x<=x2);
    rvalue w1 = ws[0];
    rvalue w2 = ws[1];
    rvalue temp = rvalue(((x-x1) * w2.w1 + (x2-x)*w1.w1)/h, ((x-x1) * w2.w2 + (x2-x)*w1.w2)/h);
    return temp;
}


rvalue gcmethod_1d::approximate_quadratic(double x, int n_of_segment, std::vector<rvalue> & ws)
{
    double x1 = mesh.get_point(n_of_segment);
    double x2 = mesh.get_point(n_of_segment+1);
    assert(x>=x1 && x<=x2);
    rvalue w1 = ws[0];
    rvalue w2 = ws[1];
    rvalue w3 = ws[2];
    double s = (x-x1)/h;
    rvalue w = w1 * (2*s*s-3*s+1) + w2 * 4*s*(1-s) + w3 * s*(2*s-1);
    return w;
}

rvalue gcmethod_1d::approximate_cubic(double x, int n_of_segment, std::vector<rvalue> & ws)
{
    double x1 = mesh.get_point(n_of_segment);
    double x2 = mesh.get_point(n_of_segment+1);
    assert(x>=x1 && x<=x2);
    rvalue w1 = ws[0];
    rvalue w2 = ws[1];
    rvalue w3 = ws[2];
    rvalue w4 = ws[3];
    double s = (x-x1)/h;
    rvalue w = w1 * (1-5.5*s+9*s*s-4.5*s*s*s) + w2 * s*(9-22.5*s+13.5*s*s) + w3 * s*(-4.5+18*s-13.5*s*s) + w4 * s * (1 - 4.5*s+4.5*s*s);
    return w;
}


void gcmethod_1d::calculate()
{
	for ( int i = 0; i < number_of_steps; i++)
	{
        if (i % saving_frequency == 0)
            save_to_vtk("out/out_" + std::to_string(i/saving_frequency) + ".vtk");
		step();
	}
    save_to_vtk("out/out_" + std::to_string(number_of_steps/saving_frequency) + ".vtk");
}

double gcmethod_1d::find_divergence_q(std::string type, double cl, double cr, double rl, double rr)
{
    double left = 0.00001;
    double right = 200;
    for (int k = 0; k < 20; k++)
    {
        if (type == "tau")
            tau = (left + right)/2;
        else if (type == "cl")
            cl = (left + right)/2;
        else if (type == "cr")
            cr = (left + right)/2;
        else if (type == "rl")
            rl = (left + right)/2;
        else if (type == "rr")
            rr = (left + right)/2;


        set_one_border_params(cl, cr, rl, rr, (mesh.get_left()+mesh.get_right())/2);
        set_initial();
        double max = 0;
        for (int i = 0; i < number_of_steps; i++)
        {
            step();

            max = 0;
            for (int j = 0; j < main0.size(); j++)
            {
                if (main0.at(j).v > max)
                    max = main0.at(j).v;
            }
            if (max > 20.0)
                break;
        }
        if (max > 20)
            right = (left + right)/2;
        else
            left = (left + right)/2;
    }
    return left;
}

void gcmethod_1d::calculate_diverge_contour()
{
    double cl = 1.0;
    double cr = 1.0;
    double rl = 1.0;
    double rr = 10.0;

    for ( double cl = 0.000001; cl <= 1.0; cl *= 10)
    {
        double diverge_t = find_divergence_q("tau", cl, cr, rl, rr);
        //std::cout << "cr : " << cr << ", tau : " << diverge_t << ", q : " << diverge_t/h << std::endl;
        std::cout << cl << " " << diverge_t/h << std::endl;
    }

    for ( double cl = 0.1; cl <= 22.0; cl > 0.9999 ? cl += 1.0 : cl += 0.1)
    {
        double diverge_t = find_divergence_q("tau", cl, cr, rl, rr);
        //std::cout << "cr : " << cr << ", tau : " << diverge_t << ", q : " << diverge_t/h << std::endl;
        std::cout << cl << " " << diverge_t/h << std::endl;
    }
}

void gcmethod_1d::save_to_vtk(std::string name)
{
    std::cout << name << std::endl;
    std::ofstream vtk_file;
    vtk_file.open(name.c_str(), std::ios::out);
	vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
	vtk_file << "DATASET POLYDATA\nPOINTS " << mesh.get_number_of_points() << " float\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
		vtk_file << mesh.get_point(i) << " " << 0.0 << " "  << 0.0 << "\n";
	vtk_file << "\nLINES " << mesh.get_number_of_segments() << " " << mesh.get_number_of_segments()*3 << "\n";
	for (int i = 0; i < mesh.get_number_of_segments(); i++)
		vtk_file << 2 << " " << i << " " << i+1 << "\n";
	vtk_file << "\nPOINT_DATA " << mesh.get_number_of_points() << "\n";;
    vtk_file << "SCALARS v FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for ( int i = 0; i < mesh.get_number_of_points(); i++ )
	{
        vtk_file <<  main0.at(i).v << "\n";
	}
    vtk_file << "SCALARS sigma FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for ( int i = 0; i < mesh.get_number_of_points(); i++ )
    {
        vtk_file <<  main0.at(i).s << "\n";
    }
    vtk_file << "SCALARS c FLOAT\n";
    vtk_file << "LOOKUP_TABLE default\n";
    vtk_file <<  (c.at(mesh.get_number_of_segments()-1) + c.at(0))/2 << "\n";
    for ( int i = 0; i < mesh.get_number_of_segments()-1; i++ )
    {
        vtk_file <<  (c.at(i) + c.at(i+1))/2 << "\n";
    }
    vtk_file <<  (c.at(mesh.get_number_of_segments()-1) + c.at(0))/2 << "\n";
    vtk_file.close();
}

void gcmethod_1d::read_from_file(std::string path)
{
    using std::string;
    using std::cout;
    using std::cin;
    using std::endl;
    boost::property_tree::ptree pt;
    try
    {
        boost::property_tree::read_ini(path, pt);
    }
    catch (boost::property_tree::ini_parser_error& error)
    {
        cout
            << error.message() << ": "
            << error.filename() << ", line "
            << error.line() << endl;
        cout << "Error! Press any key to close." << endl;
        std::cin.get();
        std::exit(1);
    }
    tau = pt.get<double>("Method.tau");
    h = pt.get<double>("Method.h");
    number_of_steps = pt.get<int>("Method.number_of_steps");
    N = pt.get<int>("Method.N");
    std::string is_monotonic_str = pt.get<std::string>("Method.is_monotonic");
    is_monotonic = (is_monotonic_str == "true" || is_monotonic_str == "TRUE" || is_monotonic_str == "True");
    saving_frequency = pt.get<int>("Method.saving_frequency");
}

void mesh_1d::read_from_file(std::string path)
{
    using std::string;
    using std::cout;
    using std::cin;
    using std::endl;
    boost::property_tree::ptree pt;
    try
    {
        boost::property_tree::read_ini(path, pt);
    }
    catch (boost::property_tree::ini_parser_error& error)
    {
        cout
            << error.message() << ": "
            << error.filename() << ", line "
            << error.line() << endl;
        cout << "Error! Press any key to close." << endl;
        std::cin.get();
        std::exit(1);
    }
    h = pt.get<double>("Method.h");
    left = pt.get<double>("Mesh.left");
    right = pt.get<double>("Mesh.right");
    N = pt.get<int>("Method.N");
}
