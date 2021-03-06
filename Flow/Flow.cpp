// Flow.cpp: определяет точку входа для консольного приложения.
//

#include <iostream>
#include <iterator>
#include "Obstacle.h"
#include "StaticSolver.h"
#include "DynamicSolver.h"
#include <Eigen/Dense.h>
//#include <boost/range/combine.hpp>
//#include <boost/range/adaptor/indexed.hpp>
//#include <boost/range/algorithm/copy.hpp>
//#include <boost/range/category.hpp>
#include <boost/program_options.hpp>
#include <tqdm.hpp>
#include <memory>
#include "util.h"

using namespace std;
using namespace Flow;
using namespace Eigen;
namespace po = boost::program_options;

template<typename ExecutionPolicy>
void exec_all(const po::variables_map& vm)
{
	const size_t N{ vm["nodes_number"].as<size_t>() };
	const string output_fname{ vm["output"].as<string>() };
	const string sep{ vm["sep"].as<string>() };
	const Rectangle<double> region{ vm["region"].as<Rectangle<double>>() };
	const Rectangle<double> obstacle_region{ vm["obstacle_region"].as<Rectangle<double>>() };
	const scalar_field_func scalar_field_name{ vm["scalar_field"].as<scalar_field_func>() };
	const size_t rows{ scalar_field_name.name == "null" ? 1 : vm["rows"].as<size_t>() };
	const size_t cols{ scalar_field_name.name == "null" ? 1 : vm["cols"].as<size_t>() };
	const size_t steps{ vm["steps"].as<size_t>() };
	const double alpha{ vm["alpha"].as<double>() };
	//const double Gamma_0{ vm["Gamma0"].as<double>() };
	const policy_name exec_policy_sf{ vm["exec_policy_sf"].as<policy_name>() };

	const auto[obs, parts_nodes] = create_obstacle(N, obstacle_region);
	const DynamicSolver<Vector2d, double, ExecutionPolicy> ds{ obs, {cos(alpha), sin(alpha)}/*, Gamma_0*/ };
	const auto sf{ scalar_field_name.get_func(ds) };

	const vector<Vector2d> grid{ create_grid<Vector2d>(region, rows, cols) };
	const auto[taus, vals, plumes] = exec_policy_sf.name == "seq" ? run(execution::seq, ds, sf, steps, grid) : 
		(exec_policy_sf.name == "par" ? run(execution::par, ds, sf, steps, grid) : 
			run(execution::par_unseq, ds, sf, steps, grid));
	save<double>(taus, vals, region, rows, cols, parts_nodes, plumes, output_fname, sep);
}

int main(int argc, char* argv[])
{
	const auto[def_output, def_alpha, def_Gamma_0, def_scalar_field, def_steps, def_nodes_number, 
		def_rows, def_cols, def_region, def_obstacle_region, def_sep, def_exec_policy, def_exec_policy_sf] = input_defauls<double>();
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("output,o", po::value<string>()->default_value(def_output), "output file name")
		("alpha,a", po::value<double>()->default_value(def_alpha), "angle of flow")
		//("Gamma0,G", po::value<double>()->default_value(def_Gamma_0), "circulation")
		("scalar_field,f", po::value<scalar_field_func>()->default_value(def_scalar_field), "scalar field function; options: c_p, der_phi_Gamma, der_phi_gamma, der_phi_plume")
		("steps,s", po::value<size_t>()->default_value(def_steps), "number of steps of algorithm")
		("nodes_number,N", po::value<size_t>()->default_value(def_nodes_number), "number of nodes on obstacle")
		("rows,r", po::value<size_t>()->default_value(def_rows), "number of rows in grid of scalar field calculating")
		("cols,c", po::value<size_t>()->default_value(def_cols), "number of columns in grid of scalar field calculating")
		("region", po::value<Rectangle<double>>()->default_value(def_region), "region of scalar field calculating")
		("obstacle_region", po::value<Rectangle<double>>()->default_value(def_obstacle_region), "region of obstacle")
		("sep", po::value<string>()->default_value(def_sep), "separator in the output format")
		("exec_policy,p", po::value<policy_name>()->default_value(def_exec_policy), "algorithm execution policy; options: seq, par, par_unseq")
		("exec_policy_sf", po::value<policy_name>()->default_value(def_exec_policy_sf), "scalar field calculation execution policy; options: seq, par, par_unseq")
	;
	po::variables_map vm;
	try
	{
		po::store(po::parse_command_line(argc, argv, desc), vm);
	}
	catch (po::validation_error e)
	{
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	po::notify(vm);

	if (vm.count("help"))
	{
		cout << desc << endl;
		return EXIT_SUCCESS;
	}
	const policy_name exec_policy{ vm["exec_policy"].as<policy_name>() };

	if (exec_policy.name == "seq")
	{
		exec_all<execution::sequenced_policy>(vm);
	}
	else if(exec_policy.name == "par")
	{
		exec_all<execution::parallel_policy>(vm);
	}
	else
	{
		exec_all<execution::parallel_unsequenced_policy>(vm);
	}
	
    return EXIT_SUCCESS;
}



