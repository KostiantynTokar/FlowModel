#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Eigen.h>
#include <boost/program_options.hpp>
#include <tqdm.hpp>
#include "DynamicSolver.h"

using namespace std;
using namespace Eigen;
using namespace Flow;
namespace po = boost::program_options;

template<typename Scalar>
struct Rectangle
{
	constexpr Rectangle(const Scalar& x_left = 0, const Scalar& y_bot = 0, const Scalar& x_right = 1, const Scalar& y_top = 1) :
		x_left{ x_left }, y_bot{ y_bot }, x_right{ x_right }, y_top{ y_top }
	{
		assert(x_left < x_right && y_bot < y_top);
	}
	Scalar x_left;
	Scalar y_bot;
	Scalar x_right;
	Scalar y_top;
};
istream& operator>>(istream& in, Rectangle<double>& r)
{
	string s;
	in >> s;
	vector<double*> r_elems{ &r.x_left, &r.y_bot, &r.x_right, &r.y_top };
	vector<pair<string, string>> seps{ {"(", ","}, {",", ")"}, {"(", ","}, {",", ")"} };
	for (auto[i, pos] = make_pair(size_t{ 0 }, size_t{ 0 }); i < 4; ++i)
	{
		const auto sep1 = s.find(seps[i].first, pos);
		const auto sep2 = s.find(seps[i].second, sep1 + 1);
		const string substr{ s.substr(sep1 + 1, sep2 - sep1 - 1) };
		*r_elems[i] = stod(substr);
		pos = sep2;
	}

	return in;
}
ostream& operator<<(ostream& out, const Rectangle<double>& r)
{
	out << "{(" << boost::lexical_cast<string>(r.x_left)
		<< "," << boost::lexical_cast<string>(r.y_bot)
		<< "),(" << boost::lexical_cast<string>(r.x_right)
		<< "," << boost::lexical_cast<string>(r.y_top) << ")}";

	return out;
}

struct scalar_field_func
{
	scalar_field_func(const string& s = "c_p") : name{ s } {}
	template<typename Point, typename Scalar, typename ExecutionPolicy>
	function<Scalar(const DynamicSolver<Point, Scalar, ExecutionPolicy>&, const Point&)>
		get_func(const DynamicSolver<Point, Scalar, ExecutionPolicy>&) const
	{
		if (name == "c_p")
		{
			return &DynamicSolver<Point, Scalar, ExecutionPolicy>::c_p;
		}
		else if (name == "der_phi_Gamma")
		{
			return &DynamicSolver<Point, Scalar, ExecutionPolicy>::der_phi_Gamma_change;
		}
		else if (name == "der_phi_gamma")
		{
			return &DynamicSolver<Point, Scalar, ExecutionPolicy>::der_phi_gamma_change;
		}
		else if(name == "der_phi_plume")
		{
			return &DynamicSolver<Point, Scalar, ExecutionPolicy>::der_phi_plume_change;
		}
		else if (name == "null")
		{
			return [](const DynamicSolver<Point, Scalar, ExecutionPolicy>&, const Point&) {return Scalar{ 0 }; };
		}
		return nullptr;
	}
	string name;
};
void validate(boost::any& v, const vector<string>& values, scalar_field_func*, int)
{
	using namespace boost::program_options;

	validators::check_first_occurrence(v);
	const string& s = validators::get_single_string(values);
	if (s == "c_p" || s == "der_phi_Gamma" || s == "der_phi_gamma" || s == "der_phi_plume" || s == "null")
	{
		v = boost::any(scalar_field_func{ s });
	}
	else
	{
		throw validation_error{ validation_error::invalid_option_value };
	}
}
istream& operator>>(istream& in, scalar_field_func& f)
{
	in >> f.name;
	return in;
}
ostream& operator<<(ostream& out, const scalar_field_func& f)
{
	return out << f.name;
}

struct policy_name
{
	policy_name(const string& s) : name{ s } {}
	string name;
};
void validate(boost::any& v, const vector<string>& values, policy_name*, int)
{
	using namespace boost::program_options;

	validators::check_first_occurrence(v);
	const string& s = validators::get_single_string(values);
	if (s == "seq" || s == "par" || s == "par_unseq")
	{
		v = boost::any(policy_name{ s });
	}
	else
	{
		throw validation_error{ validation_error::invalid_option_value };
	}
}
istream& operator>>(istream& in, policy_name& p)
{
	in >> p.name;
	return in;
}
ostream& operator<<(ostream& out, const policy_name& p)
{
	return out << p.name;
}

tuple<Obstacle<Vector2d>, vector<vector<Vector2d>>> create_obstacle(const size_t n_nodes,
	const Rectangle<double>& r);

template<typename Point, typename Scalar>
vector<Vector2d> create_grid(const Rectangle<Scalar>& r,
	const size_t rows, const size_t cols);

//ƒл€ кожного кроку повертаЇ tau, скал€рне поле та точки вихору
template<typename Point, typename Scalar, typename ExecutionPolicy, typename ExecutionPolicySF>
tuple<vector<Scalar>, vector<vector<Scalar>>, vector<vector<Point>>>
run(const ExecutionPolicySF epsf, DynamicSolver<Point, Scalar, ExecutionPolicy> ds,
	const function<Scalar(const DynamicSolver<Point, Scalar, ExecutionPolicy>&, const Point&)> scalar_field,
	const size_t n_steps,
	const vector<Point>& grid);

template<typename Scalar>
void save(const vector<Scalar>& taus, const vector<vector<Scalar>>& vals,
	const Rectangle<Scalar>& r,
	const size_t rows, const size_t cols,
	const vector<vector<Vector2d>>& parts_nodes,
	const vector<vector<Vector2d>>& plumes,
	const string fname = "output.txt", const string sep = ";");

//output, scalar_field, steps, nodes_number, rows, cols, region, obstacle_region, sep, exec_policy, exec_policy_sf
template<typename Scalar>
tuple<string, Scalar, Scalar, scalar_field_func, size_t, size_t, size_t, size_t, Rectangle<Scalar>, Rectangle<Scalar>, string, policy_name, policy_name>
input_defauls();

///////////////////////////////////////////////////////////////////////////////
////*****************************IMPLEMENTATION****************************////
///////////////////////////////////////////////////////////////////////////////

tuple<Obstacle<Vector2d>, vector<vector<Vector2d>>> create_obstacle(const size_t n_nodes,
	const Rectangle<double>& r)
{
	/*const auto n_parts_nodes = [n_nodes]()->tuple<size_t, size_t>
	{
		if (n_nodes % 2 == 0)
		{
			if ((n_nodes / 2) % 2 == 0)

			{
				return make_tuple(n_nodes / 2 + 1, n_nodes / 2 - 1);
			}
			else
			{
				return make_tuple(n_nodes / 2, n_nodes / 2);
			}
		}
		else
		{
			if ((n_nodes / 2 + 1) % 2 == 0)
			{
				return make_tuple(n_nodes / 2, n_nodes / 2 + 1);
			}
			else
			{
				return make_tuple(n_nodes / 2 + 1, n_nodes / 2);
			}
		}
	};
	const auto[n_master_nodes, n_slave_nodes] = n_parts_nodes();
	const double slave_step{ (r.y_top - r.y_bot) / n_slave_nodes };
	const double master_step{ (r.x_right - r.x_left) / (n_master_nodes - 1) };

	const double slave_x{ r.x_left + (r.x_right - r.x_left) / 2. };
	const double master_y{ r.y_top };

	vector<Vector2d> slave_nodes(n_slave_nodes);
	for (size_t i{ 0 }; i < n_slave_nodes; ++i)
	{
		slave_nodes[i] = { slave_x, r.y_bot + i * slave_step };
	}
	vector<Vector2d> master_nodes(n_master_nodes);
	for (size_t i{ 0 }; i < n_master_nodes; ++i)
	{
		master_nodes[i] = { r.x_left + i * master_step, master_y };
	}

	const Line<Vector2d> master_line{ master_nodes, {0, n_master_nodes - 1} };
	const AttachedLine<Vector2d> slave_line{ slave_nodes, {0}, n_master_nodes / 2 };
	const AttachedPart<Vector2d> slave_part{ slave_line };
	Obstacle<Vector2d> obs{ master_line, {slave_part} };

	slave_nodes.push_back(master_nodes[n_master_nodes / 2]);
	vector<vector<Vector2d>> parts_nodes{ slave_nodes, master_nodes };

	return make_tuple(obs, parts_nodes);*/

	const auto n_parts_nodes = [n_nodes]()->tuple<size_t, size_t>
	{
		if (n_nodes % 2 == 0)
		{
			if ((n_nodes / 2) % 2 == 0)

			{
				return make_tuple(n_nodes / 2 + 1, n_nodes / 2 - 1);
			}
			else
			{
				return make_tuple(n_nodes / 2, n_nodes / 2);
			}
		}
		else
		{
			if ((n_nodes / 2 + 1) % 2 == 0)
			{
				return make_tuple(n_nodes / 2, n_nodes / 2 + 1);
			}
			else
			{
				return make_tuple(n_nodes / 2 + 1, n_nodes / 2);
			}
		}
	};
	const auto[n_master_nodes, n_slaves_nodes] = n_parts_nodes();
	const size_t n_slave1_nodes{ n_slaves_nodes / 2 };
	const size_t n_slave2_nodes{ n_slaves_nodes - n_slave1_nodes };

	const double x_center{ (r.x_right + r.x_left) / 2. };
	const double y_center{ (r.y_top + r.y_bot) / 2. };

	const double slave1_step{ (r.y_top - y_center) / n_slave1_nodes };
	const double slave2_step{ (y_center - r.y_bot) / n_slave2_nodes };
	const double master_step{ (r.x_right - r.x_left) / (n_master_nodes - 1) };

	vector<Vector2d> slave1_nodes(n_slave1_nodes);
	for (size_t i{ 0 }; i < n_slave1_nodes; ++i)
	{
		slave1_nodes[i] = { x_center, r.y_bot + i * slave1_step };
	}
	vector<Vector2d> slave2_nodes(n_slave2_nodes);
	for (size_t i{ 0 }; i < n_slave2_nodes; ++i)
	{
		slave2_nodes[i] = { x_center, r.y_top - i * slave2_step };
	}
	vector<Vector2d> master_nodes(n_master_nodes);
	for (size_t i{ 0 }; i < n_master_nodes; ++i)
	{
		master_nodes[i] = { r.x_left + i * master_step, y_center };
	}

	const Line<Vector2d> master_line{ master_nodes, {0, n_master_nodes - 1} };
	const AttachedLine<Vector2d> slave1_line{ slave1_nodes, {0}, n_master_nodes / 2 };
	const AttachedLine<Vector2d> slave2_line{ slave2_nodes, {0}, n_master_nodes / 2 };
	const AttachedPart<Vector2d> slave1_part{ slave1_line };
	const AttachedPart<Vector2d> slave2_part{ slave2_line };
	Obstacle<Vector2d> obs{ master_line, {slave1_part, slave2_part} };

	slave1_nodes.push_back(master_nodes[n_master_nodes / 2]);
	slave2_nodes.push_back(master_nodes[n_master_nodes / 2]);
	vector<vector<Vector2d>> parts_nodes{ slave1_nodes, slave2_nodes, master_nodes };

	return make_tuple(obs, parts_nodes);

	/*size_t n_tmp{ n_nodes };
	std::array<size_t, 3> ns;
	ns[0] = n_tmp / 3;
	for (size_t i{ 1 }; i < 3; ++i)
	{
		ns[i] = n_tmp % 3 > 0 ? --n_tmp, ns[0] + 1 : ns[0];
	}

	const double hor_dist{ (r.x_right - r.x_left) / 2. };
	const double vert_dist{ (r.y_top - r.y_bot) / 2. };
	std::array<double, 3> steps;
	for (size_t i{ 0 }; i < 3; ++i)
	{
		const double dist{ i % 2 == 0 ? vert_dist : hor_dist };
		steps[i] = dist / (i == 2 ? ns[i] + 1 : ns[i]);
	}

	vector<Vector2d> master_nodes;

	const double y_center{ (r.y_top + r.y_bot) / 2. };
	for (size_t i{ 0 }; i < ns[0]; ++i)
	{
		master_nodes.emplace_back(r.x_left, r.y_top - i * steps[0]);
	}
	for (size_t i{ 0 }; i < ns[1]; ++i)
	{
		master_nodes.emplace_back(r.x_left + i * steps[1], y_center);
	}
	for (size_t i{ 0 }; i < ns[2]; ++i)
	{
		master_nodes.emplace_back(r.x_right, y_center - i * steps[0]);
	}

	const Line<Vector2d> master_line{ master_nodes, {0, n_nodes - 1} };
	Obstacle<Vector2d> obs{ master_line };

	vector<vector<Vector2d>> parts_nodes{ master_nodes };
	return make_tuple(obs, parts_nodes);*/

	/*const auto n_parts_nodes = [](const size_t n_nodes)->tuple<size_t, size_t>
	{
		if (n_nodes % 2 == 0)

		{
			if ((n_nodes / 2) % 2 == 0)

			{
				return make_tuple(n_nodes / 2 + 1, n_nodes / 2 - 1);
			}
			else
			{
				return make_tuple(n_nodes / 2, n_nodes / 2);
			}
		}
		else
		{
			if ((n_nodes / 2 + 1) % 2 == 0)
			{
				return make_tuple(n_nodes / 2, n_nodes / 2 + 1);
			}
			else
			{
				return make_tuple(n_nodes / 2 + 1, n_nodes / 2);
			}
		}
	};
	const auto[n_master_nodes, n_slaves_nodes] = n_parts_nodes(n_nodes);
	const auto[n_slave1_nodes, n_slave2_nodes] = n_parts_nodes(n_slaves_nodes);

	const double y_center{ (r.y_top + r.y_bot) / 2. };
	const double slave_dist{ hypot(r.x_right - r.x_left, r.y_top - y_center) };

	const double slave1_step{ slave_dist / (n_slave1_nodes) };
	const double slave2_step{ slave_dist / (n_slave2_nodes) };
	const double master_step{ (r.y_top - r.y_bot) / (n_master_nodes - 1) };

	vector<Vector2d> slave1_nodes(n_slave1_nodes);
	for (size_t i{ 0 }; i < n_slave1_nodes; ++i)
	{
		const double t{ (double)i / (double)n_slave1_nodes };
		slave1_nodes[i] = { (1 - t) * r.x_right + t * r.x_left, (1 - t) * r.y_top + t * y_center };
	}
	vector<Vector2d> slave2_nodes(n_slave2_nodes);
	for (size_t i{ 0 }; i < n_slave2_nodes; ++i)
	{
		const double t{ (double)i / (double)n_slave2_nodes };
		slave2_nodes[i] = { (1 - t) * r.x_right + t * r.x_left, (1 - t) * r.y_bot + t * y_center };
	}
	vector<Vector2d> master_nodes(n_master_nodes);
	for (size_t i{ 0 }; i < n_master_nodes; ++i)
	{
		master_nodes[i] = { r.x_left, r.y_bot + i * master_step };
	}

	const Line<Vector2d> master_line{ master_nodes, {0, n_master_nodes - 1} };
	const AttachedLine<Vector2d> slave1_line{ slave1_nodes, {0}, n_master_nodes / 2 };
	const AttachedLine<Vector2d> slave2_line{ slave2_nodes, {0}, n_master_nodes / 2 };
	const AttachedPart<Vector2d> slave1_part{ slave1_line };
	const AttachedPart<Vector2d> slave2_part{ slave2_line };
	Obstacle<Vector2d> obs{ master_line, {slave1_part, slave2_part} };

	slave1_nodes.push_back(master_nodes[n_master_nodes / 2]);
	slave2_nodes.push_back(master_nodes[n_master_nodes / 2]);
	vector<vector<Vector2d>> parts_nodes{ slave1_nodes, slave2_nodes, master_nodes };

	return make_tuple(obs, parts_nodes);*/
	
}

template<typename Point, typename Scalar>
vector<Vector2d> create_grid(const Rectangle<Scalar>& r,
	const size_t rows, const size_t cols)
{
	const Scalar dist_x{ r.x_right - r.x_left };
	const Scalar dist_y{ r.y_top - r.y_bot };
	const Scalar step_x{ dist_x / cols };
	const Scalar step_y{ dist_y / rows };

	vector<Point> grid(rows * cols);
	for (size_t i{ 0 }; i < rows; ++i)
	{
		for (size_t j{ 0 }; j < cols; ++j)
		{
			grid[i * cols + j] = { r.x_left + step_x * (j + 0.5), r.y_bot + step_y * (i + 0.5) };
		}
	}
	return grid;
}

template<typename Point, typename Scalar, typename ExecutionPolicy, typename ExecutionPolicySF>
tuple<vector<Scalar>, vector<vector<Scalar>>, vector<vector<Point>> >
run(const ExecutionPolicySF epsf, DynamicSolver<Point, Scalar, ExecutionPolicy> ds,
	const function<Scalar(const DynamicSolver<Point, Scalar, ExecutionPolicy>&, const Point&)> scalar_field,
	const size_t n_steps,
	const vector<Point>& grid)
{
	vector<Scalar> taus(n_steps);
	vector<vector<Scalar>> vals(n_steps);
	vector<vector<Point>> plumes(n_steps);
	for (size_t i : tq::trange(n_steps))
	{
		taus[i] = ds.get_last_tau();

		vals[i].resize(grid.size());
		transform(epsf,
			grid.cbegin(), grid.cend(),
			vals[i].begin(),
			[&ds, &scalar_field](const Point& p) {return scalar_field(ds, p); });

		for_each(ds.get_plumes().cbegin(), ds.get_plumes().cend(),
			[&plumes, i](const vector<Point>& plume_part)
		{
			copy(plume_part.cbegin(), plume_part.cend(), back_inserter(plumes[i]));
		});

		if (i != n_steps - 1)
		{
			ds.next_step();
		}

	}
	return make_tuple(taus, vals, plumes);
}

template<typename Scalar>
void save(const vector<Scalar>& taus, const vector<vector<Scalar>>& vals,
	const Rectangle<Scalar>& r,
	const size_t rows, const size_t cols,
	const vector<vector<Vector2d>>& parts_nodes,
	const vector<vector<Vector2d>>& plumes,
	const string fname, const string sep)
{
	ofstream ofs{ fname };
	if (!ofs)
		return;

	const auto copy_to_ostream_scalar = [&ofs, &sep](const vector<Scalar>& vals)
	{
		copy(vals.cbegin(), vals.cend(), ostream_iterator<Scalar>{ ofs, sep.c_str() });
		ofs << "\n";
	};
	const IOFormat sep_format{ StreamPrecision, DontAlignCols, sep, sep };
	const auto copy_to_ostream_point = [&ofs, &sep, &sep_format](const vector<Vector2d>& points)
	{
		transform(points.cbegin(), points.cend(), ostream_iterator<WithFormat<Vector2d>>{ofs, sep.c_str()},
			[&sep_format](const Vector2d& point)
		{
			return point.format(sep_format);
		});
		ofs << "\n";
	};

	ofs << taus.size() << sep << "\n";
	ofs << r.x_left << sep << r.y_bot << sep << r.x_right << sep << r.y_top << sep << "\n";
	ofs << rows << sep << cols << sep << "\n";
	copy_to_ostream_scalar(taus);
	for_each(vals.cbegin(), vals.cend(), copy_to_ostream_scalar);
	ofs << parts_nodes.size() << sep << "\n";
	for_each(parts_nodes.cbegin(), parts_nodes.cend(), copy_to_ostream_point);
	for_each(plumes.cbegin(), plumes.cend(), copy_to_ostream_point);

	ofs.flush();
	ofs.close();
}

template<typename Scalar>
tuple<string, Scalar, Scalar, scalar_field_func, size_t, size_t, size_t, size_t, Rectangle<Scalar>, Rectangle<Scalar>, string, policy_name, policy_name>
input_defauls()
{
#ifdef _DEBUG
	return make_tuple("output.txt", Scalar{ 0 }, Scalar{ 0 }, scalar_field_func{ "c_p" }, 5, 11, 5, 5, Rectangle<Scalar>{ 0, 0, 3, 3 }, Rectangle<Scalar>{ 1, 1, 2, 2 }, ";", policy_name{ "seq" }, policy_name{ "seq" });
#else
	return make_tuple("output.txt", Scalar{ 0 }, Scalar{ 0 }, scalar_field_func{ "c_p" }, 50, 101, 25, 50, Rectangle<Scalar>{ 0, 0, 6, 3 }, Rectangle<Scalar>{ 1, 1, 2, 2 }, ";", policy_name{ "seq" }, policy_name{ "par_unseq" });
#endif // _DEBUG
}

//void tests()
//{
//	/*vector<int> v1{ 1,2,3,4 };
//	vector<char> v2{ 'a','b','c','d' };
//	for (auto elem : boost::combine(v1, v2) | boost::adaptors::indexed())
//	{
//		cout << elem.index() << ": " << elem.value().get<0>() << elem.value().get<1>() << endl;
//	}
//	auto comb{ boost::combine(v1,v2) };
//	is_base_of_v<std::input_iterator_tag, decltype(boost::make_zip_iterator(boost::make_tuple(v1.begin(), v2.begin())))::iterator_category>;
//	std::is_base_of_v <std::input_iterator_tag, boost::range_category<vector<int>>>;
//	cout << std::is_base_of_v<std::forward_iterator_tag, decltype(v2.begin())::iterator_category> << endl;
//	cout << std::is_base_of_v<std::forward_iterator_tag, decltype(boost::begin(comb))::iterator_category> << endl;
//	std::is_base_of_v<std::input_iterator_tag, decltype((v2 | boost::adaptors::indexed()).begin())::iterator_category>*/
//
//	/*Line<Vector2d> master{ {{1.,2.}, {1.5,2.}, {2.,2.}}, {0, 2} };
//	shared_ptr<Line<Vector2d>> sp_master{ &master, [](Line<Vector2d>*) {} };
//	AttachedLine<Vector2d> slave_line{ {{1.5,1.}}, {0}, 1 };
//	shared_ptr<AttachedLine<Vector2d>> sp_line{ &slave_line, [](AttachedLine<Vector2d>*) {} };
//	AttachedPart<Vector2d> slave_part{ sp_line };
//	shared_ptr<BasePart<Vector2d>> sp_slave{ &slave_part, [](BasePart<Vector2d>*) {} };
//	Obstacle<Vector2d> obs{ sp_master, {sp_slave} };*/
//
//	Line<Vector2d> master{ {{1.,2.}, {1.5,2.}, {2.,2.}}, {0, 2} };
//	AttachedLine<Vector2d> slave_line{ {{1.5,1.}}, {0}, 1 };
//	//AttachedPart<Vector2d> slave_part{ slave_line };
//	Obstacle<Vector2d> obs{ master, { { slave_line } } };
//
//	//auto& plumes = obs.get_master().get_plumes();
//	//plumes[0].push_back({ 1.,2.05 });
//	//cout << obs.get_master().get_plumes()[0].size() << endl;
//
//	/*auto it{ obs.begin() };
//	decltype(auto) deref = *it;
//	cout << deref.get_grid().size() << endl;
//	decltype(auto) cur = it.get_cur_part_ptr();
//	cout << cur->get_master().get_grid().size() << endl;
//	decltype(auto) mas = it.get_master_part_ptr();
//	cout << mas->get_master().get_grid().size() << endl;
//	++it;
//	decltype(auto) deref2 = *it++;
//	cout << deref2.get_grid().size() << endl;
//	cout << boolalpha << (it == obs.end()) << endl;*/
//	/*vector<Line<Vector2d>> v;
//	copy(obs.cbegin(), obs.cend(), back_inserter(v));
//	cout << v[0].get_grid().size() << endl;
//	cout << v[1].get_grid().size() << endl;*/
//	/*auto point_iterable = obs.as_const_point_iterable_rep();
//	auto it1 = point_iterable.cbegin();
//	auto it2 = point_iterable.cend();
//	cout << it1->transpose() << endl;
//	++it1;
//	cout << it1->transpose() << endl;
//	cout << (++it1)->transpose() << endl;
//	cout << (++it1)->transpose() << endl;
//	cout << (++it1)->transpose() << endl;
//	++it1;
//	cout << boolalpha << (it1 == it2) << endl;*/
//	/*auto point_iterable = obs.as_const_point_iterable_rep();
//	copy(point_iterable.cbegin(), point_iterable.cend(), ostream_iterator<Vector2d>(cout, "\n\n"));*/
//	//transform(point_iterable.cbegin(), point_iterable.cend(), ostream_iterator<Matrix<double, 1, 2>>(cout, "\n"),
//	//	[](const Vector2d& x) {return x.transpose(); });
//	//StaticSolver<Vector2d, double, std::execution::parallel_policy> ss{ obs };
//
//	//StaticSolver<Vector2d> ss{ obs };
//	//cout << ss.V({ 1, 2  }) << endl;
//
//	DynamicSolver<Vector2d> ds{ obs };
//	cout << ds.V({ 1,2 }) << endl;
//	cout << endl;
//	cout << ds.c_p({ 1,1 }) << endl;
//	cout << endl;
//	ds.next_step();
//	cout << ds.V({ 1,2 }) << endl;
//	cout << endl;
//	cout << ds.c_p({ 1,1 }) << endl;
//	cout << endl;
//
//	DynamicSolver<Vector2d, double, std::execution::parallel_policy> ds_par{ obs };
//	cout << ds_par.V({ 1,2 }) << endl;
//	cout << endl;
//	cout << ds_par.c_p({ 1,1 }) << endl;
//	cout << endl;
//	ds_par.next_step();
//	cout << ds_par.V({ 1,2 }) << endl;
//	cout << endl;
//	cout << ds_par.c_p({ 1,1 }) << endl;
//	cout << endl;
//
//	DynamicSolver<Vector2d, double, std::execution::parallel_unsequenced_policy> ds_par_unseq{ obs };
//	cout << ds_par_unseq.V({ 1,2 }) << endl;
//	cout << endl;
//	cout << ds_par_unseq.c_p({ 1,1 }) << endl;
//	cout << endl;
//	ds_par_unseq.next_step();
//	cout << ds_par_unseq.V({ 1,2 }) << endl;
//	cout << endl;
//	cout << ds_par_unseq.c_p({ 1,1 }) << endl;
//	cout << endl;
//}
