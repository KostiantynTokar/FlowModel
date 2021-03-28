#pragma once

#include <vector>
#include <algorithm>
#include <execution>
#include <limits>
#include <cmath>
#include "Obstacle.h"
#include "Eigen/Dense.h"
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/combine.hpp>

namespace Flow
{
	using std::vector;
	using std::pair;
	using namespace Eigen;

	template<typename Point, typename Scalar = double, typename ExecutionPolicy = std::execution::sequenced_policy,
		std::enable_if_t<std::is_execution_policy_v<ExecutionPolicy>, int> = 0>
	class StaticSolver
	{
	public:
		using point_t = Point;
		using scalar_t = Scalar;
		using execution_policy_t = ExecutionPolicy;
		constexpr static execution_policy_t policy = {};

		StaticSolver(const Obstacle<point_t>& obs, point_t V_inf = { 1,0 }, scalar_t Gamma_0 = 0) :
			obs{ obs }, V_inf{ V_inf }, Gamma_0{ Gamma_0 }
		{
			init();
			auto Gamma_vec{ calc_Gamma() };
			index_t N = nodes.size();
			for (index_t i{ 0 }; i < N; ++i)
			{
				Gamma.push_back(Gamma_vec(i));
			}
		}
		virtual ~StaticSolver() = default;

		virtual point_t V(const point_t& p) const
		{
			point_t res{ V_inf };
			index_t N{ Gamma.size() };
			for (index_t i{ 0 }; i < N; ++i)
			{
				res += Gamma[i] * local_V(p, nodes[i]);
			}
			return res;
		}

		const Obstacle<point_t>& get_obstacle() const
		{
			return obs;
		}
		const point_t& get_V_inf() const
		{
			return V_inf;

		}
		const scalar_t& get_Gamma_0() const
		{
			return Gamma_0;
		}
		const vector<scalar_t>& get_Gamma() const
		{
			return Gamma;
		}

	protected:
		using index_t = size_t;

		Obstacle<point_t> obs;
		//Velocity at infinity
		point_t V_inf;
		//Circulation around obstacle
		scalar_t Gamma_0;

		//Pairs of nodes that are located next to each other on the obstacle
		vector<pair<point_t, point_t>> adjacent_points;
		//Vector of nodes of the obstacle
		vector<point_t> nodes;
		//Points between adjacent nodes
		vector<point_t> mid_points;
		//Normals in mid points
		vector<point_t> normals;
		//0.5 of minimum distance between adjacent nodes
		scalar_t delta;

		//Intensities in nodes
		vector<scalar_t> Gamma;

		void init()
		{
			init_nodes();
			init_adjacent_points();
			init_mid_points();
			init_normals();
			init_delta();
		}

		scalar_t R(const point_t& p, const point_t& p0) const
		{
			return std::max(delta, (p - p0).norm());
		}

		//If p0 is node and Gamma is intensity of this node,
		//than Gamma * local_V(p, p0) is velocity in p that caused by intensity in p0
		point_t local_V(const point_t& p, const point_t& p0) const
		{
			scalar_t loc_R{ R(p,p0) };
			scalar_t denom = 2 * M_PI * loc_R * loc_R;
			return point_t{ (p0[1] - p[1]) / denom, (p[0] - p0[0]) / denom };
		}

		//Translates vector of scalars, each of which corresponds to node,
		//to vector of scalars in dipole view,
		//calculating corresponding sums of input scalars, that depend on the obstacle shape
		vector<scalar_t> to_dipole(vector<scalar_t> vals) const
		{
			index_t N = vals.size();
			assert(N == nodes.size());

			vector<scalar_t> res(vals.size(), scalar_t{ 0 });
			auto end{ obs.cend() };
			index_t n{ 0 };
			auto master_part_ptr{ obs.cbegin().get_master_part_ptr() };
			index_t slave_ind{ 0 };
			vector<pair<index_t, scalar_t>> to_master;

			for (auto it{ obs.cbegin() }; it != end; ++it)
			{
				if (master_part_ptr != nullptr && slave_ind == master_part_ptr->get_slaves_ptr().size())
				{
					std::for_each(to_master.cbegin(), to_master.cend(),
						[n,&vals](const pair<index_t, scalar_t>& p)
					{
						vals[n + p.first] += p.second;
					});
					to_master.clear();
					slave_ind = 0;
					master_part_ptr = it.get_master_part_ptr();
				}

				index_t line_size = it->get_grid().size();
				std::inclusive_scan(policy,
					std::next(vals.cbegin(), n), std::next(vals.cbegin(), n + line_size),
					std::next(res.begin(), n));
				n += line_size;

				if (master_part_ptr != nullptr)
				{
					auto cur_line_ptr{ it.get_cur_attached_line() };
					index_t attach_ind{ cur_line_ptr->get_attach_ind() };
					scalar_t val_to_master{ res[n - 1] };
					to_master.push_back(std::make_pair(attach_ind, val_to_master));

					++slave_ind;
				}
			}
			return res;
		}

		//Matrix of SOLE for Gamma
		virtual Matrix<scalar_t, Dynamic, Dynamic> calc_matrix() const
		{
			index_t N{ nodes.size() };
			Matrix<double, Dynamic, Dynamic> res{ N, N };

			for (index_t i{ 0 }; i < N - 1; ++i)
			{
				for (index_t j{ 0 }; j < N; ++j)
				{
					res(i, j) = local_V(mid_points[i], nodes[j]).dot(normals[i]);
				}
			}
			for (index_t j{ 0 }; j < N; ++j)
			{
				res(N - 1, j) = scalar_t{ 1 };
			}

			return res;
		}

		//Right hand side of SOLE for Gamma
		virtual Matrix<scalar_t, Dynamic, 1> calc_rhs() const
		{
			index_t N{ nodes.size() };
			Matrix<double, Dynamic, 1> res{ N };

			for (index_t i{ 0 }; i < N - 1; ++i)
			{
				res(i) = -V_inf.dot(normals[i]);
			}
			res(N - 1) = Gamma_0;
			return res;
		}

	private:
		//Init vector of nodes of the obstacle
		virtual void init_nodes()
		{
			auto iterable{ obs.as_const_point_iterable() };
			std::copy(iterable.cbegin(), iterable.cend(), std::back_inserter(nodes));
		}

		//Find pairs of nodes that are located next to each other on the obstacle
		virtual void init_adjacent_points()
		{
			auto iterable{ obs.as_const_point_iterable_rep() };
			vector<pair<point_t, point_t>> all_adjacent_points;
			std::transform(std::next(iterable.cbegin()), iterable.cend(),
				iterable.cbegin(),
				std::back_inserter(all_adjacent_points),
				[](const point_t& next, const point_t& prev)
			{
				return std::make_pair(prev, next);
			});

			auto line_iter{ obs.cbegin() };
			decltype(line_iter->get_grid().size()) counter{ 0 };
			std::remove_copy_if(all_adjacent_points.cbegin(), all_adjacent_points.cend(),
				std::back_inserter(adjacent_points),
				[&line_iter, &counter](const pair<point_t, point_t>& p)
			{
				if (counter < line_iter->get_grid().size())
				{
					++counter;
					return false;
				}
				++line_iter;
				counter = 0;
				return true;
			});
		}

		//Calculate mid points as points between adjacent nodes
		virtual void init_mid_points()
		{
			std::transform(adjacent_points.cbegin(), adjacent_points.cend(),
				std::back_inserter(mid_points),
				[](const pair<point_t, point_t>& pp)
			{
				return (pp.first + pp.second) / 2;
			});
		}

		//Calculate normals in mid points
		virtual void init_normals()
		{
			std::transform(adjacent_points.cbegin(), adjacent_points.cend(),
				std::back_inserter(normals),
				[](const pair<point_t, point_t>& pp)->point_t
			{
				point_t tau_unnormed{ (pp.second - pp.first) };
				point_t tau{ tau_unnormed / tau_unnormed.norm() };
				return { -tau[1], tau[0] };
			});
		}

		//Calculate 0.5 of minimum distance between adjacent nodes and assign to delta
		virtual void init_delta()
		{
			auto it = std::min_element(policy,
				adjacent_points.cbegin(),
				adjacent_points.cend(),
				[](const pair<point_t, point_t>& a, const pair<point_t, point_t>& b)
			{
				return (a.first - a.second).norm() < (b.first - b.second).norm();
			});
			delta = 0.5 * (it->first - it->second).norm();
		}

		//Calculate Gamma solving SOLE
		Matrix<scalar_t, Dynamic, 1> calc_Gamma() const
		{
			Matrix<scalar_t, Dynamic, Dynamic> matr = calc_matrix();
			Matrix<scalar_t, Dynamic, 1> rhs = calc_rhs();
			return matr.colPivHouseholderQr().solve(rhs);
		}
	};
}
