#pragma once

#include <vector>
#include <algorithm>
#include <execution>
#include "Obstacle.h"
#include "StaticSolver.h"

namespace Flow
{
	using std::vector;

	template<typename Point, typename Scalar = double, typename ExecutionPolicy = std::execution::sequenced_policy,
		std::enable_if_t<std::is_execution_policy_v<ExecutionPolicy>, int> = 0>
	class DynamicSolver : public StaticSolver<Point, Scalar, ExecutionPolicy>
	{
	public:
		using point_t = Point;
		using scalar_t = Scalar;
		using execution_policy_t = ExecutionPolicy;
		using plume_t = vector<point_t>;
		using plume_vals_t = vector<scalar_t>;
		using base_t = StaticSolver<Point, Scalar, ExecutionPolicy>;

		DynamicSolver(const Obstacle<point_t>& obs, point_t V_inf = { 1,0 }/*, scalar_t Gamma_0 = 0*/, bool movable_obs = false) :
			StaticSolver<Point, Scalar, ExecutionPolicy>{ obs, V_inf, Scalar{ 0 }/*Gamma_0*/ },
			inv_matr{ base_t::calc_matrix().inverse() },
			prev_Gamma(base_t::Gamma.size(), scalar_t{ 0 }),
			last_tau{ 0 },
			movable_obs{ movable_obs }
		{
			index_t n{ 0 };
			for (auto it{ obs.cbegin() }; it != obs.cend(); ++it)
			{
				decltype(auto) line = *it;
				for (index_t ind : line.get_angles())
				{
					angle_points_indexes.push_back(n + ind);
				}
				n += line.get_grid().size();
			}
			plumes.resize(angle_points_indexes.size());
			gamma.resize(angle_points_indexes.size());
			
			for (index_t ind : angle_points_indexes)
			{
				angle_points.push_back(base_t::nodes[ind]);
			}

			d.resize(angle_points_indexes.size(), point_t{ 0,0 });
			D.resize(base_t::mid_points.size(), point_t{ 0,0 });
			recalc_c_p_params();
			recalc_obs_rectangle();
		}

		virtual point_t V(const point_t& p) const override
		{
			point_t res{ base_t::V(p) };
			res += sum_over_plumes(point_t{ 0, 0 },
				[&p, this](const point_t& plume_p, const scalar_t& val)->point_t
			{
				return val * base_t::local_V(p, plume_p);
			});
			return res;
		}

		//Term of c_p
		scalar_t der_phi_Gamma_change(const point_t& p) const
		{
			return std::transform_reduce(base_t::policy,
				D.cbegin(), D.cend(),
				base_t::nodes.cbegin(),
				scalar_t{ 0 },
				std::plus<scalar_t>{},
				[this, &p](const point_t& D_i, const point_t& node)
			{
				return D_i.dot(base_t::local_V(p, node));
			});
		}

		//Term of c_p
		scalar_t der_phi_gamma_change(const point_t& p) const
		{
			if (plumes.empty() || plumes[0].empty())
				return scalar_t{ 0 };
			return std::transform_reduce(base_t::policy,
				d.cbegin(), d.cend(),
				plumes.cbegin(),
				scalar_t{ 0 },
				std::plus<scalar_t>{},
				[this, &p](const point_t& d_i, const plume_t& plume)
			{
				return d_i.dot(base_t::local_V(p, plume.back()));
			});
		}

		//Term of c_p
		scalar_t der_phi_plume_change(const point_t& p) const
		{
			return -sum_over_plumes(scalar_t{ 0 },
				[this, &p](const point_t& plume_p, const scalar_t& val)->scalar_t
			{
				return val * base_t::local_V(p, plume_p).dot(V(plume_p));
			});
		}

		//Pressure coefficient
		scalar_t c_p(const point_t& p) const
		{
			scalar_t der_phi{ 0 };

			der_phi += der_phi_Gamma_change(p);
			der_phi += der_phi_gamma_change(p);
			der_phi += der_phi_plume_change(p);

			scalar_t V_at_p_norm{ V(p).norm() };
			scalar_t V_inf_norm{ base_t::V_inf.norm() };
			return 1 - (V_at_p_norm * V_at_p_norm + 2 * der_phi) / (V_inf_norm * V_inf_norm);
		}

		const scalar_t& get_last_tau() const
		{
			return last_tau;
		}
		const vector<plume_t>& get_plumes() const
		{
			return plumes;
		}
		const vector<plume_vals_t>& get_gamma() const
		{
			return gamma;
		}
		
		void next_step()
		{
			vector<scalar_t> new_gamma{ calc_new_gamma() };
			vector<vector<point_t>> V_at_plumes{ calc_V_at_plumes() };
			recalc_tau(V_at_plumes);
			move_plume_points(V_at_plumes);
			for (index_t i{ 0 }; i < new_gamma.size(); ++i)
			{
				gamma[i].push_back(new_gamma[i]);
			}

			std::copy(base_t::policy, base_t::Gamma.cbegin(), base_t::Gamma.cend(), prev_Gamma.begin());
			recalc_Gamma();
			recalc_c_p_params();
			if (movable_obs) recalc_obs_rectangle();
		}

	protected:
		using index_t = typename StaticSolver<Point, Scalar, ExecutionPolicy>::index_t;

		//Recalculate intensities of nodes on this step
		virtual void recalc_Gamma()
		{
			Matrix<scalar_t, Dynamic, 1> Gamma_vec{ inv_matr * calc_rhs() };
			index_t N = base_t::nodes.size();
			for (index_t i{ 0 }; i < N; ++i)
			{
				base_t::Gamma[i] = Gamma_vec(i);
			}
		}

		//Calculate intensities of new (for this step) plume points
		vector<scalar_t> calc_new_gamma() const
		{
			vector<scalar_t> res;
			for (index_t ind : angle_points_indexes)
			{
				res.push_back(base_t::Gamma[ind]);
			}
			return res;
		}

		vector<vector<point_t>> calc_V_at_plumes()
		{
			vector<vector<point_t>> res(plumes.size());
			std::transform(base_t::policy,
				plumes.cbegin(), plumes.cend(),
				angle_points.cbegin(),
				res.begin(),
				[this](const plume_t& plume, const point_t& angle_point)
			{
				vector<point_t> V_at_plume(plume.size());
				std::transform(base_t::policy,
					plume.cbegin(), plume.cend(),
					V_at_plume.begin(),
					[this](const point_t& p) {return V(p); });
				V_at_plume.push_back(V(angle_point));
				return std::move(V_at_plume);
			});
			return res;
		}

		//Calculate time step
		void recalc_tau(const vector<vector<point_t>>& V_at_plumes)
		{
			scalar_t Delta = 2 * base_t::delta;
			auto max_op{ [](scalar_t x, scalar_t y) {return std::max(x, y); } };

			scalar_t max_V = std::transform_reduce(base_t::policy,
				V_at_plumes.cbegin(), V_at_plumes.cend(),
				std::numeric_limits<scalar_t>::min(),
				max_op,
				[this](const vector<point_t>& V_at_one_plume)
			{
				auto max_it = std::max_element(base_t::policy,
					V_at_one_plume.cbegin(), V_at_one_plume.cend(),
					[](const point_t& a, const point_t& b)
				{
					return a.norm() < b.norm();
				});
				return max_it->norm();
			});

			last_tau = Delta / max_V;
		}
		
		//Apply V to plume points
		void move_plume_points(const vector<vector<point_t>>& V_at_plumes)
		{
			for (index_t i{ 0 }; i < angle_points.size(); ++i)
			{
				plumes[i].push_back(angle_points[i]);
			}

			auto euler = [this](const point_t& plume_p, const point_t& V_at_p)
			{
				return plume_p + get_last_tau() * V_at_p;
			};

			auto impenetrable_move = [this, &euler](const point_t& plume_p, const point_t& V_at_p)
			{
				point_t res{ euler(plume_p, V_at_p) };
				return imply_impenetrability(res);
			};

			auto move_plume = [this, &impenetrable_move](const plume_t& plume, const vector<point_t>& V_at_plume)
			{
				plume_t new_plume(plume.size());
				std::transform(base_t::policy,
					plume.cbegin(), plume.cend(),
					V_at_plume.cbegin(),
					new_plume.begin(),
					impenetrable_move);
				return std::move(new_plume);
			};

			std::transform(base_t::policy,
				plumes.cbegin(), plumes.cend(),
				V_at_plumes.cbegin(),
				plumes.begin(),
				move_plume);
		}

		//Impenetrability condition for point p
		point_t imply_impenetrability(const point_t& p) const
		{
			if (p[0] < bot_left[0] || p[1] < bot_left[1] || p[0] > top_right[0] || p[1] > top_right[1])
			{
				return p;
			}
			bool repeat{ false };
			point_t res{ p };
			//катет с гипотенузой 2*delta (расстояние между двумя узлами) и другим катетом delta
			const scalar_t radius{ sqrt(3) * base_t::delta };

			for (const point_t& node : base_t::nodes)
			{
				scalar_t dist{ (node - res).norm() };
				if (dist < radius)
				{
					res = (radius + 1e-5) * (res - node) / dist + node;
					repeat = true;
					break;
				}
			}
			if (repeat)
			{
				res = imply_impenetrability(res);
			}
			return res;
		}

		//Recalculate d and D;
		//have to be called on the begining of new step before calling c_p
		void recalc_c_p_params()
		{
			vector<scalar_t> q(angle_points.size());
			if (!gamma[0].empty())
			{
				std::transform(base_t::policy,
					gamma.cbegin(), gamma.cend(),
					q.begin(),
					[this](const plume_vals_t& plume_vals)
				{
					return plume_vals.back() / get_last_tau();
				});

				for (index_t i{ 0 }; i < angle_points.size(); ++i)
				{
					d[i] = q[i] * (angle_points[i] - plumes[i].back());
				}
			}

			vector<scalar_t> q_hat(base_t::nodes.size());
			if (get_last_tau() != 0)
			{
				std::transform(base_t::policy,
					base_t::Gamma.cbegin(), base_t::Gamma.cend(),
					prev_Gamma.cbegin(),
					q_hat.begin(),
					[this](const scalar_t& G, const scalar_t& G_prev)
				{
					return (G - G_prev) / get_last_tau();
				});
			}
			for (index_t i{ 0 }; i < angle_points_indexes.size(); ++i)
			{
				q_hat[angle_points_indexes[i]] += q[i];
			}

			vector<scalar_t> Q{ base_t::to_dipole(q_hat) };
			std::transform(base_t::policy,
				base_t::adjacent_points.cbegin(), base_t::adjacent_points.cend(),
				Q.cbegin(),
				D.begin(),
				[](const std::pair<point_t, point_t>& pp, const scalar_t& Q_i)
			{
				return Q_i * (pp.second - pp.first);
			});
		}

		//Apply f to all plume points and its intensities and sum;
		//need to specify return type of f, if it is a lambda.
		//Ret f(const point_t& plume_point, const scalar_t& plume_val)
		template<typename Ret, typename Func>
		Ret sum_over_plumes(Ret init, Func f) const
		{
			return std::transform_reduce(base_t::policy,
				plumes.cbegin(), plumes.cend(),
				gamma.cbegin(),
				init,
				std::plus<Ret>{},
				[this, &init, &f](const plume_t& plume, const plume_vals_t& plume_vals)
			{
				return transform_reduce(base_t::policy,
					plume.cbegin(), plume.cend(),
					plume_vals.cbegin(),
					init,
					std::plus<Ret>{},
					f);
			});
		}

		//Calculate rectangle that cover obstacle;
		//for optimization of impenetrability
		void recalc_obs_rectangle()
		{
			const auto get_x = [](const point_t& p) {return p[0]; };
			const auto get_y = [](const point_t& p) {return p[1]; };
			const auto calc_min = [](const scalar_t& a, const scalar_t& b) {return std::min(a, b); };
			const auto calc_max = [](const scalar_t& a, const scalar_t& b) {return std::max(a, b); };
			const scalar_t x0 = std::transform_reduce(std::cbegin(base_t::nodes),
				std::cend(base_t::nodes),
				std::numeric_limits<scalar_t>::max(),
				calc_min,
				get_x);
			const scalar_t y0 = std::transform_reduce(std::cbegin(base_t::nodes),
				std::cend(base_t::nodes),
				std::numeric_limits<scalar_t>::max(),
				calc_min,
				get_y);
			const scalar_t x1 = std::transform_reduce(std::cbegin(base_t::nodes),
				std::cend(base_t::nodes),
				std::numeric_limits<scalar_t>::min(),
				calc_max,
				get_x);
			const scalar_t y1 = std::transform_reduce(std::cbegin(base_t::nodes),
				std::cend(base_t::nodes),
				std::numeric_limits<scalar_t>::min(),
				calc_max,
				get_y);

			const scalar_t Delta{ 2 * base_t::delta };
			bot_left = point_t{ x0 - Delta, y0 - Delta };
			top_right = point_t{ x1 + Delta, y1 + Delta };
		}
	private:
		vector<plume_t> plumes;
		//Intensities in plume points
		vector<plume_vals_t> gamma;
		vector<index_t> angle_points_indexes;
		//Points where plumes begin
		vector<point_t> angle_points;
		//Inverse matrix of SOLE for Gamma
		Matrix<scalar_t, Dynamic, Dynamic> inv_matr;
		//Gamma on previous step for calculating derivative of phi
		vector<scalar_t> prev_Gamma;
		//Last temporal step
		scalar_t last_tau;
		//Points for calculating gamma change influence on derivative of phi
		vector<point_t> d;
		//Points for calculating Gamma change influence on derivative of phi
		vector<point_t> D;
		bool movable_obs;
		//Bottom left point of rectangle covering obstacle
		point_t bot_left;
		//Top right point of rectangle covering obstacle
		point_t top_right;

		//Right hand side of SOLE for Gamma
		virtual Matrix<scalar_t, Dynamic, 1> calc_rhs() const override
		{
			Matrix<scalar_t, Dynamic, 1> res{ StaticSolver<Point, Scalar, ExecutionPolicy>::calc_rhs() };
			index_t N{ base_t::nodes.size() };
			for (index_t i{ 0 }; i < N; ++i)
			{
				res(i) -= sum_over_plumes(scalar_t{ 0 },
					[i, N, this](const point_t& p, const scalar_t& val)
				{
					if (i != N - 1)
					{
						return val * base_t::local_V(base_t::mid_points[i], p).dot(base_t::normals[i]);
					}
					return val;
				});
			}
			return res;
		}
	};
}
