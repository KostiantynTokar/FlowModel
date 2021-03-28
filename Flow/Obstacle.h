#pragma once

#include <vector>
#include <stack>
#include <assert.h>
#include <algorithm>

namespace Flow
{
	using std::vector;
	using std::shared_ptr;
	using std::stack;

	/*template<typename Point>
	using Plume = vector<Point>;*/

	template<typename Point>
	class Line
	{
	public:
		using point_t = Point;
		using index_t = typename vector<point_t>::size_type;
		//using plume_t = Plume<point_t>;

		/*Line(const vector<point_t>& grid, const vector<index_t>& angles, const vector<plume_t>& plumes):
			grid{ grid },
			angles{ angles },
			plumes{ plumes }
		{
			assert(!grid.empty());
			for (auto angle : angles)
				assert(angle < grid.size());
			assert(angles.size() = plumes.size());
		}*/
		Line(const vector<point_t>& grid, const vector<index_t>& angles):
			grid{ grid },
			angles{ angles }
		{
			assert(!grid.empty());
			for (auto angle : angles)
				assert(angle < grid.size());
		}
		virtual ~Line() = default;

		const vector<point_t>& get_grid() const
		{
			return grid;
		}
		const vector<index_t>& get_angles() const
		{
			return angles;
		}
		/*const vector<plume_t>& get_plumes() const
		{
			return plumes;
		}
		vector<plume_t>& get_plumes()
		{
			return plumes;
		}*/

	private:
		vector<point_t> grid;
		vector<index_t> angles;
		//vector<plume_t> plumes;
	};

	template<typename Point>
	class AttachedLine : public Line<Point>
	{
	public:
		using point_t = typename Line<Point>::point_t;
		using index_t = typename Line<point_t>::index_t;
		//using plume_t = typename Line<point_t>::plume_t;

		/*AttachedLine(const vector<point_t>& grid, const vector<index_t>& angles, const vector<plume_t>& plumes,
			index_t attach_ind):
			Line<Point>(grid, angles, plumes), attach_ind{ attach_ind }
		{

		}*/
		AttachedLine(const vector<point_t>& grid, const vector<index_t>& angles,
			index_t attach_ind) :
			Line<Point>(grid, angles), attach_ind{ attach_ind }
		{

		}

		index_t get_attach_ind() const
		{
			return attach_ind;
		}

	private:
		index_t attach_ind;
	};

	template<typename Point>
	class BasePart
	{
	public:
		using point_t = Point;
		using slave_t = BasePart<point_t>;
		using slave_ptr_t = shared_ptr<slave_t>;

		BasePart(const vector<slave_ptr_t>& slaves_ptr = {}) :
			slaves_ptr{ slaves_ptr }
		{}
		template<typename ConcreteSlaveT>
		BasePart(const vector<ConcreteSlaveT>& slaves)
		{
			transform(slaves.cbegin(), slaves.cend(), back_inserter(this->slaves_ptr),
				[](const ConcreteSlaveT& slave)
			{
				return shared_ptr<slave_t>{ new ConcreteSlaveT{ slave } };
			});
		}
		virtual ~BasePart() = default;

		const vector<slave_ptr_t>& get_slaves_ptr() const
		{
			return slaves_ptr;
		}
		vector<slave_ptr_t>& get_slaves_ptr()
		{
			return slaves_ptr;
		}
		vector<slave_t> get_slaves() const
		{
			vector<slave_t> res(slaves_ptr.size());
			transform(slaves_ptr.cbegin(), slaves_ptr.cend(), res,
				[](const slave_ptr_t& slave_ptr)
			{
				return *slave_ptr;
			});
			return res;
		}
		virtual const Line<point_t>& get_master() const = 0;
		virtual Line<point_t>& get_master() = 0;
		virtual shared_ptr<const Line<point_t>> get_master_ptr() const = 0;
		virtual shared_ptr<Line<point_t>> get_master_ptr() = 0;

	private:
		vector<slave_ptr_t> slaves_ptr;
	};

	template<typename Point, template<typename> typename LineT>
	class Part : public BasePart<Point>
	{
	public:
		using point_t = Point;
		using slave_t = typename BasePart<point_t>::slave_t;
		using slave_ptr_t = typename BasePart<point_t>::slave_ptr_t;
		using master_t = LineT<point_t>;
		using master_ptr_t = shared_ptr<master_t>;
		

		Part(const master_ptr_t& master_ptr, const vector<slave_ptr_t>& slaves_ptr = {}) :
			BasePart<point_t>{ slaves_ptr }, master_ptr{ master_ptr }
		{
			const auto right_attach = [&master_ptr](const slave_ptr_t& slave_ptr)
			{
				const auto ind{ std::dynamic_pointer_cast<const AttachedLine<point_t>>(slave_ptr->get_master_ptr())->get_attach_ind() };
				return (0 < ind && ind < master_ptr->get_grid().size() - 1);
			};
			assert(std::all_of(slaves_ptr.cbegin(), slaves_ptr.cend(), right_attach));
		}
		Part(const master_t& master, const vector<Part<point_t, AttachedLine>>& slaves = {}) :
			BasePart<point_t>{ slaves },
			master_ptr{ std::make_shared<master_t>(master) }
		{
			const auto right_attach = [&master](const Part<point_t, AttachedLine>& slave)
			{
				const auto ind{ std::dynamic_pointer_cast<const AttachedLine<point_t>>(slave.get_master_ptr())->get_attach_ind() };
				return (0 < ind && ind < master.get_grid().size() - 1);
			};
			assert(std::all_of(slaves.cbegin(), slaves.cend(), right_attach));
		}

		virtual const Line<point_t>& get_master() const override
		{
			return *master_ptr;
		}
		virtual Line<point_t>& get_master() override
		{
			return *master_ptr;
		}
		virtual shared_ptr<const Line<point_t>> get_master_ptr() const override
		{
			return master_ptr;
		}
		virtual shared_ptr<Line<point_t>> get_master_ptr() override
		{
			return master_ptr;
		}

	private:
		master_ptr_t master_ptr;
	
		template<typename PartT>
		class raw_line_iterator
		{
		public:
			using is_const = std::is_const<PartT>;

			using iterator_category = std::forward_iterator_tag;
			using value_type = std::conditional_t<is_const::value, const Line<point_t>, Line<point_t>>;
			using reference = value_type&;
			using pointer = value_type*;
			using difference_type = std::ptrdiff_t;

			using index_t = typename vector<slave_ptr_t>::size_type;
			using part_ptr_t = shared_ptr<std::conditional_t<is_const::value,
				const BasePart<point_t>, BasePart<point_t>>>;
			using bt_t = std::pair<part_ptr_t, index_t>;

			raw_line_iterator(const shared_ptr<PartT>& p = nullptr):
				backtrace{},
				ind{ 0 }
			{
				cur = p;
				if (cur) find_next();
			}

			/*template<typename PartOther>
			bool operator==(const raw_line_iterator<PartOther>& other) const
			{
				return backtrace == other.backtrace && cur == other.cur && ind == other.ind;
			}
			template<typename PartOther>
			bool operator!=(const raw_line_iterator<PartOther>& other) const
			{
				return !operator==(other);
			}*/
			friend bool operator==(const raw_line_iterator& a, const raw_line_iterator& b)
			{
				return a.backtrace == b.backtrace && a.cur == b.cur && a.ind == b.ind;
			}
			friend bool operator!=(const raw_line_iterator& a, const raw_line_iterator& b)
			{
				return !(a == b);
			}

			reference operator*()
			{
				return cur->get_master();
			}
			reference operator*() const
			{ 
				return cur->get_master();
			}
			pointer operator->()
			{
				return &cur->get_master();
			}
			pointer operator->() const
			{
				return &cur->get_master();
			}

			const raw_line_iterator& operator++()
			{
				if (*this == raw_line_iterator{})
					//end - don't move
					return *this;
				find_next();
				return *this;
			}
			raw_line_iterator operator++(int)
			{
				if (*this == raw_line_iterator{})
					//end - don't move
					return *this;
				raw_line_iterator tmp{ *this };
				find_next();
				return tmp;
			}

			/*part_ptr_t& get_cur_part_ptr()
			{
				return cur;
			}
			const part_ptr_t& get_cur_part_ptr() const
			{
				return cur;
			}
			part_ptr_t& get_master_part_ptr()
			{
				if (backtrace.empty())
				{
					return nullptr;
				}
				auto&[master, tmp] = backtrace.top();
				return master;
			}
			const part_ptr_t& get_master_part_ptr() const
			{
				if (backtrace.empty())
				{
					return nullptr;
				}
				const auto&[master, tmp] = backtrace.top();
				return master;
			}
			*/
			
			part_ptr_t get_cur_part_ptr()
			{
				return cur;
			}
			
			part_ptr_t get_cur_part_ptr() const
			{
				return cur;
			}
			
			part_ptr_t get_master_part_ptr()
			{
				if (backtrace.empty())
				{
					return nullptr;
				}
				auto&[master, tmp] = backtrace.top();
				return master;
			}
			
			part_ptr_t get_master_part_ptr() const
			{
				if (backtrace.empty())
				{
					return nullptr;
				}
				const auto&[master, tmp] = backtrace.top();
				return master;
			}

			shared_ptr<std::conditional_t<is_const::value, const AttachedLine<point_t>, AttachedLine<point_t>>>
				get_cur_attached_line()
			{
				part_ptr_t cur_part_ptr{ get_cur_part_ptr() };
				auto cur_line_base_ptr{ cur_part_ptr->get_master_ptr() };
				return std::dynamic_pointer_cast<std::conditional_t<is_const::value,
					const AttachedLine<point_t>, AttachedLine<point_t>>
					>(cur_line_base_ptr);
			}

			shared_ptr<std::conditional_t<is_const::value, const AttachedLine<point_t>, AttachedLine<point_t>>>
				get_cur_attached_line() const
			{
				part_ptr_t cur_part_ptr{ get_cur_part_ptr() };
				auto cur_line_base_ptr{ cur_part_ptr->get_master_ptr() };
				return std::dynamic_pointer_cast<std::conditional_t<is_const::value,
					const AttachedLine<point_t>, AttachedLine<point_t>>
					>(cur_line_base_ptr);
			}
			
			operator raw_line_iterator<std::add_const_t<PartT>>() const
			{
				using ci = raw_line_iterator<std::add_const_t<PartT>>;
				//const shared_ptr<std::add_const_t<PartT>>& tmp = cur;
				//const shared_ptr<const BasePart<point_t>>& const_cur = cur;
				const typename ci::part_ptr_t const_cur{ cur };
				std::stack<typename ci::bt_t> const_backtrace;
				
				return { const_cur, const_backtrace, ind };
			}

		private:
			stack<bt_t> backtrace;
			part_ptr_t cur;
			index_t ind;

			raw_line_iterator(const shared_ptr<PartT>& cur, stack<bt_t> backtrace, index_t ind) :
				backtrace{ backtrace },
				ind{ ind }
			{
				this->cur = cur;
			}

			void find_next()
			{
				index_t n = cur->get_slaves_ptr().size();
				if (ind < n)
				{
					//we need to go deeper
					backtrace.push(make_pair(cur, ind));
					cur = cur->get_slaves_ptr()[ind];
					ind = 0;
					find_next();
				}
				else if (ind == n)
				{
					//return master
					++ind;
				}
				else if (ind > n && backtrace.empty())
				{
					//end
					ind = 0;
					cur = nullptr;
				}
				else
				{
					//go outside
					std::tie(cur, ind) = backtrace.top();
					backtrace.pop();
					++ind;
					find_next();
				}
			}
		};

		template<typename FwrdLineIt, bool RepeatAttachPoints = false,
			std::enable_if_t<FwrdLineIt::is_const::value, int> = 0>
		class raw_point_iterator
		{
		public:
			using is_const = typename FwrdLineIt::is_const;

			using iterator_category = std::forward_iterator_tag;
			using value_type = std::conditional_t<is_const::value, const point_t, point_t>;
			using reference = value_type & ;
			using pointer = value_type * ;
			using difference_type = std::ptrdiff_t;

			using index_t = typename vector<slave_ptr_t>::size_type;
			
			raw_point_iterator(const FwrdLineIt& first = {}, const FwrdLineIt& last = {}) :
				cur{ first },
				last{ last },
				ind{ 0 }
			{

			}

			template<typename LineItOther>
			bool operator==(const raw_point_iterator<LineItOther, RepeatAttachPoints>& other) const
			{
				return cur == other.cur && last == other.last && ind == other.ind;
			}
			template<typename LineItOther>
			bool operator!=(const raw_point_iterator<LineItOther, RepeatAttachPoints>& other) const
			{
				return !operator==(other);
			}

			reference operator*() const
			{
				return get_cur_val();
			}
			pointer operator->() const
			{
				return &get_cur_val();
			}

			const raw_point_iterator& operator++()
			{
				if (*this == raw_point_iterator{})
					//end - don't move
					return *this;
				find_next();
				return *this;
			}
			raw_point_iterator operator++(int)
			{
				if (*this == raw_point_iterator{})
					//end - don't move
					return *this;
				raw_point_iterator tmp{ *this };
				find_next();
				return tmp;
			}

		private:
			FwrdLineIt cur;
			FwrdLineIt last;
			index_t ind;

			reference get_cur_val() const
			{
				index_t n = cur->get_grid().size();
				if (ind < n)
					return cur->get_grid()[ind];
				else
				{
					decltype(auto) master_part_ptr = cur.get_master_part_ptr();
					decltype(auto) master_line = master_part_ptr->get_master();
					shared_ptr<const AttachedLine<point_t>> cur_line_ptr{ cur.get_cur_attached_line() };
					index_t attach_ind = cur_line_ptr->get_attach_ind();
					return master_line.get_grid()[attach_ind];
				}
			}
			void find_next()
			{
				index_t n = cur->get_grid().size();
				++ind;
				if (ind > n || (ind == n && !RepeatAttachPoints) || (ind == n && is_last_line()))
				{
					ind = 0;
					++cur;
					if (cur == last)
					{
						cur = {};
						last = {};
					}
				}
			}
			bool is_last_line()
			{
				auto tmp = cur;
				++tmp;
				return tmp == last;
			}
		};

	public:
		using line_iterator = raw_line_iterator<Part>;
		using const_line_iterator = raw_line_iterator<const Part>;
		using iterator = line_iterator;
		using const_iterator = const_line_iterator;

		using const_point_iterator = raw_point_iterator<const_line_iterator>;
		using const_point_iterator_rep = raw_point_iterator<const_line_iterator, true>;

		const_line_iterator cbegin() const
		{
			return shared_ptr<const Part>{this, [](auto) {}};
		}
		const_line_iterator cend() const
		{
			return {};
		}
		line_iterator begin()
		{
			return shared_ptr<Part>{this, [](auto) {}};
		}
		line_iterator end()
		{
			return {};
		}

	private:
		template<bool Rep>
		class ConstPointIterable
		{
		public:
			using it_t = std::conditional_t<Rep, const_point_iterator_rep, const_point_iterator>;
			using const_iterator = it_t;
			ConstPointIterable(const Part& p) : p{ p }
			{

			}
			it_t cbegin() const
			{
				return { p.cbegin(), p.cend() };
			}
			it_t cend() const
			{
				return {};
			}
		private:
			const Part& p;
		};

	public:
		decltype(auto) as_const_point_iterable() const
		{
			return ConstPointIterable<false>{*this};
		}
		decltype(auto) as_const_point_iterable_rep() const
		{
			return ConstPointIterable<true>{*this};
		}
	};

	template<typename Point>
	using AttachedPart = Part<Point, AttachedLine>;

	template<typename Point>
	using Obstacle = Part<Point, Line>;

};
