#ifndef DEF_SIMPLEXMETHOD_HPP
#define DEF_SIMPLEXMETHOD_HPP

#include "__include.hpp"
#include "OnDiskMatrix.hpp"

struct InfiniteSolutionsError :public exception {
	using exception::exception;
};

struct NoSolutionError :public exception {
	using exception::exception;
};



template<typename V>
class SimplexMethod {
public:
	using value_type = V;
	using const_reference = const V&;
	using reference = V&;
	using DiskMatrixType = OnDiskMatrix<value_type>;
	using SparseMatrixType = Eigen::SparseMatrix<value_type>;
	using DenseMatrixType = Eigen::Matrix<value_type, -1, -1>;
	using SolutionType = map<int, value_type>;

	// run simplex method with original matrix
	SimplexMethod(const string &filename, const DenseMatrixType &_vec_b, const DenseMatrixType &_vec_c) {
		init_not_extended(filename,_vec_b,_vec_c);
	}

	/*
	// run simplex method with matrix with artificial variables already
	SimplexMethod(const string &extended_mat_filename, const string &extended_mat_trans_filename,
		const DenseMatrixType &_vec_b, const DenseMatrixType &_vec_c) {


	}
	*/

	// returns the map of solutions, set val to be the maximum value
	SolutionType solve(value_type &val) {
		while (!run_once());
		// generate solution map
		SolutionType sol;
		auto x_b = get_x_b_vec();
		for (auto i = 0; i < x_b.rows(); ++i) {
			sol.insert(SolutionType::value_type{base[i],x_b(i,0)});
		}

		// not the most efficient way to get z though...
		val = get_z();
		return sol;
	}

protected:
	using TripletType = Eigen::Triplet<value_type>;
	unique_ptr<DiskMatrixType> ondisk_mat;
	unique_ptr<DiskMatrixType> ondisk_trans;
	DenseMatrixType vec_b;
	DenseMatrixType vec_c;
	SparseMatrixType B_inv;
	vector<int> base;
	vector<int> non_base;


	template<typename K>
	typename vector<K>::size_type guaranteed_sequencial_find(const vector<K> &vec, const K &target) {
		using size_t = vector<K>::size_type;
		size_t pos = 0;
		for (; pos < vec.size(); ++pos) {
			if (vec[pos] == target)return pos;
		}
		throw runtime_error{ "target not found." };
	}

	template<typename K>
	typename vector<K>::size_type guaranteed_find_max(const vector<K> &vec) {
		using size_t = vector<K>::size_type;
		size_t pos = 0;
		K value = vec[pos];

		for (size_t i = 1; i < vec.size(); ++i) {
			if (vec[i] > value) {
				value = vec[i];
				pos = i;
			}
		}

		return pos;
	}

	value_type get_z() {
		return (get_c_b_vec() * B_inv * vec_b)(0, 0);
	}

	DenseMatrixType get_x_b_vec() {
		return B_inv * vec_b;
	}

	DenseMatrixType get_c_n_vec() {
		DenseMatrixType ret{ 1,non_base.size() };
		auto i = 0;
		for (auto iter = non_base.begin(); iter != non_base.end(); ++iter, ++i) {
			ret(0, i) = vec_c(0, *iter);
		}
		return ret;
	}

	DenseMatrixType get_c_b_vec() {
		DenseMatrixType ret{ 1,base.size() };
		auto i = 0;
		for (auto iter = base.begin(); iter != base.end(); ++iter, ++i) {
			ret(0, i) = vec_c(0, *iter);
		}
		return ret;
	}

	DenseMatrixType get_sigma_vec() {
		DenseMatrixType ret{ 1,non_base.size() };
		

		auto c_b{ move(get_c_b_vec()) };
		auto c_n{ move(get_c_n_vec()) };

		DenseMatrixType product_row = c_b * B_inv;
		
		auto i = 0;
		for (auto iter = non_base.begin(); iter != non_base.end(); ++iter,++i) {
			auto col = ondisk_trans->read_row(*iter);;
			col.transposeInPlace();
			ret(0, i) = (product_row * col)(0, 0);
		}

		//cout << "c_n :" << c_n << endl;
		//cout << "ret: " << ret << endl;

		return c_n - ret;
	}

	void base_alteration(int out_pos, int in_pos, const DenseMatrixType &y_k) {
		//auto out_pos = guaranteed_sequencial_find(base, out_of_base);
		//auto in_pos = guaranteed_sequencial_find(non_base, into_base);

		//auto y_k = ondisk_trans->read_row(into_base);
		//y_k = (B_inv * y_k).eval();

		value_type major_element = y_k(out_pos,0);

		vector<TripletType> elements;
		
		// add elements that are in the major column
		for (auto i = 0; i < B_inv.rows(); ++i) {
			if (i == out_pos)elements.emplace_back(TripletType{ i, out_pos, (value_type)1 / major_element });
			else elements.emplace_back(TripletType{ i, out_pos, -y_k(i,0) / major_element });
		}

		// add elements in the rest of columns
		for (auto i = 0; i < B_inv.cols(); ++i) {
			if (i != out_pos)elements.emplace_back(TripletType{ i,i,(value_type)1});
		}

		// form the E matrix
		SparseMatrixType mat_e{ B_inv.rows(),B_inv.cols() };
		mat_e.setFromTriplets(elements.begin(), elements.end());

		//DenseMatrixType e_dense{ mat_e };
		//cout << e_dense << endl;

		B_inv = (mat_e * B_inv).eval();

		// set the base vectors
		auto out = base[out_pos];
		base[out_pos] = non_base[in_pos];
		non_base[in_pos] = out;

		/*
		auto print = [](int ele) {cout << ele << " "; };
		for_each(base.begin(), base.end(), print); cout << endl;
		for_each(non_base.begin(), non_base.end(), print); cout << endl;
		*/
	}

	bool run_once() {
		// optimal condition check
		bool optimal = true;
		auto sigma_vec = get_sigma_vec();

		for (auto i = 0; i < sigma_vec.cols(); ++i) {
			if (sigma_vec(0, i) > mach_eps) {
				optimal = false;
				break;
			}
		}

		if (optimal) {
			// make sure no artificial variable is in the base
			for (auto iter = base.begin(); iter != base.end(); ++iter) {
				if (*iter >= B_inv.rows())throw NoSolutionError{ "no solution" };
			}
			return true;
		}


		//cout << sigma_vec << endl;

		// find the one that should go into base
		auto into_base = 0;
		value_type max_val = sigma_vec(0, 0);
		for (auto i = 1; i < sigma_vec.cols(); ++i) {
			if (sigma_vec(0, i) > max_val) {
				into_base = i;
				max_val = sigma_vec(0, i);
			}
		}

		// determine if there is infinite solution
		bool no_sol = true;
		auto p_k = ondisk_trans->read_row(non_base[into_base]);
		p_k.transposeInPlace();
		//cout << "p k:" << p_k.transpose() << endl;
		for (auto i = 0; i < p_k.rows(); ++i) {
			if (p_k(i,0) > mach_eps) {
				no_sol = false;
				break;
			}
		}

		if (no_sol)throw InfiniteSolutionsError{ "infinite solution" };

		// find the element that should go out of base
		DenseMatrixType vec_test{ vec_b.rows(),1 };
		auto y_k = (B_inv * p_k);
		auto x_b = get_x_b_vec();
		for (auto i = 0; i < y_k.rows(); ++i) {
			if (y_k(i,0) < mach_eps)vec_test(i,0) = numeric_limits<value_type>::max();
			else vec_test(i,0) = x_b(i,0) / y_k(i,0);
		}

		//cout << vec_test.transpose() << endl;

		auto out_of_base = 0;
		value_type min_val = vec_test(0,0);
		for (auto i = 1; i < vec_test.rows(); ++i) {
			if (vec_test(i,0) < min_val) {
				min_val = vec_test(i,0);
				out_of_base = i;
			}
		}


		// apply base alteration
		//cout << "out of base index: " << out_of_base << " into base index:" << into_base << "\n";
		base_alteration(out_of_base,into_base,y_k);

		return false;
	}
private:
	void warn_no_solution(){
		cerr << "the current problem does not have a solution.";
		throw runtime_error{ "no solution." };
	}

	void fill_row(const DenseMatrixType &old_row, DenseMatrixType &new_row, int fill_pos) {
		// copy from old row
		for (auto i = 0; i < old_row.cols(); ++i) {
			new_row(0, i) = old_row(0, i);
		}

		// fill the rest with zero
		for (auto i = old_row.cols(); i < new_row.cols(); ++i) {
			new_row(0, i) = (value_type)0.0;
		}

		// set target position one
		new_row(0, old_row.cols() + fill_pos) = (value_type)1.0;
	}

	// after init, check whether the size of matrices are correct
	void vector_size_check(DiskMatrixType &original_matrix, const DenseMatrixType &_vec_b, const DenseMatrixType &_vec_c) {
		const char* size_matching_info = "matrix size does not match";
		expr_check(original_matrix.cols() >= original_matrix.rows(), "the size of matrix is invalid");
		expr_check(original_matrix.rows() == _vec_b.rows(), size_matching_info);
		expr_check(original_matrix.cols() == _vec_c.cols(), size_matching_info);
	}

	value_type find_big_M(const DenseMatrixType &vec_c) {
		auto max = numeric_limits<value_type>::min();
		for (auto i = 0; i < vec_c.cols(); ++i) {
			if (abs(vec_c(0, i)) > max)max = abs(vec_c(0, i));
		}

		// suppose 200 is a good amplification
		return (value_type)200 * max;
	}

	void init_not_extended(const string &filename, const DenseMatrixType &_vec_b, const DenseMatrixType &_vec_c) {
		// open the matrix file
		OnDiskMatrix<value_type> original_mat{ filename };

		vector_size_check(original_mat, _vec_b, _vec_c);

		// create a new matrix
		auto new_filename = filename + string{ "_extended" };
		ondisk_mat = move(unique_ptr<DiskMatrixType>{ new DiskMatrixType{ new_filename,original_mat.rows(),original_mat.cols() + original_mat.rows() } });


		// add artificial variables
		DenseMatrixType new_row{ 1,ondisk_mat->cols() };
		for (auto i = 0; i < original_mat.rows(); ++i) {
			auto old_row = move(original_mat.read_row(i));
			fill_row(old_row, new_row, i);
			//cout << old_row << endl;
			//cout << new_row << endl;
			ondisk_mat->write_row(new_row, i);
		}

		
		// generate transpose matrix
		auto trans_filename = filename + string{ "_t" };
		ondisk_mat->generate_transpose_matrix(trans_filename);
		ondisk_trans = move(unique_ptr<DiskMatrixType>{new DiskMatrixType{ trans_filename }});

		// copy b
		vec_b = _vec_b;
		
		// init c
		vec_c = DenseMatrixType{ 1,ondisk_mat->cols() };
		for (auto i = 0; i < _vec_c.cols(); ++i) {
			vec_c(0, i) = _vec_c(0, i);
		}
		// fill big M value
		auto big_M = find_big_M(_vec_c);
		for (auto i = _vec_c.cols(); i < vec_c.cols(); ++i) {
			vec_c(0, i) = -big_M;
		}

		// checking vec_c
		//cout << fixed << setprecision(5) << vec_c << "\n";

		// set base and non-base pointers
		for (int i = 0; i < _vec_c.cols(); ++i) {
			non_base.push_back(i);
		}
		for (int i = (int)_vec_c.cols(); i < (int)vec_c.cols(); ++i) {
			base.push_back(i);
		}

		// set B_inv matrix
		vector<TripletType> elements;
		B_inv = move(SparseMatrixType{ ondisk_mat->rows() ,ondisk_mat->rows() });
		for (auto i = 0; i < ondisk_mat->rows(); ++i) {
			elements.emplace_back(TripletType(i, i, (value_type)1));
		}
		B_inv.setFromTriplets(elements.begin(), elements.end());
	}
};

#endif // !DEF_SIMPLEXMETHOD_HPP
